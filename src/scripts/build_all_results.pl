#!/usr/bin/perl

################################################################################
# Global Variables and Modules
################################################################################
use lib "../lib";
use lib "../lib/perl";

use strict;
use threads;
use threads::shared;
use File::Spec::Functions;
use File::Basename;
use File::Find::Rule;
use Benchmark;
use Getopt::Long;
use Cwd;
use Data::Dumper;
use POSIX;

use Config::ImageSet qw(ParseConfig);

my %opt;
$opt{debug} = 0;
$opt{lsf} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l", "exp_filter=s", "no_email",
	"sync=s") or die;

my $lsf_return = system("which bsub > /dev/null 2> /dev/null");

#Remember a return code of 0 means success
if ($lsf_return == 0 && not $opt{lsf}) {
	die "LSF appears to be installed on this machine, don't you want to use it?" 
}	

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

my $t1 = new Benchmark;
$|  = 1;

################################################################################
# Main
################################################################################

#This data structure acts as the overall control mechanism for which programs
#are run to perform the analysis. The first layer of the array are another set
#of arrays with sets of commands that can be executed simultaneously. The next
#layer holds all of those commands with the appropriate directory to execute the
#commands in.
my @overall_command_seq = (
	[ "../find_cell_features",    "./run_matlab_over_field.pl -script flat_field_correct_gel" ],
	[ "../find_cell_features",    "./run_matlab_over_field.pl -script find_cell_mask" ],
	[ "../find_cell_features",    "./run_matlab_over_field.pl -script apply_bleaching_correction" ],
	[ "../find_cell_features",    "./run_matlab_over_field.pl -script register_with_matlab", 1],
	[ "../find_cell_features",    "./run_matlab_over_field.pl -script find_cell_mask" ],
	[ "../find_cell_features",    "./run_matlab_over_field.pl -script find_image_set_min_max" ],
	[ "../find_cell_features",    "./run_matlab_over_field.pl -script find_puncta_thresh" ],
	[ "../find_cell_features",    "./run_matlab_over_field.pl -script find_puncta", 2],
	[ "../find_cell_features",    "./run_matlab_over_field.pl -script find_puncta_properties", 2],
	[ "../find_cell_features",    "./run_matlab_over_field.pl -script ../analyze_cell_features/track_puncta" ],
	[ "../find_cell_features",    "./run_matlab_over_field.pl -script ../analyze_cell_features/gather_tracking_results" ],
	[ "../find_cell_features",    "./run_matlab_over_field.pl -script ../analyze_cell_features/find_pre_birth_diffs" ],
	[ "../analyze_cell_features", "./find_invadopodia_puncta.pl", 1],
	[ "../find_cell_features",    "./run_matlab_over_field.pl -script ../visualize_cell_features/build_all_visualizations", 1],
);

my $cfg_suffix = basename($opt{cfg});
$cfg_suffix =~ s/.*\.(.*)/$1/;

#config file processing
my @config_files = File::Find::Rule->file()->name( "*.$cfg_suffix" )->in( ($cfg{data_folder}) );
@config_files = sort @config_files;
if (exists($opt{exp_filter})) {
   @config_files = grep $_ =~ /$opt{exp_filter}/, @config_files;
}
die "No config files left after filtering." if (scalar(@config_files) == 0);

#######################################
# Program Running
#######################################
my $starting_dir = getcwd;
for (@overall_command_seq) {
	my $command_start_bench = new Benchmark;
	my @command_seq = @{$_};
	$command_seq[1] .= " -lsf" if $opt{lsf};

	print "Starting on $command_seq[1]\n";
	&execute_command_seq(\@command_seq, $starting_dir);

	#If debugging is not on, we want to wait till the current jobs finish
	#and then check the file complements of the experiments for completeness
	&wait_till_LSF_jobs_finish if ($opt{lsf} && not($opt{debug}));
	
	if (not $opt{debug}) {
		my $command_end_bench = new Benchmark;
		my $td = timediff($command_end_bench, $command_start_bench);
		print "The command took:",timestr($td),"\n";
	}
	
	#we are done with the current command, so we print a few spacer lines to
	#indicate we have moved onto a new command
	print "\n\n";
}

if (not($opt{debug})) {
	my $t_bsub = new Benchmark;
	my $td = timediff($t_bsub, $t1);
	
	my $time_diff_str = "\"Took:" . timestr($td) . "\"";
	
	my $command;
	if ($opt{sync}) {
		$command = "results/sync_to.pl -server $opt{sync}";
		system($command);
	}
	
	$command = "bsub -J \"Job Finished: $opt{cfg}\" echo $time_diff_str";
	if (not $opt{no_email}) {
		system($command) if $opt{lsf};
	}
}

my $t2 = new Benchmark;
my $td = timediff($t2, $t1);
print "\nThe pipeline took:",timestr($td),"\n\n";

################################################################################
# Functions
################################################################################

sub execute_command_seq {
    my @command_seq  = @{ $_[0] };
    my $starting_dir = $_[1];
    my @these_config_files = @config_files;
    if (scalar(@_) > 2) {
        @these_config_files = @{$_[2]};
    }
	
	my $dir     = $command_seq[0];
	my $command = $command_seq[1];
	my $cmd_bundling = 5;
	if (scalar(@command_seq) > 2) {
		$cmd_bundling = $command_seq[2];
	}

	chdir $dir;
	while (@these_config_files) {
		my $cfg_set = pop @these_config_files;
		for (2..$cmd_bundling) {
			next if (not(@these_config_files));
			my $cfg_file = pop @these_config_files;
			$cfg_set .= ",$cfg_file";
		}

		my $config_command = "$command -cfg $cfg_set";
		my $return_code = 0;
		if ($opt{debug}) {
			print $config_command, "\n";
		} else {
			$return_code = system($config_command);
		}
	}
	chdir $starting_dir;
	print "\n";
}

#######################################
# LSF
#######################################

sub wait_till_LSF_jobs_finish {
    #After each step of the pipeline, we want to wait till all the individual
    #jobs are completed, which will be checked three times
	my $total_checks = 1;
	my $first_check = &running_LSF_jobs;
	
	if (! $first_check) {
		print "No LSF jobs detected, jumping to next command.\n";
		return 0;
	} else {
		for (1 .. $total_checks) {
			print "LSF finished check number $_/$total_checks\n";
			my $sleep_time = 10;
			do {
				sleep($sleep_time);
			} while (&running_LSF_jobs);
		}
	}
	print "\n";
	return 1;
}

sub running_LSF_jobs {
    my $bjobs_command = "bjobs";
    if (defined $cfg{job_group}) {
        $bjobs_command .= " -g $cfg{job_group} 2>/dev/null";
    }

    my @lines = `$bjobs_command`;
    
    if (scalar(@lines) <= 1) {
        return 0;
    } else {
        return 1;
    }
}
