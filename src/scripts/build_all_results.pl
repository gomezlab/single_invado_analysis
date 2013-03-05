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

use Config::Adhesions qw(ParseConfig);

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
	[ [ "../find_cell_features",      "./run_matlab_over_field.pl -script find_cell_mask" ], ],
	[ [ "../find_cell_features",      "./run_matlab_over_field.pl -script apply_bleaching_correction" ], ],
	[ [ "../find_cell_features",      "./run_matlab_over_field.pl -script register_with_matlab" ], ],
	[ [ "../find_cell_features",      "./run_matlab_over_field.pl -script find_cell_mask" ], ],
	[ [ "../find_cell_features",      "./run_matlab_over_field.pl -script find_image_set_min_max" ], ],
	[ [ "../find_cell_features",      "./run_matlab_over_field.pl -script find_puncta_thresh" ], ],
	[ [ "../find_cell_features",      "./run_matlab_over_field.pl -script find_puncta" ], ],
	[ [ "../find_cell_features",      "./run_matlab_over_field.pl -script find_puncta_properties" ], ],
	[ [ "../find_cell_features",      "./run_matlab_over_field.pl -script ../analyze_cell_features/track_invado" ], ],
	[ [ "../find_cell_features",      "./run_matlab_over_field.pl -script ../analyze_cell_features/find_pre_birth_diffs" ], ],
	[ [ "../analyze_cell_features",   "./gather_tracking_results.pl" ], ],
	[ [ "../analyze_cell_features",   "./find_invadopodia_puncta.pl" ], ],
	[ [ "../visualize_cell_features", "./collect_sm_puncta_vis.pl" ], ],
	[ [ "../find_cell_features",      "./run_matlab_over_field.pl -script ../visualize_cell_features/make_tracking_visualization" ], ],
	[ [ "../find_cell_features",      "./run_matlab_over_field.pl -script ../visualize_cell_features/build_dual_highlight_image" ], ],
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
	@command_seq = map { [ $_->[0], $_->[1] . " -lsf" ] } @command_seq if $opt{lsf};
	print "Starting on $command_seq[0][1]\n";
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
	
    foreach my $set (@command_seq) {
        my $dir     = $set->[0];
        my $command = $set->[1];
        my $executed_scripts_count = 0;
        foreach my $cfg_file (@these_config_files) {
            my $config_command = "$command -cfg $cfg_file";
            chdir $dir;
            my $return_code = 0;
			$executed_scripts_count++;
			print "Done submitting: " if $executed_scripts_count == 1;
            
			if ($opt{debug}) {
                print $config_command, "\n";
            } else {
                $return_code = system($config_command);
				if ($return_code) {
                	print "PROBLEM WITH: $config_command\n";
					print "RETURN CODE: $return_code\n";
				}
				if ($executed_scripts_count % ceil(scalar(@these_config_files)/10) == 0) {
					print sprintf('%.0f%% ',100*($executed_scripts_count/scalar(@these_config_files)));
				}
            }
            chdir $starting_dir;

            #if the return code was any number beside zero, indicating a problem
            #with the program exit, remove that config file from the run and
            #continue
            if ($return_code) {
				print "REMOVING: $cfg_file\n";
                @config_files = grep $cfg_file ne $_, @config_files;
            }
        }
		print "\n";
    }
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
