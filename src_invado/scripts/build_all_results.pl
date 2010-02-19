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
use File::Find;
use Benchmark;
use Getopt::Long;
use Cwd;
use Data::Dumper;

use Config::Adhesions qw(ParseConfig);

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l", "skip_vis|skip_visualization", 
        "only_vis|vis_only|only_visualization", "exp_filter=s") or die;

if (-e '/opt/lsf/bin/bjobs' && not $opt{lsf}) {
	die "LSF appear to be installed on this machine, don't you want to use it?" 
}	

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

die "The visualization options can't be specified without the lsf option (lsf)"
  if ($opt{skip_vis} || $opt{only_vis}) && not($opt{lsf});

die "Only one of the options skip_vis or only_vis can be selected."
  if $opt{skip_vis} && $opt{only_vis};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

################################################################################
# Main
################################################################################
$t1 = new Benchmark;
$|  = 1;

my $cfg_suffix = basename($opt{cfg});
$cfg_suffix =~ s/.*\.(.*)/$1/;

my @config_files = File::Find::Rule->file()->name( "*.$cfg_suffix" )->in( ($cfg{data_folder}) );
@config_files = sort @config_files;
if (exists($opt{exp_filter})) {
   @config_files = grep $_ =~ /$opt{exp_filter}/, @config_files;
}

my @runtime_files = map catfile(dirname($_), "run.txt"), @config_files;

#collect all the config information 
my %all_configs;
for my $file (@config_files) {
	my %temp_opt = %opt;
	$temp_opt{cfg} = $file;
	%{$all_configs{$file}} = ParseConfig(\%temp_opt);
}

#This data structure acts as the overall control mechanism for which programs
#are run to perform the analysis. The first layer of the array are another set
#of arrays with sets of commands that can be executed simultaneously. The next
#layer holds all of those commands with the appropriate directory to execute the
#commands in.
my @overall_command_seq = (
	[ [ "../find_cell_features",      "./setup_results_folder.pl" ], ],
	[ [ "../find_cell_features",      "./register_gel_image_set.pl" ], ],
	[ [ "../find_cell_features",      "./find_cascaded_registrations.pl" ], ],
	[ [ "../find_cell_features",      "./apply_registration_set.pl" ], ],
	[ [ "../find_cell_features",      "./collect_mask_image_set.pl" ], ],
	[ [ "../find_cell_features",      "./find_min_max.pl" ], ],
	[ [ "../find_cell_features",      "./find_image_thresholds.pl" ], ],
	[ [ "../find_cell_features",      "./collect_degradation_image_set.pl" ], ],
	[ [ "../find_cell_features",      "./collect_fa_image_set.pl" ], ],
	[ [ "../find_cell_features",      "./collect_fa_properties.pl" ], ],
	[ [ "../analyze_cell_features",   "./build_tracking_data.pl" ], ],
	[ [ "../analyze_cell_features",   "./track_adhesions.pl" ], ],
	[ [ "../analyze_cell_features",   "./gather_tracking_results.pl" ], ],
	[ [ "../analyze_cell_features",   "./build_R_models.pl" ], ],
	[ [ "../visualize_cell_features", "./collect_visualizations.pl" ], ],
);

if ($opt{lsf}) {
    if ($opt{skip_vis}) {
        @overall_command_seq = @overall_command_seq[ 0 .. $#overall_command_seq - 1 ];
    } elsif ($opt{only_vis}) {
        @overall_command_seq = @overall_command_seq[-1];
    }
    
    if (not($opt{debug})) {
        my $job_queue_thread = threads->create('shift_idle_jobs','');
    }
    
    my $starting_dir = getcwd;
    for (@overall_command_seq) {
        my @command_seq = @{$_};
        my @command_seq = map { [ $_->[0], $_->[1] . " -lsf" ] } @command_seq;
        print "Starting on $command_seq[0][1]\n";
        &execute_command_seq(\@command_seq, $starting_dir);

        #If debugging is not on, we want to wait till the current jobs finish
        #and then check the file complements of the experiments for completeness
        if (not($opt{debug})) {
            &wait_till_LSF_jobs_finish;
            print "Checking for all output files on command $command_seq[0][1]\n";
            my @exp_to_retry = &check_file_sets(\@config_files);

            for (1..3) {
                #if no experiments are left to retry, we break out and continue
                #to the next command set
                if (not(@exp_to_retry)) {
                    last;
                    next;
                }
                print "Retrying these experiments:\n".
                      join("\n", @exp_to_retry) . "\n";
                &execute_command_seq(\@command_seq, $starting_dir, \@exp_to_retry);
                &wait_till_LSF_jobs_finish;
                @exp_to_retry = &check_file_sets(\@config_files);
            }
            if (@exp_to_retry) {
                die "Problem with collecting full file complement on experiments after three retries:\n" .
                    join("\n", @exp_to_retry);
            } else {
                print "Output file set complete, moving on.\n\n\n";
            }
        }
    }

    if (not($opt{debug})) {
        #Find and clean up the error files produced during program execution
        our @error_files;
        &File::Find::find(\&remove_unimportant_errors, ($cfg{results_folder}));
    	
        system("bsub -J \"Job Finished: $opt{cfg}\" tail " . join(" ", @error_files));
    }
} else {
    my $starting_dir = getcwd;
    for (@overall_command_seq) {
		my @command_seq = @{$_};
		$command_seq[0][1] =~ m#/(.*)\.pl#;
		
		print "Starting on $1\n";
        
		my $command_start_bench = new Benchmark;
        
        &execute_command_seq(\@command_seq, $starting_dir);
		
		my $command_end_bench = new Benchmark;
		my $td = timediff($command_end_bench, $command_start_bench);
		print "The command took:",timestr($td),"\n\n";
	}
}

$t2 = new Benchmark;

################################################################################
# Functions
################################################################################

#######################################
# LSF
#######################################

sub wait_till_LSF_jobs_finish {
    #After each step of the pipeline, we want to wait till all the individual
    #jobs are completed, which will be checked three times
    for (1 .. 3) {
        print "On check number $_\n";
        my $sleep_time = 5;
        do {
            sleep($sleep_time);
            $sleep_time++;
        } while (&running_LSF_jobs);
    }
}

sub running_LSF_jobs {
    my $bjobs_command = "bjobs";
    if (defined $cfg{job_group}) {
        $bjobs_command .= " -g $cfg{job_group}";
    }

    my @lines = `$bjobs_command`;
    if (scalar(@lines) <= 1) {
        return 0;
    } else {
        return 1;
    }
}

sub shift_idle_jobs {
    my $bjobs_command = "bjobs";
    if (defined $cfg{job_group}) {
        $bjobs_command .= " -g $cfg{job_group}";
    }
    
    my @lines = `$bjobs_command`;
    my @running_lines = `$bjobs_command -r`;
    
    if (scalar(@lines) > 1) {
        #dealing with the strange case where all the idle jobs refuse to be
        #shifted into a running queue position, so we slowly shift all the
        #pending jobs into the week queue, where they will certainly get a
        #running jobs spot
        if (scalar(@lines) > scalar(@running_lines)) {
            &move_job_to_week_queue($lines[-1]);
        }
    }
    sleep(60);
    &shift_idle_jobs();
}

sub move_job_to_week_queue {
    my $line = pop @_;
    if ($line =~ /^(\d+)/) {
       system("bmod -q week $1");
    }
}

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
        foreach my $cfg_file (@these_config_files) {
            my $config_command = "$command -cfg $cfg_file";
			# print $all_configs{$cfg_file}{exp_data_folder}, "\n\n";
            chdir $dir;
            my $return_code = 0;
            if ($opt{debug}) {
				print "Working in directory: $dir\n";
                print $config_command, "\n";
            } else {
                print "RUNNING: $config_command\n";
                $return_code = system($config_command);
				print "RETURN CODE: $return_code\n";
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
    }
}

sub check_file_sets {
    my @config_files = @{$_[0]};

    my @exp_to_retry;
    foreach my $this_config_file (@config_files) {
        my $return_code = system "./check_file_complement.pl -cfg $this_config_file";
        
        #if the return code is anything besides zero, add that config file back
        #to the exp_to_retry list
        push @exp_to_retry, $this_config_file if ($return_code);
    }

    return @exp_to_retry;
}

sub remove_unimportant_errors {
    if ($File::Find::name =~ /error.txt/) {
        open INPUT, "$_" or die "$!";
        my @errors = <INPUT>;
        close INPUT;

        my @cleaned_errors;
        foreach my $line (@errors) {
            if ($line =~ /Pending job threshold reached./) {
                next;
            }
            push @cleaned_errors, $line;
        }

        unlink $_;

        open OUTPUT, ">$_";
        print OUTPUT @cleaned_errors;
        close OUTPUT;
    }
}
