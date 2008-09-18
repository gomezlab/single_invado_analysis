#!/usr/bin/perl

################################################################################
# Global Variables and Modules
################################################################################
use lib "../lib";
use lib "../lib/perl";

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
GetOptions(\%opt, "cfg|c=s", "debug|d", "emerald|e", "skip_vis|skip_visualization", "only_vis|only_visualization")
  or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

die "The visualization options (skip_vis or only_vis) can't be specified without the emerald option (emerald)"
  if ($opt{skip_vis} || $opt{only_vis}) && not($opt{emerald});

die "Only one of the options skip_vis or only_vis can be selected."
  if $opt{skip_vis} && $opt{only_vis};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

################################################################################
# Main
################################################################################
$t1 = new Benchmark;
$|  = 1;

my @config_files = sort <$cfg{data_folder}/*/*.cfg>;
@config_files = ($config_files[0]) if $opt{debug};
my @runtime_files = map catfile(dirname($_), "run.txt"), @config_files;

if ($opt{emerald}) {
    my $starting_dir        = getcwd;
    my @overall_command_seq = (
        [ [ "../find_cell_features",      "./setup_results_folder.pl" ], ],
        [ [ "../find_cell_features",      "./collect_mask_image_set.pl" ], ],
        [ [ "../find_cell_features",      "./collect_fa_image_set.pl" ], ],
        [ [ "../find_cell_features",      "./collect_fa_properties.pl" ], ],
        [ [ "../analyze_cell_features",   "./build_tracking_data.pl" ], ],
        [ [ "../analyze_cell_features",   "./track_adhesions.pl" ], ],
        [ [ "../analyze_cell_features",   "./gather_tracking_results.pl" ], ],
        [ [ "../analyze_cell_features",   "./filter_tracking_matrix.pl" ], ],
        [ [ "../visualize_cell_features", "./collect_visualizations.pl" ], ],
    );
    if ($opt{skip_vis}) {
        @overall_command_seq = @overall_command_seq[ 0 .. $#overall_command_seq - 1 ];
    } elsif ($opt{only_vis}) {
        @overall_command_seq = @overall_command_seq[-1];
    }

    for (@overall_command_seq) {
        my @command_seq = @{$_};
        my @command_seq = map { [ $_->[0], $_->[1] . " -emerald" ] } @command_seq;
        &execute_command_seq(\@command_seq, $starting_dir);
        &wait_till_LSF_jobs_finish;
    }
    find(\&remove_unimportant_errors, ($cfg{results_folder}));
} else {
    unlink(<$cfg{data_folder}/time_series_*/stat*>);

    my $max_processes = 4;

    my @processes : shared;

    @processes =
      map { "nice -20 ./build_results.pl -cfg $config_files[$_] -d > $runtime_files[$_]" } (0 .. $#config_files);

    my @started;
    while (@processes) {
        while (@processes && &gather_running_status(@started) < $max_processes) {
            my $command = shift @processes;
            my $config  = shift @config_files;
            push @started, $config;

            print "Executing $command\n" if $opt{debug};
            threads->create('execute_process', $command)->detach;
        }
    }

    while (&gather_running_status(@started)) {
        sleep 100;
    }
}

$t2 = new Benchmark;
print "Runtime: ", timestr(timediff($t2, $t1)), "\n";

open OUTPUT, ">most_recent_runtime.txt";
print OUTPUT "Runtime: ", timestr(timediff($t2, $t1)), "\n";
close OUTPUT;

################################################################################
# Functions
################################################################################

#######################################
# LSF
#######################################

sub wait_till_LSF_jobs_finish {
    for (1 .. 6) {
        my $sleep_time = 10;
        do {
            sleep($sleep_time);
            $sleep_time++;
        } while (&running_LSF_jobs);
    }
}

sub running_LSF_jobs {
    my @lines = `bjobs`;
    if (scalar(@lines) > 1) {
        return scalar(@lines) - 1;
    } else {
        0;
    }
}

sub execute_command_seq {
    my @command_seq  = @{ $_[0] };
    my $starting_dir = $_[1];
    foreach my $set (@command_seq) {
        my $dir     = $set->[0];
        my $command = $set->[1];
        foreach (@config_files) {
            my $config_command = "$command -cfg $_";
            chdir $dir;
            system $config_command;
            chdir $starting_dir;
        }
    }
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

#######################################
# non-LSF
#######################################

sub print_sys_line {
    my $process_string = shift;
    sleep rand(10) + 5;
    print "From Func: $process_string\n";
}

sub execute_process {
    my $process_string = shift;
    system($process_string);
    print "Finished $process_string\n";
}

sub gather_running_status {
    my @started_configs = @_;

    my $running = 0;
    foreach (@started_configs) {
        my $status_file = catfile(dirname($_), 'status.txt');
        if (-e $status_file) {
            open STATUS, $status_file;
            my $line = <STATUS>;
            close STATUS;

            chomp($line);
            if (not($line =~ /DONE/)) {
                $running++;
            }
        } else {
            $running++;
        }
    }
    return $running;
}
