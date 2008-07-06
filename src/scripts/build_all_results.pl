#!/usr/bin/perl

################################################################################
# Global Variables and Modules
################################################################################
use threads;
use threads::shared;
use File::Spec::Functions;
use File::Basename;
use Benchmark;
use Getopt::Long;
use Cwd;
use Data::Dumper;

use lib "../lib";
use Config::Adhesions qw(ParseConfig);

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "emerald");
die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

################################################################################
# Main
################################################################################
$t1 = new Benchmark;
$| = 1;

unlink(<$cfg{data_folder}/time_series_*/stat*>);

my $max_processes = 4;

my @config_files = <$cfg{data_folder}/time_series_*/*.cfg>;

my @data_folders = map catfile(dirname($_),"run.txt"), @config_files;
my @time_series_list = map { /(time_series_\d*)/; $1; } @config_files;

if ($opt{emerald}) {
    my $starting_dir = getcwd;
    my @command_seq = (["../find_cell_features", "./setup_results_folder.pl -emerald"],
                       ["../find_cell_features", "./collect_mask_set.pl -emerald"],
                       ["../find_cell_features", "./collect_fa_image_set.pl -emerald"],
                      );
    foreach my $set (@command_seq) {
        my $dir = $set->[0];
        my $command = $set->[1];
        foreach (@config_files) {
            $command .= " -cfg $_";
            chdir $dir; 
            system $command; 
            chdir $starting_dir;
        }
    }
} else {
    my @processes :shared;

    @processes = map { "nice -20 ./build_results.pl -cfg $config_files[$_] -d > $data_folders[$_]" } (0..$#config_files);

    my @started;
    while (@processes) {
	    while (@processes && &gather_running_status(@started) < $max_processes) {
		    my $command = shift @processes;
    		my $config = shift @config_files;
	    	push @started, $config;

		    print "Executing $command\n" if $opt{debug};
    		threads->create('execute_process',$command)->detach;
	    }
    }

    while(&gather_running_status(@started)) {
	    sleep 100;
    }
}

$t2 = new Benchmark;
print "Runtime: ",timestr(timediff($t2,$t1)), "\n";

################################################################################
# Functions
################################################################################

sub print_sys_line {
	my $process_string = shift;
	sleep rand(10)+5;
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
		my $status_file = catfile(dirname($_),'status.txt');
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
