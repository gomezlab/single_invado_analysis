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
use File::Copy;
use Benchmark;
use Getopt::Long;
use Cwd;
use Data::Dumper;

use Config::Adhesions qw(ParseConfig);

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d") or die;

my $base_data_dir = '/Volumes/Data/projects/fa_webapp/data/';
my $base_results_dir = '/Volumes/Data/projects/fa_webapp/results/';
my $run_file = '/tmp/run/FA_processing';

if (-e $run_file) {
    die "Processing occuring, exiting" if $opt{debug};
    exit;
} else {
    &create_run_file;
}

################################################################################
# Main
################################################################################
my @results_nums = <$base_results_dir*>;
@results_nums = map basename($_), @results_nums;

my @dirs = <$base_data_dir*>;
@dirs = grep {
    if (/config/) {
        0;
    } else {
        1;
    }
} @dirs;
if (scalar(@dirs) == 0) {
    &remove_run_file;
    die "No Files to work on" if $opt{debug};
    exit;
}

my %current_dir_ages;
foreach my $this_dir (@dirs) {
    if (-d $this_dir) {
        if (-e catfile($this_dir, "working")) {
            next;
        }
        if (grep $this_dir =~ /$_/, @results_nums) {
            next;
        }
        my($dev, $ino, $mode, $nlink, $uid, $gid, $rdev ,
           $size, $atime, $mtime, $ctime, $blksize, $blocks) = stat($this_dir);
        $current_dir_ages{$this_dir} = $mtime;
    }
}
if (scalar keys %current_dir_ages == 0) {
    &remove_run_file;
    die "No finished data folders to work on" if $opt{debug};
    exit;
}

my @sorted_dir_names = sort { $current_dir_ages{$a} <=> $current_dir_ages{$b} } keys %current_dir_ages;

my $oldest_dir_num = basename($sorted_dir_names[0]);

my $command = "./build_all_results.pl -cfg ../../data/config/webapp_default.cfg -exp_filter $oldest_dir_num\n";
if ($opt{debug}) {
    print $command;
} else {
    system $command;
}

chdir('../../results');


&remove_run_file;
################################################################################
# Functions
################################################################################

sub remove_run_file {
    unlink($run_file) or die "$!";
}

sub create_run_file {
    open OUTPUT, ">$run_file" or die "$!";
    #the $$ variable is the PID of the currently running program
    print OUTPUT $$;
    close OUTPUT;
}
