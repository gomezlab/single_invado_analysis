#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use threads;
use threads::shared;
use File::Basename;
use File::Find::Rule;
use File::Spec::Functions;
use Cwd 'abs_path';
use Getopt::Long;
use Benchmark qw(:all);

use lib "../lib";
use Config::Adhesions qw(ParseConfig);

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my $t_start = Benchmark->new();

my %opt;
$opt{debug} = 0;
$opt{extra} = "";
$opt{thread_count} = 1;
GetOptions(\%opt, "cfg|config=s", "debug|d", "program|p=s", "extra|e=s", 
                  "run_all_debug", "exp_filter=s", "filter_negate", "thread_count=s") or die;

die "Can't find cfg file specified on the command line" if not(exists $opt{cfg});
die "Can't find program to execute on the command line" if not(exists $opt{program});

print "Collecting Configuration\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

###############################################################################
# Main Program
###############################################################################
chdir(dirname($opt{program}));

my $program_base = basename($opt{program});

my $debug_string = ($opt{debug}) ? "-d" : "";

my $cfg_suffix = basename($opt{cfg});
$cfg_suffix =~ s/.*\.(.*)/$1/;

my @config_files = File::Find::Rule->file()->name( "*.$cfg_suffix" )->in( ($cfg{data_folder}) );
@config_files = sort @config_files;
if (exists($opt{exp_filter})) {
	if ($opt{filter_negate}) {
		@config_files = grep not($_ =~ /$opt{exp_filter}/), @config_files;
	} else {
		@config_files = grep $_ =~ /$opt{exp_filter}/, @config_files;
	}
} 

my @thread_list;
while (@config_files) {
    if ($opt{debug}) {
		my $this_config = shift @config_files;
        print("./$program_base -cfg $this_config $debug_string $opt{extra}\n");
    } else {
		my @running = threads->list(threads::running);
		if (scalar(@running) >= $opt{thread_count}) {
			next;
		} else {
			my $this_config = shift @config_files;
			my $new_thread = threads->create('run_command',$program_base, $this_config, $debug_string, $opt{extra});
			push @thread_list, $new_thread;
		}
    }
}

my @running = threads->list(threads::running);
while (@running) {
	@running = threads->list(threads::running);
	print "Found running " . scalar(@running) . "\n"; sleep(1);
}

my $join_count = 0;

foreach (@thread_list) {
	$_->detach();
	$join_count++;
}

my $t_end = Benchmark->new;
my $td = timediff($t_end, $t_start);
print "Full processing time was " . timestr($td, 'noc') . "\n";

###############################################################################
# Functions
###############################################################################

sub run_command {
	my $program_base = $_[0];
	my $this_config = $_[1];
	my $debug_string = $_[2];
	my $extra_commands = $_[3];

	# print "Starting on $this_config\n";
	my $t0 = Benchmark->new;
	system("./$program_base -cfg $this_config $debug_string $opt{extra}");
	my $t1 = Benchmark->new;
	my  $td = timediff($t1, $t0);
	print "Done with $this_config, took " . timestr($td, 'noc') . "\n";

	return 0;
}

###############################################################################
# Documentation
###############################################################################

=head1 NAME

run_perl_program_for_all_threads.pl

=head1 SYNOPSIS

run_program_for_all_threads.pl -cfg config.cfg -p program 

=head1 Description

This program provides the means to run a program for each of the experiments in
a set. This is helpful when a single step of a processing pipeline needs to be
executed again. Config files (those with the same suffix as the config file
provided on the command line) are searched for in subdirectories of the
$cfg{data_folder} directory.

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=item * program or p: the path to a program that should be run using the
standard "-cfg exp_config.cfg" parameter

=back

Optional parameter(s):

=over 

=item * debug or d: if present do not run the program specified in the program
parameter, instead print the commands that would be run

=item * thread_count: number of threads to use when for processing the jobs,
defaults to one

=back

=head1 EXAMPLES

run_program_for_all_threads.pl -cfg FA_config -p ../analyze_cell_features/gather_tracking_results.pl

OR

run_program_for_all_threads.pl -cfg FA_config -p
../analyze_cell_features/gather_tracking_results.pl -e '-skip_pix_props'

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 6/9/2008 

=cut
