#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Basename;
use File::Spec::Functions;
use Cwd 'abs_path';
use Getopt::Long;

use lib "../lib";
use Config::Adhesions qw(ParseConfig);

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
$opt{extra} = "";
GetOptions(\%opt, "cfg|config=s", "debug|d", "program|p=s", "extra|e=s", "run_all_debug") or die;

die "Can't find cfg file ($opt{cfg}) specified on the command line" if not(exists $opt{cfg});
die "Can't find program to execute on the command line" if not(exists $opt{program});

print "Collecting Configuration\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

###############################################################################
#Main Program
###############################################################################
chdir(dirname($opt{program}));

my $program_base = basename($opt{program});

my $debug_string = ($opt{debug}) ? "-d" : "";

my $cfg_suffix = basename($opt{cfg});
$cfg_suffix =~ s/.*\.(.*)/$1/;

my @exp = <$cfg{data_folder}/*/*$cfg_suffix>;

foreach (@exp) {
    next if /config\/default/; 
    if ($opt{run_all_debug}) {
        print("./$program_base -cfg $_ $debug_string $opt{extra}\n");
    } else {
        print("$_\n");
        system("./$program_base -cfg $_ $debug_string $opt{extra}");
    }
}

###############################################################################
#Documentation
###############################################################################

=head1 NAME

run_perl_program_for_all.pl

=head1 SYNOPSIS

run_program_for_all.pl -cfg config.cfg -p program 

=head1 Description

This program provides the means to run a program for each of the experiments in a set. This is helpful when a single step of a processing pipeline needs to be executed again. Config files (those with the same suffix as the config file provided on the command line) are searched for in subdirectories of the $cfg{data_folder} directory.

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=item * program or p: the path to a program that should be run using the standard "-cfg exp_config.cfg" parameter

=back

Optional parameter(s):

=over 

=item * run_all_debug: if present do not run the program specified in the program parameter, instead print the commands that would be run

=item * debug or d: if present print debuging information during program execution, including debug information for all subsequent program calls

=back

=head1 EXAMPLES

run_program_for_all.pl -cfg FA_config -p ../analyze_cell_features/gather_tracking_results.pl

OR

run_program_for_all.pl -cfg FA_config -p ../analyze_cell_features/gather_tracking_results.pl -e '-skip_pix_props'

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 6/9/2008 

=cut
