#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################
use lib "../lib";
use lib "../lib/perl";

use strict;
use File::Path;
use File::Basename;
use File::Find;
use Image::ExifTool;
use Math::Matlab::Local;
use Getopt::Long;
use Data::Dumper;
use File::Spec::Functions;
use Benchmark;
use POSIX;

use Config::Adhesions qw(ParseConfig);
use Image::Stack;
use Math::Matlab::Extra;
use Emerald;
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg=s", "debug|d", "movie_debug", "config_only|only_config", 
                  "lsf|l", "single_ad_folders") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my %cfg = ParseConfig(\%opt);

###############################################################################
#Main Program
###############################################################################

my @matlab_code = &build_single_ad_commands;

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'sm_vis');
$opt{error_file} = catfile($opt{errors_folder}, 'visualization', 'error.txt');
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

if (not($opt{config_only})) {
    &FA_job::run_matlab_progam(\@matlab_code,\%opt);
}

###############################################################################
#Functions
###############################################################################

sub build_single_ad_commands {
	my $puncta_file = catfile($cfg{exp_results_folder},$cfg{adhesion_props_folder},'invado_data.csv');
	my $not_puncta_file = catfile($cfg{exp_results_folder},$cfg{adhesion_props_folder},'not_invado_data.csv');
	
	my @commands;
	if (-e $puncta_file) {
		push @commands, "make_single_puncta_frames($cfg{exp_results_folder},'puncta_file',$puncta_file)";
	} else {
		print "Couldn't find puncta file ($puncta_file)";
	}
	
	if (-e $not_puncta_file) {
		push @commands, "make_single_puncta_frames($cfg{exp_results_folder},'puncta_file',$not_puncta_file)";
	} else {
		print "Couldn't find puncta file ($not_puncta_file)";
	}
    
    return @commands;
}

###############################################################################
#Documentation
###############################################################################

=head1 NAME

collect_visualizations.pl - build the visualizations of the focal adhesion
movies

=head1 SYNOPSIS

collect_mask_set.pl -cfg FA_config

=head1 DESCRIPTION

This program builds a series of movies based on files available in the tracking
matrices folder. Each file in the tracking matrices folder which ends with
'.csv' and does not contain 'no_movie' is used to build a visualization of the
tracked focal adhesions.

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=back

Optional parameter(s):

=over 

=item * debug or d: print debuging information during program execution

=item * movie_debug: pass along the debug flag to the MATLAB visualization
program, causing only a small subset of the tracked adhesions to be visualized
in a single frame

=item * config_only: only write the MATLAB config files out, do not execute the
MATLAB program

=item * emerald: build and execute long commands throught the LSF job system

=back

=head1 EXAMPLES

collect_visualizations.pl -cfg FA_config

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 6/30/2008

=cut
