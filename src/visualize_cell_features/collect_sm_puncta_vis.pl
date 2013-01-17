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
use Math::Matlab::Local;
use Getopt::Long;
use Data::Dumper;
use File::Spec::Functions;
use Benchmark;
use POSIX;

use Config::Adhesions qw(ParseConfig);
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

# $opt{memsize} = 1024*10;
$opt{resource} = "mem96"; 
$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'sm_vis');
$opt{error_file} = catfile($opt{error_folder}, 'error.txt');
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
		push @commands, "make_single_puncta_frames('$cfg{exp_results_folder}','puncta_file','$puncta_file')";
	} else {
		print "Couldn't find puncta file ($puncta_file)";
	}
	
	if (-e $not_puncta_file) {
		push @commands, "make_single_puncta_frames('$cfg{exp_results_folder}','puncta_file','$not_puncta_file')";
	} else {
		print "Couldn't find puncta file ($not_puncta_file)";
	}
    
    return @commands;
}

###############################################################################
#Documentation
###############################################################################
