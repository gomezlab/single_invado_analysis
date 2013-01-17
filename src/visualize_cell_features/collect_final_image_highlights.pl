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
use File::Find::Rule;
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
GetOptions(\%opt, "cfg=s", "debug|d", "lsf|l") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my %cfg = ParseConfig(\%opt);

###############################################################################
#Main Program
###############################################################################
my @matlab_code;

########################################
# Full Filter Set 
########################################

#search for needed files
my @config_files = File::Find::Rule->file()->name( "*.m" )->in( catdir($cfg{exp_results_folder}, $cfg{movie_output_folder}) );
die "Found more than one config file:\n", join("\n\t",@config_files) if (scalar(@config_files) > 1);

my $invado_file = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, 'invado_data.csv');
die "Unable to find: $invado_file" if (! -e $invado_file);

my $area_lineage_ts = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, 
	$cfg{lineage_ts_folder},'Area.csv');
die "Unable to find: $area_lineage_ts" if (! -e $area_lineage_ts);

push @matlab_code, "make_final_puncta_highlight('$config_files[0]','$invado_file','$area_lineage_ts','output_dir','all_filt')";

########################################
# Filter Only By Local Diff
########################################

my $invado_file = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, 'local_invado_data.csv');
if (-e $invado_file) {
	push @matlab_code, "make_final_puncta_highlight('$config_files[0]','$invado_file','$area_lineage_ts','output_dir','local_diff')";
}

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'final_image_highlight');
$opt{error_file} = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'final_image_highlight', 'error.txt');
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

if (not($opt{config_only})) {
    &FA_job::run_matlab_progam(\@matlab_code,\%opt);
}

###############################################################################
#Functions
###############################################################################

###############################################################################
#Documentation
###############################################################################
