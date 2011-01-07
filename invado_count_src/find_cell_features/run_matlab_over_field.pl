#!/usr/bin/env perl

################################################################################
# Global Variables and Modules
################################################################################

use lib "../lib";
use lib "../lib/perl";

use strict;
use File::Path;
use File::Spec::Functions;
use File::Basename;
use Image::ExifTool;
use Getopt::Long;
use Data::Dumper;

use Config::Adhesions;
use Math::Matlab::Extra;
use Emerald;
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "script=s", "debug|d", "lsf|l", "queue=s", "memsize=s") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};
die "Can't find script specified on the command line" if not exists $opt{script};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

if ($opt{script} =~ /(.*)\.m/) {
	$opt{script} = $1;
}

################################################################################
# Main Program
################################################################################
my @matlab_code = ("$opt{script}('$cfg{individual_results_folder}')\n");

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, $opt{script});
$opt{error_file} = catfile($opt{error_folder}, 'error.txt');
$opt{runtime} = "1";
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

&FA_job::run_matlab_progam(\@matlab_code,\%opt);