#!/usr/bin/env perl

################################################################################
# Global Variables and Modules
################################################################################

use lib "../lib";
use lib "../lib/perl";

use strict;
use File::Path;
use File::Find;
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
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

################################################################################
# Main Program
################################################################################
if ($opt{lsf}) {
    #$0 - the name of the program currently running, used to protect against
    #future file name changes
    my @command = "$0 -cfg $opt{cfg}";
    $opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'global_min_max');
    if (defined $cfg{job_group}) {
        $opt{job_group} = $cfg{job_group};
    }
	$opt{queue} = "week";

    &FA_job::send_general_program(\@command,\%opt);
    
    exit;
}

our @gel_min_max = (99999999,-99999999);
our @puncta_min_max = (99999999,-99999999);

find(\&find_min_maxes,(dirname($cfg{exp_results_folder})));

if (not $opt{debug}) {
	find(\&output_min_maxes,(dirname($cfg{exp_results_folder})));
} else {
	print join(" ", @gel_min_max), "\n";
	print join(" ", @puncta_min_max), "\n";
}

################################################################################
# Functions
################################################################################

sub find_min_maxes {
	if ($_ eq 'gel_image_range.csv') {
		open INPUT, "gel_image_range.csv";
		my @min_max = split(",",<INPUT>);
		close INPUT;
		
		$gel_min_max[0] = $min_max[0] if ($min_max[0] < $gel_min_max[0]);
		$gel_min_max[1] = $min_max[1] if ($min_max[1] > $gel_min_max[1]);
	}
	if ($_ eq 'puncta_image_range.csv') {
		open INPUT, "puncta_image_range.csv";
		my @min_max = split(",",<INPUT>);
		close INPUT;
		
		if ($min_max[0] < $puncta_min_max[0]) {
			$puncta_min_max[0] = $min_max[0];
			print "Found new min: $min_max[0]\n" if $opt{debug};
		}

		$puncta_min_max[1] = $min_max[1] if ($min_max[1] > $puncta_min_max[1]);
	}
}

sub output_min_maxes {
	if ($_ eq 'gel_image_range.csv') {
		open OUTPUT, ">gel_image_range.csv";
		print OUTPUT join(",",@gel_min_max);
		close OUTPUT;
	}
	if ($_ eq 'puncta_image_range.csv') {
		open OUTPUT, ">puncta_image_range.csv";
		print OUTPUT join(",",@puncta_min_max);
		close OUTPUT;
	}
}
