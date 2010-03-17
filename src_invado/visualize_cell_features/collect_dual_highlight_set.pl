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
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

################################################################################
# Main Program
################################################################################

my @image_folders = <$cfg{individual_results_folder}/*>;


if ($opt{debug}) {
    if (scalar(@image_folders) > 1) {
        print "Focal image files found: $image_folders[0] - $image_folders[$#image_folders]\n";
    } elsif (scalar(@image_folders) == 0) {
        warn "Couldn't find any focal image files in $cfg{individual_results_folder} subfolders\n\n";
    } else {
        print "Focal image file found: $image_folders[0]\n";
    }
}

my @matlab_code = &create_all_matlab_commands(@image_folders);

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'dual_highlight');
$opt{error_file} = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'dual_highlight', 'error.txt');
$opt{runtime} = "1";
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

&FA_job::run_matlab_progam(\@matlab_code,\%opt);

################################################################################
#Functions
################################################################################

sub create_all_matlab_commands {
	my @image_folders = @_;
    my @matlab_code;

    foreach my $folder_name (@image_folders) {
        my $extra_opt = "";

		$matlab_code[0] .= "build_dual_highlight_image('$folder_name'$extra_opt)\n";
    }

    return @matlab_code;
}

################################################################################
#Documentation
################################################################################
