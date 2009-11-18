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
my @image_files   = <$cfg{individual_results_folder}/*/$cfg{adhesion_image_file}>;
die "Expected to find the same number of image files as folders in the results directory ($cfg{individual_results_folder})."
  if (scalar(@image_files) != scalar(@image_folders));

if ($opt{debug}) {
    if (scalar(@image_files) > 1) {
        print "Focal image files found: $image_files[0] - $image_files[$#image_files]\n";
    } elsif (scalar(@image_files) == 0) {
        warn "Couldn't find any focal image files in $cfg{individual_results_folder} subfolders\n\n";
    } else {
        print "Focal image file found: $image_folders[0]\n";
    }
}

my @matlab_code = &create_all_matlab_commands(@image_files);

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'puncta_props');
$opt{error_file} = catfile($opt{error_folder}, 'puncta_props', 'error.txt');
$opt{runtime} = "1";
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

&FA_job::run_matlab_progam(\@matlab_code,\%opt);

################################################################################
#Functions
################################################################################

sub create_all_matlab_commands {
	my @image_files = @_;
    my @matlab_code;

    foreach my $file (@image_files) {
		
		my $gel_binary_file = catfile(dirname($file), "degradation_binary.png");
		my $puncta_file = catfile(dirname($file), "puncta_labeled.png");

        my $extra_opt = "";
        if (not(-e $gel_binary_file)) {
			die "Unable to find the gel binary file in: $gel_binary_file";
        }
        if (not(-e $puncta_file)) {
			die "Unable to find the labeled puncta file in: $file";
        }
        
		$matlab_code[0] .= "find_adhesion_properties('$file','$puncta_file','$gel_binary_file'$extra_opt)\n";
    }

    return @matlab_code;
}
