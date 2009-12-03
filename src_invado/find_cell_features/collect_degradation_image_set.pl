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
my @image_files   = sort <$cfg{individual_results_folder}/*/registered_gel.png>;
die "Expected to find the same number of image files as folders in the results directory ($cfg{individual_results_folder})." if (scalar(@image_files) != scalar(@image_folders));

if ($opt{debug}) {
    if (scalar(@image_files) > 1) {
        print "Gel image files found: $image_files[0] - $image_files[$#image_files]\n";
    } elsif (scalar(@image_files) == 0) {
        warn "Couldn't find any gel image files in $cfg{individual_results_folder} subfolders\n\n";
    } else {
        print "Gel image file found: $image_folders[0]\n";
    }
}

my @matlab_code = &create_all_matlab_commands(@image_files);

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'Gel_degrade');
$opt{error_file} = catfile($opt{error_folder}, 'error.txt');
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

    foreach my $file_name (@image_files) {
        my $extra_opt = "";

		my $directory_name = dirname($file_name);
		my $registration_binary = catfile($directory_name,'binary_shift.png');
		if (-e $registration_binary) {
			$extra_opt = ",'registration_binary','$registration_binary'";
		} else {
			die "Unable to find registration binary in: $file_name";
		}

		$matlab_code[0] .= "find_degraded_regions('$file_name','$image_files[0]'$extra_opt)\n";
    }

    return @matlab_code;
}

################################################################################
#Documentation
################################################################################
