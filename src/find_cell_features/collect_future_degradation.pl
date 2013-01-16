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
my @image_files   = <$cfg{individual_results_folder}/*/puncta_labeled.png>;
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

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'future_degrade');
$opt{error_file} = catfile($opt{error_folder}, 'error.txt');
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
		my $gel_file = catfile(dirname($file), "registered_gel.png");
		my $binary_shift = catfile(dirname($file), "binary_shift.png");
        
		die "Can't find puncta image file." if (not -e $file);
		die "Can't find gel image file." if (not -e $gel_file);
		die "Can't find binary shift image." if (not -e $binary_shift);

		$matlab_code[0] .= "find_future_degrade_values('$file','$gel_file','$binary_shift')\n";
    }

    return @matlab_code;
}
