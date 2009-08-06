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

my @matlab_code = &create_all_matlab_commands;

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'FA');
$opt{error_file} = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'FA', 'error.txt');
$opt{runtime} = "1";
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

&FA_job::run_matlab_progam(\@matlab_code,\%opt);

################################################################################
#Functions
################################################################################

sub create_all_matlab_commands {
    my @matlab_code;

    my @image_files = <$cfg{individual_results_folder}/*/$cfg{adhesion_image_file}>;
    foreach my $file_name (@image_files) {
        my $cell_mask = catfile(dirname($file_name), $cfg{cell_mask_file});

        my $extra_opt = "";
        if (defined $cfg{filter_thresh}) {
            $extra_opt .= ",'filter_thresh',$cfg{filter_thresh}";
        }
        if (defined $cfg{scale_filter_thresh}) {
            $extra_opt .= ",'scale_filter_thresh',$cfg{scale_filter_thresh}";
        }

        if (-e $cell_mask) {
            $matlab_code[0] .= "find_focal_adhesions('$file_name','cell_mask','$cell_mask'$extra_opt)\n";
        } else {
            $matlab_code[0] .= "find_focal_adhesions('$file_name'$extra_opt)\n";
        }
    }

    return @matlab_code;
}

################################################################################
#Documentation
################################################################################

=head1 NAME

collect_fa_image_set.pl - Executes the MATLAB programs designed collect the
masks which define the focal adhesions

=head1 SYNOPSIS

collect_fa_image_set.pl -cfg FA_config

=head1 Description

This program is used to create all the binary mask files which define the
location of the focal adhesions, which will be used in subsequent steps. The
primary logic of the program is in a set of MATLAB scripts which do all the
image analysis/writing and also collects properties of the focal adhesions and
writes those to file.

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=back

Optional parameter(s):

=over 

=item * debug or d: print debuging information during program execution

=item * emerald: submit jobs through the emerald queuing system

=back

=head1 EXAMPLES

collect_fa_image_set.pl -cfg FA_config

OR

collect_fa_image_set.pl -cfg FA_config -d

=head1 SEE ALSO

collect_mask_set.pl: similar program designed to collect the binary mask that
locates the intracellular area

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 7/3/2008 

=cut
