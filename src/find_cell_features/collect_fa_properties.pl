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

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "emerald|e", "emerald_debug|e_d") or die;

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

my $error_folder = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'FA_props');
my $error_file = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'FA_props', 'error.txt');

mkpath($error_folder);
my %emerald_opt = ("folder" => $error_folder, "runtime" => "1");
if ($opt{emerald} || $opt{emerald_debug}) {
    my @lsf_command = &Emerald::create_LSF_Matlab_commands(\@matlab_code, \%emerald_opt);
    if ($opt{emerald_debug}) {
        print join("\n", @lsf_command);
    } else {
        &Emerald::send_LSF_commands(\@lsf_command);
    }
} else {
    if ($opt{debug}) {
        print join("\n", @matlab_code);
    } else {
        &Math::Matlab::Extra::execute_commands(\@matlab_code, $error_file);
    }
}

################################################################################
#Functions
################################################################################

sub create_all_matlab_commands {
    my @matlab_code;

    my @raw_image_files = grep -e $_, <$cfg{individual_results_folder}/*/$cfg{adhesion_image_file}>;
    my @adhesion_image_files = grep -e $_, <$cfg{individual_results_folder}/*/adhesions.png>;
    die "Expected to find equal number of raw data adhesion image and adhesion mask images." if (scalar(@raw_image_files) != scalar(@adhesion_image_files));
    foreach (0..$#raw_image_files) {
        my $raw_image_file = $raw_image_files[$_];
        my $adhesion_image_file = $adhesion_image_files[$_];
        my $cell_mask = catfile(dirname($raw_image_file), $cfg{cell_mask_file});
        
        my $protrusion_file = catfile($cfg{exp_results_folder}, $cfg{protrusion_data_file});
        
        my $extra_opt = "";
        if (-e $cell_mask) {
            $extra_opt .= ",'cell_mask','$cell_mask'";
        }
        if (-e $protrusion_file) {
            $extra_opt .= ",'protrusion_file','$protrusion_file'";
        }
        
        my $this_command = "find_adhesion_properties('$raw_image_file','$adhesion_image_file','i_num',";
        $this_command .= $_ + 1 . "$extra_opt)\n";

        $matlab_code[0] .= $this_command;
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
