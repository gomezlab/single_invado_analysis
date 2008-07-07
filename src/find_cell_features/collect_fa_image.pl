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
GetOptions(\%opt, "cfg|c=s", "debug|d", "folder=s");
die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

################################################################################
# Main Program
################################################################################

my @matlab_code = &create_matlab_code;

my $error_folder = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'FA');
my $error_file = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'FA', 'error.txt');
my %emerald_opt = ("folder", $error_folder);
mkpath($error_folder);

my @commands = &Emerald::create_LSF_Matlab_commands(\@matlab_code,\%emerald_opt);
&Emerald::send_LSF_commands(\@commands);

################################################################################
#Functions
################################################################################

sub create_matlab_code {
    my $command;

    my @image_files = <$opt{folder}/$cfg{adhesion_image_file}>;
    
    if (scalar(@image_files) == 0) {
    } elsif (scalar(@image_files) > 1) {
       die "More than one image file found (". scalar(@image_files) . "), quiting";
    } else {
        my $cell_mask = catfile(dirname($image_files[0]), $cfg{cell_mask_file});
        
        my $extra_opt = "";
        if (defined $cfg{filter_thresh}) {
            $extra_opt .= ",'filter_thresh',$cfg{filter_thresh}";
        }

        if (-e $cell_mask) {
            $command = "find_focal_adhesions('$image_files[0]','cell_mask','$cell_mask'$extra_opt)\n";
        } else {
            $command = "find_focal_adhesions('$image_files[0]'$extra_opt)\n";
        }
    }
    return $command;
}

################################################################################
#Documentation
################################################################################
