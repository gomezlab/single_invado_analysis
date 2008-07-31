#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################

use lib "../lib";
use lib "../lib/perl";

use strict;
use File::Path;
use Data::Dumper;
use File::Spec::Functions;
use File::Basename;
use Image::ExifTool;
use Math::Matlab::Local;
use Getopt::Long;

use Config::Adhesions;
use Image::Stack;
use Math::Matlab::Extra;
use Emerald;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "folder=s", "debug|d", "emerald|e");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

###############################################################################
# Main Program
###############################################################################

my @matlab_code = &create_matlab_command;

my $error_folder = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'mask_set');
my $error_file = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'mask_set', 'error.txt');
my %emerald_opt = ("folder", $error_folder);
mkpath($error_folder);

if (scalar(@matlab_code) >= 1) {
    my @commands = &Emerald::create_LSF_Matlab_commands(\@matlab_code,\%emerald_opt);
    &Emerald::send_LSF_commands(\@commands);
}
###############################################################################
#Functions
###############################################################################

sub create_matlab_command {
    my @command;
    
    my @image_files = grep -e $_, <$opt{folder}/$cfg{raw_mask_file}>;

    if (scalar(@image_files) == 0) {
    } elsif (scalar(@image_files) > 1) {
       die "More than one image file found (". scalar(@image_files) . "), quiting";
    } else {
        my $out_file = catfile(dirname($image_files[0]), $cfg{cell_mask_file});
        push @command, "find_cell_mask('$image_files[0]','$out_file')";
    }
    return @command;
}


###############################################################################
#Documentation
###############################################################################
