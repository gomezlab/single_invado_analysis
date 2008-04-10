#!/usr/bin/env perl

=head1 NAME

collect_fa_image_set.pl - Executes the MATLAB programs designed collect the
binary masks which define the focal adhesions

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

=item * fa_debug: only execute the focal adhesion finding MATLAB code for one of
the images

=back

=head1 EXAMPLES

collect_fa_image_set.pl -cfg FA_config

OR

collect_fa_image_set.pl -cfg FA_config -d

OR

collect_fa_image_set.pl -c FA_config -debug -fa_debug

=head1 SEE ALSO

collect_mask_set.pl: similar program designed to collect the binary mask that
locates the intracellular area

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 4/10/2008 

=cut

################################################################################
# Global Variables and Modules
################################################################################

use strict;
use File::Path;
use Image::ExifTool;
use Math::Matlab::Local;
use Getopt::Long;
use File::Spec::Functions;

use lib "../lib";
use Config::Adhesions;
use Image::Stack;
use Math::Matlab::Extra;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "fa_debug");
die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

my $matlab_wrapper;
if (defined $cfg{matlab_executable}) {
    $matlab_wrapper = Math::Matlab::Local->new({ cmd => "$cfg{matlab_executable} -nodisplay -nojvm -nosplash", });
} else {
    $matlab_wrapper = Math::Matlab::Local->new();
}

################################################################################
# Main Program
################################################################################

my @focal_image_files = <$cfg{individual_results_folder}/*/$cfg{adhesion_image_file}>;

if ($opt{debug}) {
    if (scalar(@focal_image_files) > 1) {
        print "Adhesion files found: $focal_image_files[0] - $focal_image_files[$#focal_image_files]\n";
    } else {
        print "Adhesion file found: $focal_image_files[0]\n";
    }
}

my @matlab_code = &create_matlab_code;

my $error_file = catfile($cfg{exp_results_folder}, $cfg{matlab_errors_folder}, $cfg{adhesion_errors_file});
&Math::Matlab::Extra::execute_commands($matlab_wrapper, \@matlab_code, $error_file);

################################################################################
#Functions
################################################################################

sub create_matlab_code {
    my @matlab_code;

    foreach my $file_name (@focal_image_files) {
        my $i_num;
        if ($file_name =~ /$cfg{individual_results_folder}\/(\d+)\//) {
            $i_num = $1;
        } else {
            die "Skipping file: $file_name\n", "Unable to find image number.";
            next;
        }
        next if ($i_num > 1 && $opt{fa_debug});

        next if grep $i_num == $_, @{ $cfg{exclude_image_nums} };

        my $cell_mask = catfile(dirname($file_name), "cell_mask.png");

        $matlab_code[0] .= "find_focal_adhesions('$file_name','$cell_mask')\n";
    }

    return @matlab_code;
}
