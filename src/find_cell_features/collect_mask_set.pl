#!/usr/bin/env perl


###############################################################################
# Global Variables and Modules
###############################################################################

use strict;
use File::Path;
use File::Spec::Functions;
use File::Basename;
use Image::ExifTool;
use Math::Matlab::Local;
use Getopt::Long;

use lib "../lib";
use Config::Adhesions;
use Image::Stack;
use Math::Matlab::Extra;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

###############################################################################
# Main Program
###############################################################################

my @image_files = <$cfg{individual_results_folder}/*/$cfg{cell_mask_file}>;

if ($opt{debug}) {
    if (scalar(@image_files) > 1) {
        print "Cell mask files found: $image_files[0] - $image_files[$#image_files]\n";
    } elsif ( scalar(@image_files) == 0) {
        die "Couldn't find any cell mask files in $cfg{individual_results_folder} subfolders\n\n";
    } else {
        print "Cell mask file found: $image_files[0]\n";
    }
}

my @matlab_code = &create_matlab_code;

my $error_file = catdir($cfg{exp_results_folder}, $cfg{matlab_errors_folder}, $cfg{cell_mask_errors_file});

&Math::Matlab::Extra::execute_commands(\@matlab_code, $error_file);

###############################################################################
#Functions
###############################################################################

sub create_matlab_code {
    my @matlab_code;

    foreach my $file_name (@image_files) {
        my $i_num;
        if ($file_name =~ /$cfg{individual_results_folder}\/(\d+)\//) {
            $i_num = $1;
        } else {
            die "Skipping file: $file_name\n", "Unable to find image number.";
            next;
        }

        next if grep $i_num == $_, @{ $cfg{exclude_image_nums} };

        my $padded_num = sprintf("%0" . length(scalar(@image_files)) . "d", $i_num);

        my $out_file = catfile(dirname($file_name), "cell_mask.png");
        $matlab_code[0] .= "find_cell_mask('$file_name','$out_file')\n";
    }

    return @matlab_code;
}

###############################################################################
#Documentation
###############################################################################

=head1 NAME

collect_mask_set.pl - Executes the MATLAB programs designed collect the binary
masks which define the intracellular region of a cell

=head1 SYNOPSIS

collect_mask_set.pl -cfg FA_config

=head1 DESCRIPTION

This program is used to create all the binary mask files which define the
location of the intracellular space, which will be used in subsequent steps. The
primary logic of the program is in a set of MATLAB scripts which do all the
threshold selection and image creation/writing. This perl script is more or less
an intelligent wrapper that sets up all the MATLAB commands.  

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=back

Optional parameter(s):

=over 

=item * debug or d: print debuging information during program execution

=back

=head1 EXAMPLES

collect_mask_set.pl -cfg FA_config

OR

collect_mask_set.pl -cfg FA_config -d

=head1 SEE ALSO

collect_fa_image_set.pl: similar program designed to collect the binary mask
that locates focal adhesions

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 4/10/2008

=cut
