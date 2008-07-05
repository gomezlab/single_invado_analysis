#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################

use lib "../lib";
use lib "../lib/perl";

use strict;
use File::Path;
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
$opt{emerald} = 0;
$opt{emerald_stdout} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "emerald", "emerald_stdout");

if ($opt{emerald} == 1 && $opt{emerald_stdout} == 1) {
    die "Please specify only one of emerald or emerald_stdout"; 
}
die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

###############################################################################
# Main Program
###############################################################################

my @image_folders = <$cfg{individual_results_folder}/*>;
my @image_files = <$cfg{individual_results_folder}/*/$cfg{raw_mask_file}>;

if ($opt{debug}) {
    if (scalar(@image_files) > 1) {
        print "Cell mask files found: $image_files[0] - $image_files[$#image_files]\n";
    } elsif ( scalar(@image_files) == 0) {
        die "Couldn't find any cell mask files in $cfg{individual_results_folder} subfolders\n\n";
    } else {
        print "Cell mask file found: $image_folders[0]\n";
    }
}

my @matlab_code = &create_matlab_code;

my $error_folder = catdir($cfg{exp_results_folder}, $cfg{matlab_errors_folder}, 'mask_set');
my $error_file = catfile($cfg{exp_results_folder}, $cfg{matlab_errors_folder}, 'mask_set', 'error.txt');

mkpath($error_folder);
my %emerald_opt = ("folder", $error_folder);
if ($opt{emerald}) {
    my @matlab_code = sort @matlab_code;
    my @commands = &Emerald::create_emerald_Matlab_commands(\@matlab_code,\%emerald_opt);
    &Emerald::send_emerald_commands(\@commands);
}
elsif ($opt{emerald_stdout}) {
    for (sort @image_folders) {
        my @command = "./collect_mask_image.pl -cfg $opt{cfg} -folder $_\n";
        @command = &Emerald::create_general_emerald_command(\@command);
        print @command, "\n";
    }
} else {
    &Math::Matlab::Extra::execute_commands(\@matlab_code, $error_file);
}

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
            warn "Skipping file: $file_name\n", "Unable to find image number.";
            next;
        }

        next if grep $i_num == $_, @{ $cfg{exclude_image_nums} };

        my $out_file = catfile(dirname($file_name), $cfg{cell_mask_file});
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

=item * emerald: submit jobs through the emerald queuing system

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
