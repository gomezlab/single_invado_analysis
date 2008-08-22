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
GetOptions(\%opt, "cfg|c=s", "debug|d", "fa_debug", "folder=s", "emerald|e", "emerald_debug|e_d") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

################################################################################
# Main Program
################################################################################

my @image_folders = <$cfg{individual_results_folder}/*>;
my @image_files   = <$cfg{individual_results_folder}/*/$cfg{adhesion_image_file}>;

if ($opt{debug}) {
    if (scalar(@image_files) > 1) {
        print "Focal image files found: $image_files[0] - $image_files[$#image_files]\n";
    } elsif (scalar(@image_files) == 0) {
        warn "Couldn't find any focal image files in $cfg{individual_results_folder} subfolders\n\n";
    } else {
        print "Focal image file found: $image_folders[0]\n";
    }
}

my @matlab_code;
if (exists($opt{folder})) {
    @matlab_code = &create_single_matlab_command;
} elsif ($opt{emerald}) {
} else {
    @matlab_code = &create_all_matlab_commands;
}

my $error_folder = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'FA');
my $error_file = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'FA', 'error.txt');

mkpath($error_folder);
my %emerald_opt = ("folder" => $error_folder, "runtime" => "0:5");
if (exists($opt{folder})) {
    $emerald_opt{"runtime"} = "4";
    my @command = &Emerald::create_LSF_Matlab_commands(\@matlab_code, \%emerald_opt);
    if ($opt{emerald_debug}) {
        print join("\n", @command);
    } else {
        &Emerald::send_LSF_commands(\@command);
    }
} elsif ($opt{emerald} || $opt{emerald_debug}) {
    my @command;
    for (sort @image_folders) {
        push @command, "$0 -cfg $opt{cfg} -folder $_\n";
    }
    @command = &Emerald::create_general_LSF_commands(\@command, \%emerald_opt);
    if ($opt{emerald_debug}) {
        print join("\n", @command);
    } else {
        &Emerald::send_LSF_commands(\@command);
    }
} else {
    &Math::Matlab::Extra::execute_commands(\@matlab_code, $error_file);
}

################################################################################
#Functions
################################################################################

sub create_all_matlab_commands {
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

        my $cell_mask = catfile(dirname($file_name), "cell_mask.png");

        my $extra_opt = "";
        if (defined $cfg{filter_thresh}) {
            $extra_opt .= ",'filter_thresh',$cfg{filter_thresh}";
        }

        if (-e $cell_mask) {
            $matlab_code[0] .= "find_focal_adhesions('$file_name','cell_mask','$cell_mask'$extra_opt)\n";
        } else {
            $matlab_code[0] .= "find_focal_adhesions('$file_name'$extra_opt)\n";
        }
    }

    return @matlab_code;
}

sub create_single_matlab_command {
    my $command;

    my @image_files = grep -e $_, <$opt{folder}/$cfg{adhesion_image_file}>;

    if (scalar(@image_files) == 0) {
    } elsif (scalar(@image_files) > 1) {
        die "More than one image file found (" . scalar(@image_files) . "), quiting";
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
