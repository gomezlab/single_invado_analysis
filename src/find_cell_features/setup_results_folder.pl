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
use File::Copy;
use Image::ExifTool;
use Math::Matlab::Local;
use Getopt::Long;
use Data::Dumper;

use Config::Adhesions qw(ParseConfig);
use Image::Stack;
use Math::Matlab::Extra;
use Emerald;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "emerald", "emerald_stdout");
die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

################################################################################
# Main Program
################################################################################

mkpath($cfg{individual_results_folder});

my @image_sets = ([qw(cell_mask_folder raw_mask_file)], 
                  [qw(adhesion_image_folder adhesion_image_file)]);
my @matlab_code;
my $all_images_empty = 1;

foreach (@image_sets) {
    my $folder   = $cfg{ $_->[0] };
    my $out_file = $cfg{ $_->[1] };
    
    next if (not(defined($folder)));
    
    #Remove an ' marks used in config file to keep the folder name together
    $folder =~ s/\'//g;

    my @image_files = sort <$cfg{exp_data_folder}/$folder/*>;
    my @image_files = map { $_=~s/\'//g; $_; } @image_files;

    #Move all the files with spaces in their names, MATLAB on emerald doesn't
    #like them
    @image_files = &remove_file_name_spaces(@image_files);
    $all_images_empty = 0 if (@image_files);
    
    if ($opt{debug}) {
        if (scalar(@image_files) > 1) {
            print "Image files found: $image_files[0] - $image_files[$#image_files]\n";
        } elsif (scalar(@image_files) == 0) {
            print "No image files found matching $cfg{exp_data_folder}/$folder, moving onto next image set.\n";
            next;
        } else {
            print "Image file found: $image_files[0]\n";
        }
        print "For Config Variable: ", $_->[0], "\n\n";
    } else {
        next if (not @image_files);
    }

    push @matlab_code, &create_matlab_code(\@image_files, $folder, $out_file);
}
die "Unable to find any images to include in the new experiment" if $all_images_empty;

my $error_folder = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'setup');
my $error_file = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'setup', 'error.txt');
mkpath($error_folder);

my %emerald_opt = ("folder", $error_folder);
if ($opt{emerald}) {
    my @commands = &Emerald::create_LSF_Matlab_commands(\@matlab_code, \%emerald_opt);
    &Emerald::send_LSF_commands(\@commands);
} else {
    &Math::Matlab::Extra::execute_commands(\@matlab_code, $error_file);
}

################################################################################
#Functions
################################################################################

sub create_matlab_code {
    my @image_files = @{ $_[0] };
    my $folder      = $_[1];
    my $out_file    = $_[2];

    my @image_stack_count = map { Image::Stack::get_image_stack_number($_) } @image_files;

    my @matlab_code;
    if (grep { $_ > 1 } @image_stack_count) {
        if (scalar(@image_files) > 1) {
            die "Found more than one image stack in: ", join(", ", @image_files), "\n",
              "Expected single image stack or multiple non-stacked files\n";
        }
        @matlab_code = &create_matlab_code_stack(\@image_files, $out_file);
    } else {
        @matlab_code = &create_matlab_code_single(\@image_files, $out_file);
    }
    return @matlab_code;
}

sub create_matlab_code_stack {
    my @image_files = @{ $_[0] };
    my $out_file    = $_[1];

    my @matlab_code;

    my $total_stack_images = Image::Stack::get_image_stack_number($image_files[0]);
    foreach my $i_num (1 .. $total_stack_images) {
        next if grep $i_num == $_, @{ $cfg{exclude_image_nums} };

        my $padded_num = sprintf("%0" . length($total_stack_images) . "d", $i_num);

        my $output_path = catdir($cfg{individual_results_folder}, $padded_num);
        mkpath($output_path);
        my $final_out_file = catfile($output_path, $out_file);
        $matlab_code[0] .=
          "write_normalized_image('$image_files[0]','$final_out_file','I_num',$i_num);\n";
    }
    return @matlab_code;
}

sub create_matlab_code_single {
    my @image_files = @{ $_[0] };
    my $out_file    = $_[1];

    my @matlab_code;
    

    foreach my $file_name (@image_files) {
        my $i_num;
        my $original_i_num;
        if ($file_name =~ /.*?(\d+)\./) {
            $original_i_num = $1;
            $i_num = (grep $file_name eq $image_files[$_ - 1], (1 .. $#image_files + 1))[0];
        } else {
            warn "Unable to find image number in: $file_name, skipping this image.";
            next;
        }

        next if grep $i_num == $_, @{ $cfg{exclude_image_nums} };
        next if grep $original_i_num == $_, @{ $cfg{exclude_image_nums} };

        my $padded_num = sprintf("%0" . length(scalar(@image_files)) . "d", $i_num);

        my $output_path = catdir($cfg{individual_results_folder}, $padded_num);
        mkpath($output_path);
        my $final_out_file = catfile($output_path, $out_file);
        $matlab_code[0] .= "write_normalized_image('$file_name','$final_out_file');\n";
    }
    return @matlab_code;
}

sub remove_file_name_spaces {
    my @new_files;
    for (@_) {
        my $dir = dirname($_);
        my $file_name = basename($_);
        
        if ($file_name =~ /\s/) {
            if ($file_name =~ /(\d+\..*)/) {
                if( move($_,catfile($dir,$1)) ) {
                    push @new_files, catfile($dir,$1);
                }
            } else {
                warn "Unable to determine image number for file: $_";
            }
        } else {
            push @new_files, $_;
        }

    }
    return @new_files;
}

################################################################################
#Documentation
################################################################################

=head1 NAME

setup_results_folder.pl - Move all the raw data files into the proper locations
in the results folder

=head1 SYNOPSIS

setup_results_folder.pl -cfg FA_config

=head1 Description

Since the the data sets being used come in multiple forms, mostly split and
stacked image sequences, it became easier to write a single function to move all
the data files into a standard results directory. Then all the downstream
programs would have an easier time trying to find specific images. This program
scans through a the data directory for an experiment and moves all the files
into the correct locations.

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

setup_results_folder.pl -cfg FA_config

OR

setup_results_folder.pl -cfg FA_config -d

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 4/10/2008 

=cut
