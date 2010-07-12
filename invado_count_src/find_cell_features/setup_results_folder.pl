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
use File::Which qw(which where);
use Cwd;

use Config::Adhesions qw(ParseConfig);
use Image::Stack;
use Math::Matlab::Extra;
use Emerald;
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l")
  or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

################################################################################
# Main Program
################################################################################

mkpath($cfg{individual_results_folder});

my @image_sets = ([qw(raw_mask_folder raw_mask_file)], [qw(adhesion_image_folder adhesion_image_file)],
				  [qw(gel_image_folder gel_image_file)]);

my @python_commands;
my $all_images_empty = 1;

foreach (@image_sets) {
    my $folder   = $cfg{ $_->[0] };
    my $out_file = $cfg{ $_->[1] };
    
	next if (not(defined($folder)));

    #Remove an ' marks used in config file to keep the folder name together
    $folder =~ s/\'//g;
	
	if ($opt{debug}) {
		print "Searching: $cfg{exp_data_folder}/$folder/\n";
	}
	
	#search for the image files currently present in the image folder
    my @image_files = sort <$cfg{exp_data_folder}/$folder/*>;

	#if we only find one we have to split the stacked tiff
	if (scalar(@image_files) == 1) {
		my @tiffsplit_which = which('tiffsplit');
		if (scalar(@tiffsplit_which) == 0) {
			die "Only found a single image file ($image_files[0]), which needs " . 
				"to be split, but couldn't find tiffsplit, please install tiffsplit " . 
				"and try again."
		} 
		
		my $tiff_stack_dirname = dirname($image_files[0]);

		if ($opt{debug}) {
			print("tiffsplit $image_files[0] $tiff_stack_dirname/\n");
			print("unlink $image_files[0]\n");
		} else {
			system("tiffsplit $image_files[0] $tiff_stack_dirname/ 2>/dev/null");
			system("unlink $image_files[0]\n");
		}
	
		&move_files_into_numbered_list($tiff_stack_dirname);
	}
    
	# after splitting any stacked tiffs we refresh the file listing in the data
	# directory
	@image_files = sort <$cfg{exp_data_folder}/$folder/*>;
    # my @image_files = map { $_ =~ s/\'//g; $_; } @image_files;
	
    $all_images_empty = 0 if (@image_files);

    if ($opt{debug}) {
        if (scalar(@image_files) > 0) {
            print "Image files found: $image_files[0] - $image_files[$#image_files]\n";
        } elsif (scalar(@image_files) == 0) {
            print "No image files found matching $cfg{exp_data_folder}/$folder, moving onto next image set.\n";
		}
        print "For Config Variable: ", $_->[0], "\n\n";
    }

	#jump to the next set of images if there weren't any image files found
	next if (not @image_files);

	push @python_commands, &build_python_commands(\@image_files, $out_file);
}
die "Unable to find any images to include in the new experiment" if $all_images_empty;

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'setup');
$opt{error_file} = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'setup', 'error.txt');
$opt{runtime} = "0:5";
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

&FA_job::send_general_program(\@python_commands,\%opt);

################################################################################
#Functions
################################################################################

sub move_files_into_numbered_list {
	my $dir = shift;
	
	my @files = sort <$dir/*>;
	my $number_of_digits = scalar(split("", scalar(@files)));
	my @output_files = map {
		catfile(dirname($files[$_]), sprintf("%0$number_of_digits" . "d", $_ + 1) . ".tif")
	} (0 .. $#files);
	
	for (0 .. $#files) {
		move($files[$_], $output_files[$_]);
	}
}

sub build_python_commands {
	my @src_image_files = @{$_[0]};
	my $output_file_name = $_[1];

	my $base_target = catdir($cfg{individual_results_folder});

	my @commands;
	my $image_num = 0;
	my $padding_length = scalar(split("", @src_image_files));
	foreach (@src_image_files) {
		if ($_ =~ /\/(\d+).tif/) {
			$image_num = $1;
		} else {
			$image_num++;
			# warn "Unable to find the image number in $_";
		}
		$image_num = sprintf("%0$padding_length" ."d", $image_num);

		my $target_name = catfile($base_target, $image_num, $output_file_name);

		if ($opt{debug}) {
			# print "mkpath(dirname($target_name))\n";
		} else {
			mkpath(dirname($target_name));
		}

		push @commands, "./write_normalized_image.py $_ $target_name";
	}
	return @commands;
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

=item * lsf: submit jobs through the emerald queuing system

=back

=head1 EXAMPLES

setup_results_folder.pl -cfg FA_config

OR

setup_results_folder.pl -cfg FA_config -d

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 5/24/2010

=cut
