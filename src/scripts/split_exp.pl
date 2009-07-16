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
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l", "split_num=i")
  or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};
die "The split number must be specified on the command line" if not exists $opt{split_num};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

################################################################################
# Main Program
################################################################################

my @image_sets = ([qw(raw_mask_folder raw_mask_file)], [qw(adhesion_image_folder adhesion_image_file)]);
for (@image_sets) {
	&move_target_image_set(@{$_});
}

################################################################################
#Functions
################################################################################

sub move_target_image_set {
	my $target_dir = shift @_;
	my $source_file = shift @_;
		
	if (not exists $cfg{$target_dir} || not exists $cfg{$source_file}) {
		die;
	}

	my @image_files = sort <$cfg{individual_results_folder}/*/$cfg{$source_file}>;
	my @image_dirs = sort <$cfg{individual_results_folder}/*>;
	my @image_nums = map basename($_), @image_dirs;
	
	if (scalar(@image_files) == 0 && scalar(@image_dirs) != 0) {
		next;
	}

	if (scalar(@image_files) != scalar(@image_dirs) || 
		scalar(@image_nums) != scalar(@image_files)) {
		die "Expected same number of image files as image dirs, instead found " .
			 scalar(@image_files) . " and " . scalar(@image_dirs);
	}

	my $max_digit_count = length($image_nums[-1]);

	my @split_index = grep $image_nums[$_] == $opt{split_num}, (0 .. $#image_nums);
	die if (scalar(@split_index) > 1);

	mkpath(catdir($cfg{data_folder}, $cfg{exp_name} . "_pre"));
	for (0 .. $split_index[0]) {
		my $image_name = sprintf("%0" . $max_digit_count . "d", $_ + 1);
		my $target_file = catdir($cfg{data_folder}, $cfg{exp_name} . "_pre", 
								 $cfg{$target_dir}, $image_name . ".png");
		mkpath(dirname($target_file));
		copy($image_files[$_], $target_file);
	}
	mkpath(catdir($cfg{data_folder}, $cfg{exp_name} . "_post"));
	for (($split_index[0] + 1) .. $#image_files) {
		my $image_name = sprintf("%0" . $max_digit_count . "d", $_ + 1 - $opt{split_num});
		my $target_file = catdir($cfg{data_folder}, $cfg{exp_name} . "_post", 
								 $cfg{$target_dir}, $image_name . ".png");
		mkpath(dirname($target_file));
		copy($image_files[$_], $target_file);
	}
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
