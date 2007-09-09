#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################

use strict;
use File::Path;
use Image::ExifTool;
use Math::Matlab::Local;

my $experimental_data_folder = "../../data/";
#$experimental_data_folder = "/Volumes/Data/projects/focal_adhesions/data/";
my $experiment_name = "time_series_1";

my $results_folder = "../../results/";

my $adhesion_image_prefix = "EGFP-Paxillin_";

my $matlab_code_path = "..";

my $matlab_wrapper = Math::Matlab::Local->new({
	cmd => '/usr/bin/matlab -nodisplay -nojvm -nosplash',
});
###############################################################################
# Main Program
###############################################################################

mkpath($results_folder . $experiment_name . "/individual_pictures");

my $adhesion_image_location = $experimental_data_folder . $experiment_name . "/" . $adhesion_image_prefix;
my @adhesion_image_files = <$adhesion_image_location*>;

print @adhesion_image_files;
foreach my $file_name (@adhesion_image_files) {
	my $total_images = &get_image_stack_number($file_name);
	foreach my $image_num (1..$total_images) {
		my $padded_num = sprintf("%0" . length($total_images) . "d", $image_num);
		my $output_path = $results_folder . $experiment_name . "/individual_pictures/$padded_num";
		mkpath($output_path);
		my $matlab_code = $matlab_code . "find_focal_adhesions('$file_name',$image_num,'$output_path')\n";
		$matlab_wrapper->execute($matlab_code);
	}
}

###############################################################################
# Functions
###############################################################################

sub get_image_stack_number {
	my $image_file = pop;
	
	my $image_info = new Image::ExifTool;
	$image_info->ExtractInfo($image_file);
	my @tag_list = $image_info->GetFoundTags($image_file);
	my $image_count = 0;
	foreach (@tag_list) {
		if (/\((\d+)\)/) {
			if ($1 > $image_count) {
				$image_count = $1;
			}
		}
	}

	$image_count++;
	return $image_count;
}
