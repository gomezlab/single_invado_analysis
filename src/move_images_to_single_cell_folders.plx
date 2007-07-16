#!/usr/bin/perl -w

use strict;
use File::Copy;
use File::Path;

###############################################################################
# Program Constants
###############################################################################

my $data_directory = '../data/original\ data/split/';
my $output_directory = '../data/each_cell/';
my %time_to_image_count = ( "20min" => 1, "40min" => 2, "60min" => 3, 
							"180min" => 4, "240min" => 5);
my $total_cell_number = 0;
my $debug = 0;

###############################################################################
# Main Program
###############################################################################

foreach (<$data_directory*/*/*.png>) {
	if (/$data_directory(.*?)\/.*\/(.*?)(\d+)\.png/) {
		$3 =~ /0+(\d+)/;
		my $this_num = $1;
		if ($this_num > $total_cell_number) {
			$total_cell_number = $this_num;
		}
	}
}

foreach (<$data_directory*/*/*.png>) {
	if (/$data_directory(.*?)\/.*\/(.*?)(\d+)\.png/) {
		my $new_file_name = $2 . sprintf("%0" . length(keys %time_to_image_count) . "d", $time_to_image_count{$1}) . ".png";

		my $new_folder_name = $output_directory . "$3/";
		my $full_file_name = $new_folder_name . $new_file_name;
		
		print "$full_file_name - $_\n" if $debug;

		if (! -e $new_folder_name) {	
			mkpath($new_folder_name) or die "$! - $_";
		}
		if (! -e $full_file_name) {
			copy($_,$full_file_name) or die "$! - $_";
		}
	}
}
