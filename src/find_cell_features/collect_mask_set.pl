#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################

use strict;
use File::Path;
use Image::ExifTool;
use Math::Matlab::Local;
use Config::General qw/ParseConfig/;
use Getopt::Long;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg=s", "debug");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my %cfg = &get_config;

my $matlab_wrapper = Math::Matlab::Local->new({ cmd => '/usr/bin/matlab -nodisplay -nojvm -nosplash', });

###############################################################################
# Main Program
###############################################################################

mkpath($cfg{exp_result_single_folder});

my @cell_mask_files    = <$cfg{exp_data_folder}/$cfg{cell_mask_image_prefix}*>;

if ($opt{debug}) {
    print "Cell mask files found: ", join(" ", @cell_mask_files), "\n";
}

my $matlab_code = "";
foreach my $file_name (@cell_mask_files) {
    my $total_images = &get_image_stack_number($file_name);
    foreach my $i_num (1 .. $total_images) {
        next if grep $i_num == $_, @{ $cfg{exclude_image_nums} };

		my $padded_num = sprintf("%0" . length($total_images) . "d", $i_num);
        my $output_path = "$cfg{exp_result_single_folder}/$padded_num";
        mkpath($output_path);
        $matlab_code = $matlab_code . "create_cell_mask_image('$file_name',$i_num,'$output_path')\n";
    }
}

my $error_folder = "$cfg{exp_result_folder}/$cfg{matlab_errors_folder}";
if (not($matlab_wrapper->execute($matlab_code))) {
    mkpath($error_folder);
    open ERR_OUT, ">$error_folder/$cfg{cell_mask_errors_filename}";
    print ERR_OUT $matlab_wrapper->err_msg;
	print ERR_OUT "\n\nMATLAB COMMANDS\n\n$matlab_code";
    print $matlab_wrapper->err_msg if $opt{debug};
    close ERR_OUT;
	
	$matlab_wrapper->remove_files;
}

###############################################################################
# Functions
###############################################################################

#######################################
# Config Collection
#######################################
sub get_config {
    my %cfg = ParseConfig(
        -ConfigFile            => $opt{cfg},
        -MergeDuplicateOptions => 1,
    );

    #check to see if a file for the frames that should be excluded from the
    #analysis is included, if it is, collect the data from it, otherwise, set
    #exclude_image_nums to 0
    if (defined $cfg{exclude_file}) {
        open EX_INPUT, $cfg{exclude_file} or die "Can't open the specified exclude file: $cfg{exclude_file}";
        my $temp_line = <EX_INPUT>;
        close EX_INPUT;

        chomp($temp_line);
        @{ $cfg{exclude_image_nums} } = split(",", $temp_line);
    } else {
        $cfg{exclude_image_nums} = 0;
    }
    if ($opt{debug}) {
        print "Image numbers to be excluded:", join(", ", @{ $cfg{exclude_image_nums} }), "\n";
    }

	#Check for the set of variables in the config file that must be specified
    if (not defined $cfg{results_folder}) {
        die "ERROR: The location of the folder that contains the results of the focal ",
          "adhesion identification must be specified in the config file with the ",
		  "variable \"results_folder\".";
    }

	if (not defined $cfg{exp_name}) {
		die "ERROR: The name of the experiment must be specified in the config ",
		  "file with the variable \"exp_name\".";
	}

    if (not defined $cfg{single_image_folder}) {
        die "ERROR: The location in the results folder where each of the individual ",
          "images and the associated data about those images must be specified in ",
          "the config file with the variable \"single_image_folder\".";
    }

	#Compute a few config variables from the provided values:
    $cfg{exp_data_folder} 	  	   = "$cfg{data_folder}/$cfg{exp_name}";
	$cfg{exp_result_folder}   	   = "$cfg{results_folder}/$cfg{exp_name}";
	$cfg{exp_result_single_folder} = "$cfg{exp_result_folder}/$cfg{single_image_folder}";
	
	return %cfg;
}

#######################################
# Other
#######################################

sub get_image_stack_number {
    my $image_file = pop;

    my $image_info = new Image::ExifTool;
    $image_info->ExtractInfo($image_file);
    my @tag_list    = $image_info->GetFoundTags($image_file);
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
