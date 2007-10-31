#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################

use strict;
use File::Path;
use Image::ExifTool;
use Math::Matlab::Local;
use Config::General;
use Getopt::Long;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d");
die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my %cfg = Config::General::ParseConfig($opt{cfg});

my $matlab_wrapper = Math::Matlab::Local->new({ cmd => '/usr/bin/matlab -nodisplay -nojvm -nosplash', });

###############################################################################
# Main Program
###############################################################################

my $path = $cfg{results_folder} . "/" . $cfg{exp_name} . "/individual_pictures";

mkpath($path);

my $adhesion_image_location = $cfg{exp_data_folder} . "/" . $cfg{exp_name} . "/" . $cfg{adhesion_image_prefix};
my @adhesion_image_files    = <$adhesion_image_location*>;

if ($opt{debug}) {
    print "Adhesion image files found: ", join(" ", @adhesion_image_files), "\n";
}

foreach my $file_name (@adhesion_image_files) {
    my $total_images = &get_image_stack_number($file_name);
    foreach my $image_num (1 .. $total_images) {
        if ($opt{debug} && $image_num > 1) {
            next;
        }
        if ($opt{debug} && $image_num % 10 == 0) {
            print "Working on image number $image_num of $total_images\n";
        }

        my $padded_num = sprintf("%0" . length($total_images) . "d", $image_num);
        my $output_path = $cfg{results_folder} . "/" . $cfg{exp_name} . "/individual_pictures/$padded_num";
        mkpath($output_path);
        my $matlab_code = "find_focal_adhesions('$file_name','I_num',$image_num,'out_dir','$output_path')\n";
        if (not($matlab_wrapper->execute($matlab_code))) {
            my $error_folder = $cfg{results_folder} . "/" . $cfg{exp_name} . "/errors/";
            mkpath($error_folder);
            open ERR_OUT, ">$error_folder" . "fa_error_$image_num.txt";
            print ERR_OUT $matlab_wrapper->err_msg;
            close ERR_OUT;
        }
    }
}

###############################################################################
# Functions
###############################################################################

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
