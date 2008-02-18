#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################

use strict;
use File::Path;
use Image::ExifTool;
use Math::Matlab::Local;
use Getopt::Long;
use File::Spec::Functions;

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

my @needed_vars =
  qw(data_folder results_folder exp_name single_image_folder matlab_errors_folder 
     adhesion_image_file adhesion_errors_file);
my $ad_conf = new Config::Adhesions(\%opt, \@needed_vars);
my %cfg = $ad_conf->get_cfg_hash;

my $matlab_wrapper;
if (defined $cfg{matlab_executable}) {
    $matlab_wrapper = Math::Matlab::Local->new({ cmd => "$cfg{matlab_executable} -nodisplay -nojvm -nosplash", });
} else {
    $matlab_wrapper = Math::Matlab::Local->new();
}

###############################################################################
# Main Program
###############################################################################

my @image_files = <$cfg{individual_results_folder}/*/$cfg{adhesion_image_file}>;

if ($opt{debug}) {
    if (scalar(@image_files) > 1) {
        print "Adhesion files found: $image_files[0] - $image_files[$#image_files]\n";
    } else {
        print "Adhesion file found: $image_files[0]\n";
    }
}

my @matlab_code = &create_matlab_code;

my $error_file = catfile($cfg{exp_results_folder},$cfg{matlab_errors_folder},$cfg{adhesion_image_file});
&Math::Matlab::Extra::execute_commands($matlab_wrapper,\@matlab_code,$error_file);

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
            die "Skipping file: $file_name\n",
                "Unable to find image number.";
            next;
        }

        next if grep $i_num == $_, @{ $cfg{exclude_image_nums} };
        
        my $padded_num = sprintf("%0" . length(scalar(@image_files)) . "d", $i_num);

        my $output_dir = dirname($file_name);
        my $cell_mask = catfile(dirname($file_name),"cell_mask.png");

        $matlab_code[0] .= "find_focal_adhesions('$file_name','$cell_mask')\n";
    }
    
    return @matlab_code;
}
