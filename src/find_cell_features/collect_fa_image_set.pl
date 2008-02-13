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
  qw(data_folder results_folder exp_name single_image_folder matlab_errors_folder adhesion_errors_file);
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

mkpath($cfg{individual_results_folder});

my @image_files = <$cfg{exp_data_folder}/$cfg{adhesion_image_prefix}*>;

if ($opt{debug}) {
    if (scalar(@image_files) > 1) {
        print "Adehsion image files found: $image_files[0] - $image_files[$#image_files]\n";
    } else {
        print "Adehsion image file found: $image_files[0]\n";
    }
}

my @matlab_code = &create_matlab_code;
die "@matlab_code";
my $error_file = catdir($cfg{exp_results_folder},$cfg{matlab_errors_folder},$cfg{adhesion_errors_file});

&Math::Matlab::Extra::execute_commands($matlab_wrapper,\@matlab_code,$error_file);

###############################################################################
#Functions
###############################################################################

sub create_matlab_code {
    my @image_stack_count = map { Image::Stack::get_image_stack_number($_) } @image_files;
    
    my @matlab_code;
    if (grep {$_ > 1} @image_stack_count) {
        if (scalar(@image_files) > 1) {
            die "Found more than one image stack in: ", join(", ",@image_files), "\n", 
                "Expected single image stack or multiple non-stacked files\n";
        }
        @matlab_code = &create_matlab_code_stack;
    } else {
        @matlab_code = &create_matlab_code_single;
    }
    return @matlab_code;
}

sub create_matlab_code_stack {
    my @matlab_code;
    
    my $total_stack_images = Image::Stack::get_image_stack_number($image_files[0]);
    foreach my $i_num (1 .. $total_stack_images) {
        next if grep $i_num == $_, @{ $cfg{exclude_image_nums} };
            
        my $padded_num = sprintf("%0" . length($total_stack_images) . "d", $i_num);

        my $output_path = catdir($cfg{individual_results_folder},$padded_num);
        mkpath($output_path);
        $matlab_code[0] .= "find_focal_adhesions('$image_files[0]','I_num',$i_num,'out_dir','$output_path')\n";
    }
    return @matlab_code;
}

sub create_matlab_code_single {
    my @matlab_code;

    foreach my $file_name (@image_files) {
        my $i_num;
        if ($file_name =~ /$cfg{adhesion_image_prefix}(\d+)\./) {
            $i_num = $1;
        } else {
            warn "Unable to find image number in: $file_name, skipping this image.";
            next;
        }

        next if grep $i_num == $_, @{ $cfg{exclude_image_nums} };
        
        my $padded_num = sprintf("%0" . length(scalar(@image_files)) . "d", $i_num);

        my $output_path = catdir($cfg{individual_results_folder},$padded_num);
        mkpath($output_path);
        $matlab_code[0] .= "find_focal_adhesions('$file_name','out_dir','$output_path')\n";
    }
    return @matlab_code;
}
