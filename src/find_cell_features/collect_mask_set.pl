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

#Perl built-in variable that controls buffering print output, 1 turns off 
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg=s", "debug");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my @needed_vars =
  qw(data_folder results_folder exp_name single_image_folder matlab_errors_folder cell_mask_errors_filename);
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

my @cell_mask_files = <$cfg{exp_data_folder}/$cfg{cell_mask_image_prefix}*>;

if ($opt{debug}) {
    print "Cell mask files found: ", join(" ", @cell_mask_files), "\n";
}

my $matlab_code = "";
foreach my $file_name (@cell_mask_files) {
    my $total_images = Image::Stack::get_image_stack_number($file_name);
    foreach my $i_num (1 .. $total_images) {
        next if grep $i_num == $_, @{ $cfg{exclude_image_nums} };

        my $padded_num = sprintf("%0" . length($total_images) . "d", $i_num);
        my $output_path = catdir($cfg{individual_results_folder},$padded_num);
        mkpath($output_path);
        $matlab_code = $matlab_code . "find_cell_mask('$file_name','I_num',$i_num,'out_dir','$output_path')\n";
    }
}
die;
my $error_folder = catdir($cfg{exp_result_folder},$cfg{matlab_errors_folder});
if (not($matlab_wrapper->execute($matlab_code))) {
    mkpath($error_folder);
    open ERR_OUT, ">$error_folder/$cfg{cell_mask_errors_filename}";
    print ERR_OUT $matlab_wrapper->err_msg;
    print ERR_OUT "\n\nMATLAB COMMANDS\n\n$matlab_code";
    close ERR_OUT;

    print $matlab_wrapper->err_msg if $opt{debug};

    $matlab_wrapper->remove_files;
}
