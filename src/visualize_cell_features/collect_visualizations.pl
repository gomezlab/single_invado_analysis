#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################

use strict;
use File::Path;
use File::Basename;
use Image::ExifTool;
use Math::Matlab::Local;
use Getopt::Long;
use File::Spec::Functions;

use lib "../lib";
use Config::Adhesions;
use Image::Stack;
use Math::Matlab::Extra;
use Image::Data::Collection;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg=s", "debug");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my @needed_vars =
  qw(data_folder results_folder single_image_folder folder_divider exp_name 
     single_image_folder matlab_errors_folder vis_config_file vis_errors_file 
     extr_val_file bounding_box_file);
my $ad_conf = new Config::Adhesions(\%opt, \@needed_vars);
my %cfg = $ad_conf->get_cfg_hash;

my $matlab_wrapper;
if (defined $cfg{matlab_executable}) {
    $matlab_wrapper = Math::Matlab::Local->new({ cmd => "$cfg{matlab_executable} -nodisplay -nojvm -nosplash", });
} else {
    $matlab_wrapper = Math::Matlab::Local->new();
}

###############################################################################
#Main Program
###############################################################################

&write_matlab_config;

my @matlab_code;
$matlab_code[0] .= "make_movie_frames('" . catfile($cfg{exp_data_folder}, $cfg{vis_config_file}) . "')";

my $error_file = catdir($cfg{exp_results_folder}, $cfg{matlab_errors_folder}, $cfg{vis_errors_file});
&Math::Matlab::Extra::execute_commands($matlab_wrapper,\@matlab_code,$error_file);

###############################################################################
#Functions
###############################################################################

sub build_matlab_visualization_config {
    my $adhesion_image = basename <$cfg{exp_data_folder}/$cfg{adhesion_image_prefix}*>;
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime time;
    my @timestamp = join("/", ($mon + 1, $day, $year + 1900)) . " $hour:$min";
    
    my @temp = &Image::Data::Collection::gather_sorted_image_numbers(\%cfg);
    my $i_num = scalar(@temp);

    my @config_lines = (
        "%Config file produced by collect_visualizations.pl\n",
        "%@timestamp\n\n",
        "%General Parameters\n",
        "exp_name = '$cfg{exp_name}'\n",
        "base_data_folder = fullfile('", join("\',\'", split($cfg{folder_divider}, $cfg{data_folder})),
        "', exp_name)\n",
        "base_results_folder = fullfile('", join("\',\'", split($cfg{folder_divider}, $cfg{results_folder})),
        "', exp_name)\n\n",

        "i_count = $i_num\n\n",

        "I_folder = fullfile(base_results_folder, '$cfg{single_image_folder}')\n\n",

        "adhesions_filename = 'adhesions.png'\n",
        "edge_filename = 'cell_mask.png'\n",

        "tracking_seq_file = fullfile(base_results_folder,'$cfg{tracking_output_file}')\n\n",

        "out_path = fullfile(base_results_folder,'$cfg{movie_output_folder}')\n",
        "out_prefix = {'", join("\',\'", split(/\s/, $cfg{movie_output_prefix})), "'}\n\n",

        "excluded_frames_file = fullfile(base_data_folder,'$cfg{exclude_file}')\n",
        "extr_val_file = fullfile(base_data_folder,'$cfg{extr_val_file}')\n",
        "bounding_box_file = fullfile(base_data_folder,'$cfg{bounding_box_file}')\n",
        "path_folders = '$cfg{path_folders}'\n\n",

        "image_padding_min = $cfg{image_padding_min}\n\n",

        "%Comparison Specific Settings\n\n",
    );
    return @config_lines;
}

sub write_matlab_config {
    my @config = &build_matlab_visualization_config;
    open VIS_CFG_OUT, ">" . catfile($cfg{exp_data_folder}, $cfg{vis_config_file})
      or die "Unsuccessfully tried to open visualization config file: "
      . catfile(qw($cfg{exp_data_folder} $cfg{vis_config_file}));
    print VIS_CFG_OUT @config;
    close VIS_CFG_OUT;
}
