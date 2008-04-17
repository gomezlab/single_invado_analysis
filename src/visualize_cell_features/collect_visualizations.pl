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
use Data::Dumper;
use File::Spec::Functions;
use Benchmark;

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
GetOptions(\%opt, "cfg=s", "debug|d", "movie_debug");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my @needed_vars = qw(data_folder results_folder single_image_folder folder_divider exp_name
  single_image_folder matlab_errors_folder vis_config_file vis_errors_file
  extr_val_file bounding_box_file);
my $ad_conf = new Config::Adhesions(\%opt, \@needed_vars);
my %cfg = $ad_conf->get_cfg_hash;

###############################################################################
#Main Program
###############################################################################

my $movie_debug_string = $opt{movie_debug} ? ",'debug',1" : "";

my @movie_folders = split(/\s/, $cfg{movie_output_folders});

#set the first parameter set to be empty to use all the defaults
my @movie_params = (
    { movie_path => $movie_folders[0], },
    {
        tracking_file => catfile($cfg{exp_results_folder}, $cfg{tracking_folder}, 'filtered', 'longevity', '5.csv'),
        movie_path    => $movie_folders[1],
    },
    {
        tracking_file => catfile($cfg{exp_results_folder}, $cfg{tracking_folder}, 'filtered', 'dead', 'dead.csv'),
        movie_path    => $movie_folders[2],
    },
#    {
#        tracking_file =>
#          catfile($cfg{exp_results_folder}, $cfg{tracking_folder}, 'filtered', 'special', 'high_speed.csv'),
#        movie_path  => 'movies/special',
#    },
);

my @matlab_code;
foreach (@movie_params) {
    my %params = %{$_};

    if (not exists $params{'config_file'}) {
        $params{'config_file'} = catfile($cfg{exp_results_folder}, $params{movie_path}, $cfg{vis_config_file});
    }

    mkpath(dirname($params{'config_file'}));

    &write_matlab_config(%params);
    my $error_file = catdir($cfg{exp_results_folder}, $cfg{matlab_errors_folder}, $cfg{vis_errors_file});
    my @matlab_code = "make_movie_frames('" . $params{'config_file'} . "'$movie_debug_string)";

    my $t1 = new Benchmark;
    &Math::Matlab::Extra::execute_commands(\@matlab_code, $error_file);
    my $t2 = new Benchmark;
    print "Movie: $params{movie_path}\n", timestr(timediff($t2, $t1), "nop"), "\n" if $opt{debug};
}

###############################################################################
#Functions
###############################################################################

sub build_matlab_visualization_config {
    my %params = @_;

    if (not exists $params{'tracking_file'}) {
        $params{'tracking_file'} = catdir($cfg{exp_results_folder}, $cfg{tracking_folder}, $cfg{tracking_output_file});
    }

    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime time;
    my @timestamp = join("/", ($mon + 1, $day, $year + 1900)) . " $hour:$min";

    my @temp  = &Image::Data::Collection::gather_sorted_image_numbers(\%cfg);
    my $i_num = scalar(@temp);

    my @config_lines = (
        "%Config file produced by collect_visualizations.pl\n",
        "%@timestamp\n\n",
        "%General Parameters\n",
        "exp_name = '$cfg{exp_name}';\n",
        "base_results_folder = fullfile('", join("\',\'", split($cfg{folder_divider}, $cfg{results_folder})),
        "', exp_name);\n\n",
        "base_data_folder = fullfile('", join("\',\'", split($cfg{folder_divider}, $cfg{data_folder})),
        "', exp_name);\n\n",

        "i_count = $i_num;\n\n",

        "I_folder = fullfile(base_results_folder, '$cfg{single_image_folder}');\n\n",

        "focal_image = 'focal_image.png';\n",
        "adhesions_filename = 'adhesions.png';\n",
        "edge_filename = 'cell_mask.png';\n",

        "tracking_seq_file = '$params{tracking_file}';\n\n",

        "out_path = fullfile(base_results_folder,'$params{movie_path}');\n",
        "out_prefix = {'", join("\',\'", split(/\s/, $cfg{movie_output_prefix})), "'};\n\n",

        "excluded_frames_file = fullfile(base_data_folder,'$cfg{exclude_file}');\n",
        "bounding_box_file = fullfile(base_results_folder,'$params{movie_path}','$cfg{bounding_box_file}');\n",
        "path_folders = '$cfg{path_folders}';\n\n",

        "image_padding_min = $cfg{image_padding_min};\n\n",
    );

    if (exists($cfg{pixel_size})) {
        push @config_lines, "pixel_size = $cfg{pixel_size};\n";
    }

    return @config_lines;
}

sub write_matlab_config {
    my %params = @_;

    my @config = &build_matlab_visualization_config(@_);
    open VIS_CFG_OUT, ">" . $params{'config_file'}
      or die "Unsuccessfully tried to open visualization config file: $params{config_file}";
    print VIS_CFG_OUT @config;
    close VIS_CFG_OUT;
}
