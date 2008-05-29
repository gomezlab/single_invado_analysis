#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################

use strict;
use File::Path;
use File::Basename;
use File::Find;
use Image::ExifTool;
use Math::Matlab::Local;
use Getopt::Long;
use Data::Dumper;
use File::Spec::Functions;
use Benchmark;

use lib "../lib";
use Config::Adhesions qw(ParseConfig);
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

my %cfg = ParseConfig(\%opt);

###############################################################################
#Main Program
###############################################################################

my $movie_debug_string = $opt{movie_debug} ? ",'debug',1" : "";

our @files;
find(\&include_in_vis, (catdir($cfg{exp_results_folder}, $cfg{tracking_folder})));

my @movie_params = map { 
    my $base_dir = catdir($cfg{exp_results_folder}, $cfg{tracking_folder});
    my $movie_path = $_;
    $movie_path =~ s/$base_dir(.*)\.csv/$1/;
    {
        tracking_file => $_, 
        movie_path => catdir($cfg{movie_output_folder},$movie_path),
    }
} @files;

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

sub include_in_vis {
    if (   $File::Find::name =~ /.csv/ 
        && not($File::Find::name =~ /no_movie/) ) {
       push @files, $File::Find::name;
    }
}

sub build_matlab_visualization_config {
    my %params = @_;

    if (not exists $params{'tracking_file'}) {
        $params{'tracking_file'} = catdir($cfg{exp_results_folder}, $cfg{tracking_folder}, $cfg{tracking_output_file});
    }
    
    my $excluded_image_nums = "["  . join(",", @{$cfg{exclude_image_nums}}) . "]";

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

        "excluded_image_nums = $excluded_image_nums;\n",
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
