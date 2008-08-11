#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################
use lib "../lib";
use lib "../lib/perl";

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

use Config::Adhesions qw(ParseConfig);
use Image::Stack;
use Math::Matlab::Extra;
use Emerald;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg=s", "debug|d", "movie_debug", "config_only", "emerald|e");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my %cfg = ParseConfig(\%opt);

###############################################################################
#Main Program
###############################################################################
my $movie_debug_string = $opt{movie_debug} ? ",'debug',1" : "";

our @files;
find(\&include_in_vis, catdir($cfg{exp_results_folder}, $cfg{tracking_folder}));

my @movie_params = map {
    my $base_dir = catdir($cfg{exp_results_folder}, $cfg{tracking_folder});
    my $movie_path = $_;
    $movie_path =~ s/(.*)\.csv/$1/;
    {
        tracking_file => $_,
        movie_path    => catdir($cfg{movie_output_folder}, $movie_path),
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
    
    push @matlab_code, "make_movie_frames('" . $params{'config_file'} . "'$movie_debug_string)";
}

my $error_folder = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'visualization');
my $error_file = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'visualization', 'error.txt');
mkpath($error_folder);

if ($opt{emerald}) {
    my %emerald_opt = ("folder", $error_folder);
    my @commands = &Emerald::create_LSF_Matlab_commands(\@matlab_code, \%emerald_opt);
    &Emerald::send_LSF_commands(\@commands);
} elsif (not $opt{config_only}) {
    my $t1 = new Benchmark;
    &Math::Matlab::Extra::execute_commands(\@matlab_code, $error_file);
    my $t2 = new Benchmark;
    print timestr(timediff($t2, $t1), "nop"), "\n" if $opt{debug};
}

###############################################################################
#Functions
###############################################################################
sub include_in_vis {
    if ($File::Find::name =~ /.csv/
        && not($File::Find::name =~ /no_movie/)) {
        my $folder = catdir($cfg{exp_results_folder}, $cfg{tracking_folder});
        if ($File::Find::name =~ /$folder\/(.*)/) {
            push @files, $1;
        }
    }
}

sub write_matlab_config {
    my %params = @_;

    my @config = &build_matlab_visualization_config(@_);
    open VIS_CFG_OUT, ">" . $params{'config_file'}
      or die "Unsuccessfully tried to open visualization config file: $params{config_file}";
    print VIS_CFG_OUT @config;
    close VIS_CFG_OUT;
}

sub build_matlab_visualization_config {
    my %params = @_;

    if (not exists $params{'tracking_file'}) {
        $params{'tracking_file'} = catdir($cfg{exp_results_folder}, $cfg{tracking_folder}, $cfg{tracking_output_file});
    }

    my $excluded_image_nums = "[" . join(",", @{ $cfg{exclude_image_nums} }) . "]";

    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime time;
    my @timestamp = join("/", ($mon + 1, $day, $year + 1900)) . " $hour:$min";

    my @config_lines = (
        "%Config file produced by collect_visualizations.pl\n",
        "%@timestamp\n\n",
        "%General Parameters\n",
        "exp_name = '$cfg{exp_name}';\n",
        "base_results_folder = fullfile('", join("\',\'", split($cfg{folder_divider}, $cfg{results_folder})),
        "', exp_name);\n\n",

        "I_folder = fullfile(base_results_folder, '$cfg{single_image_folder}');\n\n",

        "focal_image = '$cfg{adhesion_image_file}';\n",
        "adhesions_filename = 'adhesions.png';\n",
        "edge_filename = '$cfg{cell_mask_file}';\n",

        "tracking_seq_file = fullfile(base_results_folder, '$cfg{tracking_folder}', '$params{tracking_file}');\n\n",

        "out_path = fullfile(base_results_folder,'$params{movie_path}');\n",
        "out_prefix = {'", join("\',\'", @{ $cfg{movie_output_prefix} }), "'};\n\n",

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

###############################################################################
#Documentation
###############################################################################

=head1 NAME

collect_visualizations.pl - build the visualizations of the focal adhesion
movies

=head1 SYNOPSIS

collect_mask_set.pl -cfg FA_config

=head1 DESCRIPTION

This program builds a series of movies based on files available in the tracking
matrices folder. Each file in the tracking matrices folder which ends with
'.csv' and does not contain 'no_movie' is used to build a visualization of the
tracked focal adhesions.

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=back

Optional parameter(s):

=over 

=item * debug or d: print debuging information during program execution

=item * movie_debug: pass along the debug flag to the MATLAB visualization
program, causing only a small subset of the tracked adhesions to be visualized
in a single frame

=item * config_only: only write the MATLAB config files out, do not execute the
MATLAB program

=item * emerald: build and execute long commands throught the LSF job system

=back

=head1 EXAMPLES

collect_visualizations.pl -cfg FA_config

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 6/30/2008

=cut
