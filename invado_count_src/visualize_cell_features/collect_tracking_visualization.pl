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
use POSIX;

use Config::Adhesions qw(ParseConfig);
use Image::Stack;
use Math::Matlab::Extra;
use Emerald;
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg=s", "debug|d", "movie_debug", "config_only|only_config", 
                  "lsf|l", "single_ad_folders") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my %cfg = ParseConfig(\%opt);

###############################################################################
#Main Program
###############################################################################

my $config_file = catfile($cfg{exp_results_folder}, $cfg{movie_output_folder},'vis_config.m');
mkpath(dirname($config_file));

&write_matlab_config($config_file);

my @matlab_code = "make_movie_frames('$config_file')";

# push @matlab_code, &build_single_ad_commands(%single_ad_params);
# if ($opt{single_ad_folders}) {
#     push @matlab_code, &build_single_ad_folder_command(%single_ad_params);
# }

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'visualization');
$opt{error_file} = catfile($opt{errors_folder}, 'error.txt');
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

if (not($opt{config_only})) {
    &FA_job::run_matlab_progam(\@matlab_code,\%opt);
}

###############################################################################
#Functions
###############################################################################
sub write_matlab_config {
	my $config_file = shift;

    my @config = &build_matlab_visualization_config(@_);
    open VIS_CFG_OUT, ">" . $config_file
      or die "Unsuccessfully tried to open visualization config file: $config_file";
    print VIS_CFG_OUT join("\n",@config);
    close VIS_CFG_OUT;
}

sub build_matlab_visualization_config {
    my %params = @_;

    my $excluded_image_nums = "[" . join(",", @{ $cfg{exclude_image_nums} }) . "]";

    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime time;
    my $timestamp = join("/", ($mon + 1, $day, $year + 1900)) . " $hour:$min";

    my @config_lines = (
        "%Config file produced by $0",
        "%$timestamp\n",
        "%General Parameters",
        "exp_name = '$cfg{exp_name}';",
        "base_results_folder = fullfile('" . join("\',\'", split($cfg{folder_divider}, $cfg{results_folder})) .  "', exp_name);\n",
        "vis_folder = fullfile(base_results_folder,'$cfg{movie_output_folder}');\n",

        "I_folder = fullfile(base_results_folder, '$cfg{single_image_folder}');\n",
        
        "lin_time_series_folder = fullfile(base_results_folder, '$cfg{adhesion_props_folder}', '$cfg{lineage_ts_folder}');\n",

        # "puncta_image_filename = 'registered_focal_image.png';",
        # "adhesions_filename = 'puncta_labeled.png';",
        # "adhesions_perim_filename = 'puncta_labeled_perim.png';",
        # "edge_filename = 'cell_mask.png';",
        # "gel_image_filename = 'registered_gel.png';",

        "tracking_seq_file = fullfile(base_results_folder, '$cfg{tracking_folder}', '$cfg{tracking_output_file}');\n",

        "out_path = fullfile(vis_folder,'tracking');",
        "out_prefix = {'" . join("\',\'", @{ $cfg{movie_output_prefix} }) . "'};\n",
        
        "out_path_single = fullfile(vis_folder,'single_ad');",
    
        "path_folders = '$cfg{path_folders}';\n",

        "image_padding_min = $cfg{padding_min};",
        "single_image_padding_min = $cfg{single_ad_padding_min};\n",
    );
    
    if (exists($cfg{pixel_size})) {
        push @config_lines, "pixel_size = $cfg{pixel_size};\n";
    }

    return @config_lines;
}

sub build_single_ad_folder_command {
    my %single_ad_params = @_;

    my @commands;
    push @commands, "make_single_ad_folders('" . $single_ad_params{'config_file'} . "')";
    return @commands;
}

sub build_single_ad_commands {
    my %single_ad_params = @_;

    open TRACKING_FILE, catfile($cfg{exp_results_folder}, $cfg{tracking_folder}, $cfg{tracking_output_file}) or die;
    my @tracking_file = <TRACKING_FILE>;
    my $line_count = scalar(@tracking_file);
    close TRACKING_FILE;

    my $ad_per_run = 5000;
    my @commands;
    my $assembly_file = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, 'invado_data.csv');
    if (-e $assembly_file) {
        my $assembly_option = "'adhesion_file','$assembly_file'";
        for (0 .. (ceil($line_count/$ad_per_run) - 1)) {
            my $start_row = $_ * $ad_per_run + 1;
            my $end_row = $start_row + $ad_per_run - 1;
            $end_row = $line_count if $end_row > $line_count;
            push @commands, "make_single_ad_frames('" . $single_ad_params{'config_file'} . "','start_row',$start_row,'end_row',$end_row,$assembly_option)";
        }
    }
    
    return @commands;
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
