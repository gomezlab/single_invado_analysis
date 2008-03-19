#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use warnings;
use File::Spec;

use base qw(Config::General);

my %derived_vars = (
    individual_results_folder => [qw(results_folder exp_name single_image_folder)],
    exp_results_folder        => [qw(results_folder exp_name)],
    exp_data_folder           => [qw(data_folder exp_name)],
);

###############################################################################
# Module Definition
###############################################################################
package Config::Adhesions;

sub new {
    my $class        = $_[0];
    my %opt          = %{ $_[1] };
    my @needed_files = @{ $_[2] };

    my $conf = new Config::General(
        -ConfigFile            => $opt{cfg},
        -MergeDuplicateOptions => 1,
        -IncludeRelative       => 1,
    );
    my %cfg = $conf->getall;
    %{ $cfg{opt} }          = %opt;
    @{ $cfg{needed_files} } = @needed_files;
    %{ $cfg{var_desc} }     = $class->var_descriptions;

    my $cfg_ref = \%cfg;

    bless $cfg_ref, $class;

    $cfg_ref->check_for_missing_vars;
    $cfg_ref->build_derived_parameters;
    $cfg_ref->collect_cfg_info_from_files;

    return $cfg_ref;
}

sub collect_cfg_info_from_files {
    my $cfg = shift;

    #check to see if a file for the frames that should be excluded from the
    #analysis is included, if it is, collect the data from it, otherwise, set
    #exclude_image_nums to 0

    if (defined $cfg->{exclude_file}) {
        my $exclude_file = File::Spec->catfile($cfg->{data_folder}, $cfg->{exp_name}, $cfg->{exclude_file});
        open EX_INPUT, $exclude_file or die "Can't open the specified exclude file: $exclude_file";
        my $temp_line = <EX_INPUT>;
        close EX_INPUT;

        if (not($temp_line)) {
            @{ $cfg->{exclude_image_nums} } = (0);
        } else {
            chomp($temp_line);
            @{ $cfg->{exclude_image_nums} } = split(",", $temp_line);
        }
    } else {
        @{ $cfg->{exclude_image_nums} } = (0);
    }
}

sub build_derived_parameters {
    my $cfg = $_[0];
    foreach my $this_key (keys %derived_vars) {
        my $all_present = 1;
        foreach my $var_name (@{ $derived_vars{$this_key} }) {
            $all_present = 0 if not defined $cfg->{$var_name};
        }

        if ($all_present) {
            $cfg->{$this_key} = File::Spec->catdir(map $cfg->{$_}, @{ $derived_vars{$this_key} });
        }
    }
}

sub check_for_missing_vars {
    my $cfg           = shift;
    my $missing_files = 0;
    foreach (@{ $cfg->{needed_files} }) {
        if (not defined $cfg->{$_}) {
            warn "Missing variable(s) in config file ($cfg->{opt}->{cfg}):\n\n" if not($missing_files);
            warn "$_: $cfg->{var_desc}->{$_}\n";
            $missing_files = 1;
        }
    }
    die "Please fill in values for the above variables and rerun the program." if $missing_files;
    delete $cfg->{var_desc};
}

sub var_descriptions {
    my %desc;
    ########################################################################
    # General Parameters
    ########################################################################

    $desc{data_folder} =
      "Specifies where the data folder is located. This folder should contain folders for individual experiments.";

    $desc{results_folder} = "Specifies where the general results folder is located.";

    $desc{exp_name} = "Specifies the folder where the experimenal data is located.";

    $desc{single_image_folder} =
      "Specifies the subfolder in results_folder which contains the results from each single image.";

    $desc{raw_data_folder} =
      "Specifies the subfolder in each single_image_folder where the parameters extracted from MATLAB are stored.";

    $desc{exclude_file} =
      "Specifies the file name of the data file in the data_folder which contains the image numbers that are to be excluded in analysis.";

    ########################################################################
    # Feature Collection Parameters
    ########################################################################

    $desc{matlab_executable} = "Specifies the location of the MATLAB executable.";

    $desc{cell_mask_image_prefix} = "Specifies the prefix of the file that will be used to find the cell mask.";

    $desc{adhesion_image_prefix} = "Specifies the prefix of the file that will be used to find the adhesions.";

    $desc{matlab_errors_folder} = "Specifies the folder where a record of the MATLAB errors will be stored.";

    $desc{cell_mask_errors_filename} =
      "Specifies the file name to use when storing the errors from the cell mask finding";

    $desc{adhesion_errors_filename} =
      "Specifies the file name to use when storing the errors from the adhesion finding.";

    ###############################################################################
    # Plotting Parameters
    ###############################################################################

    $desc{file_ext} = "Specifies the kind of output plots desired for each plot created.";

    $desc{pixel_size_file} = "Specifies the file in data_folder which contains the size of the pixels in meters.";

    $desc{target_unit_size} = "Specifies the unit size that should be used on the plots.";

    ###############################################################################
    # Tracking/Lineage Analysis Parameters
    ###############################################################################

    #######################################
    # General Parameters
    #######################################

    $desc{general_data_files} =
      "Specifies the files in the raw_data_folder that are needed for both the tracking and analysis of the adhesions.";

    #######################################
    # Tracking Parameters
    #######################################

    $desc{tracking_files} = "Specifies the files in the raw_data_folder that are needed specifically for tracking.";

    $desc{tracking_probs_folder} =
      "Specifies the folder under results_folder where the tracking problems will be stored.";

    $desc{tracking_output_file} =
      "Specifies the file name under the results_folder that will be used to store the tracking matrix.";

    ####################
    # Live Tracking Parameters
    ####################

    $desc{pix_sim_indeter_percent} =
      "Specifies the percentage of first pixel similarity measure the second pixel similarity measure must be to trigger live adhesion tracking by the closest centroid.";

    ####################
    # Merge Parameters
    ####################

    $desc{merge_shift_percent} =
      "Specifies the percentage of adhesion size difference that causes a merge decision to be made by the centroid shift of the merging adhesions.";

    #######################################
    # Lineage Analysis
    #######################################

    $desc{lineage_analysis_data_files} =
      "Specifies the files in the raw_data_folder that are needed specifically for analysis.";

    $desc{lineage_props_folder} =
      "Specifies the folder under results_folder/exp_name where the lineage properties will be saved.";

    $desc{single_lineage_props_file} =
      "Specifies the file name to use when outputing the properties of each lineage, csv format";

    ###############################################################################
    #Visualization Settings
    ###############################################################################

    #######################################
    #General Settings
    #######################################

    $desc{vis_config_file} =
      "Specifies the file name to use for the matlab visualization config, file will reside in results_folder/exp_name";

    $desc{vis_errors_file} =
      "Specifies the file name to use for storing any errors encountered during the creation of the visualization, will reside in results_folder/exp_name/matlab_errors_folder";

    $desc{extr_val_file} =
      "Specifies the filename to use when searching for the extreme values in the image data, should be in data_folder/exp_name";

    $desc{bounding_box_file} =
      "Specifies the filename to use when searching for the extreme values in the image data, should be in data_folder/exp_name";

    $desc{path_folders} =
      "Specifies the folder which contains (or contains subfolders which contain) all the matlab programs needed for visualization";

    $desc{image_padding_min} =
      "Specifies the minimum number of pixels between the edge of the features on the extremes of an image and the edge of the image";

    #######################################
    #Movie Frame Settings
    #######################################

    $desc{movie_output_folder} =
      "Specifies the folder where the movie frames will be saved, folder will be located in \"results_folder/exp_name\"";

    $desc{movie_output_prefix} =
      "Specifies the folders where each image written by the visualization code will placed, each folder will be based in \"results_folder/exp_name/movie_output_folder\"";

    return %desc;
}

sub get_cfg_hash {
    my $self = shift;
    return map { $_ => ${$self}{$_} } keys %{$self};
}

1;
