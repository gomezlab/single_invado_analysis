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
        open EX_INPUT, File::Spec->catdir($cfg->{data_folder}, $cfg->{exp_name}, $cfg->{exclude_file})
          or die "Can't open the specified exclude file: ", File::Spec->catdir($cfg->{data_folder}, $cfg->{exp_name}, $cfg->{exclude_file});
        my $temp_line = <EX_INPUT>;
        close EX_INPUT;

        chomp($temp_line);
        @{ $cfg->{exclude_image_nums} } = split(",", $temp_line);
    } else {
        @{ $cfg->{exclude_image_nums} } = (0);
    }
    if ($cfg->{opt}->{debug}) {
        if (${ $cfg->{exclude_image_nums} }[0] == 0) {
            print "Image numbers to be excluded: none\n";
        } else {
            print "Image numbers to be excluded:", join(", ", @{ $cfg->{exclude_image_nums} }), "\n";
        }
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
}

sub var_descriptions {
    my %desc;
    ########################################################################
    # General Parameters
    ########################################################################

    #Experimenal Data Folder
    $desc{data_folder} =
      "Specifies where the data folder is located." . " This folder should contain folders for individual experiments.";

    #Where the results located
    $desc{results_folder} = "Specifies where the general results folder is located.";

    #Experiment Name
    $desc{exp_name} = "Specifies the folder where the experimenal data is located.";

    #Individual Image Folder
    $desc{single_image_folder} =
      "Specifies the subfolder in results_folder which contains the results from each single image.";

    #The folder in each individual image folder with raw data
    $desc{raw_data_folder} =
      "Specifies the subfolder in each single_image_folder where the parameters extracted from MATLAB are stored.";

    #The filename of the file which indicates which frames are to be excluded,
    #should be located in data_folder/exclude_file
    $desc{exclude_file} =
      "Specifies the file name of the data file in the data_folder which contains the image numbers that are to be excluded in analysis.";

    ########################################################################
    # Feature Collection Parameters
    ########################################################################

    #MATLAB executable to use
    $desc{matlab_executable} = "Specifies the location of the MATLAB executable.";

    #Output Folder, if absent, use experiment name
    $desc{feature_output_folder} =
      "If present, value is used as the subfolder under results_folder where the feature data is stored, otherwise, exp_name is used. ";

    #Prefix of the stack file to use for identifying the cell mask
    $desc{cell_mask_image_prefix} = "Specifies the prefix of the file that will be used to find the cell mask.";

    #Prefix of the stack file to use for identifying the focal adhesions
    $desc{adhesion_image_prefix} = "Specifies the prefix of the file that will be used to find the adhesions.";

    #Folder to use for MATLAB errors
    $desc{matlab_errors_folder} = "Specifies the folder where a record of the MATLAB errors will be stored.";

    #Filename for MATLAB errors from cell mask collection
    $desc{cell_mask_errors_filename} =
      "Specifies the file name to use when storing the errors from the cell mask finding";

    #Filename for MATLAB errors from cell mask collection
    $desc{adhesion_errors_filename} =
      "Specifies the file name to use when storing the errors from the adhesion finding.";

    ###############################################################################
    # Plotting Parameters
    ###############################################################################

    #File format of output plots
    $desc{file_ext} = "Specifies the kind of output plots desired for each plot created.";

    #Where the pixel size file can be found
    $desc{pixel_size_file} = "Specifies the file in data_folder which contains the size of the pixels in meters.";

    #Target units for measurements, e.g.
    #	-assume that the pixels are each 100 nm (100E-9 m) on each edge
    #		-set through the "pixel_size_file" file
    #	-assume that the target units is 1 micro meter (1E-6 m)
    #		-set through the "target_unit_size"
    #	-Centroid_distance_from_edge = 140 pixels
    #	 140 pix*(100E-9 m/1 pix)*(1 micro meter/1E-6 m) = 14 micro meters
    $desc{target_unit_size} = "Specifies the unit size that should be used on the plots.";

    ###############################################################################
    # Tracking/Lineage Analysis Parameters
    ###############################################################################

    #######################################
    # General Parameters
    #######################################

    #List of data files that will be needed in both the tracking and analysis
    $desc{general_data_files} =
      "Specifies the files in the raw_data_folder that are needed for both the tracking and analysis of the adhesions.";

    #######################################
    # Tracking Parameters
    #######################################

    #Data files needed just for tracking
    $desc{tracking_files} = "Specifies the files in the raw_data_folder that are needed specifically for tracking.";

    #Folder name to be used for output of tracking problems, will become a
    #subfolder under the folder identified in results_folder
    $desc{tracking_probs_folder} =
      "Specifies the folder under results_folder where the tracking problems will be stored.";

    #Filename to use for tracking matrix output, will be placed in the the folder
    #identified in results_folder
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

    #allowable percent difference in adheison size to use centroid shift as merge
    #determining factor, defaults to 0.25
    $desc{merge_shift_percent} =
      "Specifies the percentage of adhesion size difference that causes a merge decision to be made by the centroid shift of the merging adhesions.";

    #######################################
    # Lineage Analysis
    #######################################

    #There are data files that are needed for the lineage analysis that aren't
    #included in the general files list, they should be included in this list
    $desc{lineage_analysis_data_files} =
      "Specifies the files in the raw_data_folder that are needed specifically for analysis.";

    #Folder to use to hold the analysis of the adhesion lineages, will become a
    #subfolder under the folder identified in results_folder
    $desc{lineage_props_folder} =
      "Specifies the folder under results_folder/exp_name where the lineage properties will be saved.";

    return %desc;
}

sub get_cfg_hash {
    my $self = shift;
    return map { $_ => ${$self}{$_} } keys %{$self};
}

1;
