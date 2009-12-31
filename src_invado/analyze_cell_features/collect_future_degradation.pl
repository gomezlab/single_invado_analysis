#!/usr/bin/env perl

################################################################################
# Global Variables and Modules
################################################################################

use lib "../lib";
use lib "../lib/perl";

use strict;
use File::Path;
use File::Spec::Functions;
use File::Basename;
use Image::ExifTool;
use Getopt::Long;
use Data::Dumper;
use Text::CSV::Simple;
use Statistics::Descriptive;

use Config::Adhesions;
use Math::Matlab::Extra;
use Emerald;
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

################################################################################
# Main Program
################################################################################

########################################
#Collect all lineage numbers to test
########################################
my $parser = Text::CSV::Simple->new;
my @invado_nums = $parser->read_file(catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, 'invado_nums.csv'));
@invado_nums = @{$invado_nums[0]};

my @non_invado_nums = $parser->read_file(catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, 'non_invado_nums.csv'));
@non_invado_nums = @{$non_invado_nums[0]};

my @all_lin_nums = (@invado_nums, @non_invado_nums);
@all_lin_nums = sort {$a <=> $b} @all_lin_nums;

########################################
#Find max area indexes
########################################
my $area_data_file = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{lineage_ts_folder}, 'Area.csv');
my @start_cols = &find_max_area_indexes(\$area_data_file, \@all_lin_nums);

########################################
#Prepare and send off the matlab command
########################################
my $lin_num_matlab = "[" . join(",",@all_lin_nums)  . "]";
my $start_cols_matlab = "[" . join(",",@start_cols)  . "]";

my $tracking_file = catfile($cfg{exp_results_folder}, $cfg{tracking_folder}, $cfg{tracking_output_file});

my $output_file = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{lineage_ts_folder}, 'Local_diff_future.csv');

my @matlab_code = ("find_future_degrade_values('$cfg{individual_results_folder}','$tracking_file',$lin_num_matlab, $start_cols_matlab, '$output_file')");

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'future_degrade');
$opt{error_file} = catfile($opt{error_folder}, 'error.txt');
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

&FA_job::run_matlab_progam(\@matlab_code,\%opt);

################################################################################
# Main Program
################################################################################

sub find_max_area_indexes {
	my $area_data_file = ${$_[0]};
	my @all_lin_nums = @{$_[1]};

	my @area_data = $parser->read_file($area_data_file);

	my @start_cols;
	foreach (@all_lin_nums) {
		my @areas = @{$area_data[$_ - 1]};
		
		my @only_nums_areas = grep {
			if ($_ eq "NaN") {
				0;
			} else {
				1;
			}
		} @areas;
		
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@only_nums_areas);
		
		#search for the index that corresponds to the last position of the highest
		#area index, we also add one to that index because we will pass this to
		#matlab (indexing begins at 0 there)
		my $last_max_area_index = (grep $areas[$_] == $stat->max, (0..$#areas))[-1] + 1;
			
		push @start_cols, $last_max_area_index;
	}

	die "Difference in lengths of lineage numbers and start columns." if ($#all_lin_nums != $#start_cols);
	return @start_cols;
}
