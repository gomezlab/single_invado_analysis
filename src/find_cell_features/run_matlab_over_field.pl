#!/usr/bin/env perl

################################################################################
# Global Variables and Modules
################################################################################

use lib "../lib";

use strict;
use File::Path;
use File::Spec::Functions;
use File::Spec;
use File::Basename;
use Getopt::Long;
use Data::Dumper;

use Config::ImageSet;
use Math::Matlab::Extra;
use Emerald;
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "script=s", "debug|d", "lsf|l", "queue=s",
	"resource|R=s", "try_wrap") or die;

die "Can't find script specified on the command line" if not exists $opt{script};

my @cfg_files = split(",",$opt{cfg});
my @cfg_set;

foreach (@cfg_files) {
	die "Can't find cfg file specified on the command line" if not exists $opt{cfg};
	my $ad_conf = new Config::ImageSet({cfg => $_});
	my %cfg     = $ad_conf->get_cfg_hash;
	push @cfg_set, \%cfg;
}
my %sample_cfg = %{$cfg_set[0]};

if ($opt{script} =~ /(.*)\.m/) {
	$opt{script} = $1;
}

$opt{script_dir} = dirname($opt{script});
chdir $opt{script_dir};
if ($opt{script} =~ /$opt{script_dir}\/(.*)/) {
	$opt{script} = $1;
}
$opt{root_mwd} = File::Spec->rel2abs($opt{script_dir});

################################################################################
# Main Program
################################################################################

my $matlab_code;
for (@cfg_set) {
	my %cfg = %{$_};
	my $extra = &build_extra_command_line_opts(%cfg);
	if ($opt{try_wrap}) {
		#add try wrap to make sure errors in one script don't kill the next command
		$matlab_code = $matlab_code . "try,$opt{script}('$cfg{exp_results_folder}'$extra),end,"
	} else {
		$matlab_code = $matlab_code . "$opt{script}('$cfg{exp_results_folder}'$extra); "
	}
}
my @matlab_code = ($matlab_code);

$opt{error_folder} = catdir($sample_cfg{exp_results_folder}, $sample_cfg{errors_folder}, $opt{script});
$opt{error_file} = catfile($opt{error_folder}, 'error.txt');
if (defined $sample_cfg{job_group}) {
    $opt{job_group} = $sample_cfg{job_group};
}

%opt = &add_extra_options(%opt);

&FA_job::run_matlab_progam(\@matlab_code,\%opt);

################################################################################
# Functions
################################################################################

sub build_extra_command_line_opts {
	my %this_config = @_;
	my $extra = '';
	
	if ($opt{script} eq "apply_bleaching_correction") {
		if (defined $this_config{gel_norm_level}) {
	    	$extra .= ",'gel_norm_level',$this_config{gel_norm_level}";
		}
	}
	
	if ($opt{script} eq "find_puncta_thresh") {
	}

	if ($opt{script} eq "find_puncta") {
        if (defined $this_config{stdev_thresh}) {
			my @split_stdev_vals = split(/\s+/,$this_config{stdev_thresh});
            $extra .= ",'stdev_thresh',[" . join(",",@split_stdev_vals) . "]";
        }
        if (defined $this_config{min_puncta_seed_size}) {
	    	$extra .= ",'min_puncta_seed_size',$this_config{min_puncta_seed_size}";
        }
        if (defined $this_config{min_puncta_size}) {
	    	$extra .= ",'min_puncta_size',$this_config{min_puncta_size}";
	    }
        if (defined $this_config{max_puncta_size}) {
	    	$extra .= ",'max_puncta_size',$this_config{max_puncta_size}";
	    }
        if (defined $this_config{max_ratio}) {
	    	$extra .= ",'max_ratio',$this_config{max_ratio}";
	    }
	}

	if ($opt{script} eq "gather_tracking_results") {
        if (defined $this_config{pixel_size}) {
	    	$extra .= ",'pixel_size',$this_config{pixel_size}";
        }
	}

	return $extra;
}

sub add_extra_options {
	my %opt = @_;

	if ($opt{script} =~ /find_pre_birth_diffs/) {
	}
	
	return %opt;
}
