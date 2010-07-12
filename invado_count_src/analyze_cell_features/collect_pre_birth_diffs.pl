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
use File::Find::Rule;
use Getopt::Long;
use Data::Dumper;

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
my @image_folders = <$cfg{individual_results_folder}/*>;

#search for needed files
my @config_files = File::Find::Rule->file()->name( "*.m" )->in( catdir($cfg{exp_results_folder}, $cfg{movie_output_folder}) );
die "Found more than one config file:\n", join("\n\t",@config_files) if (scalar(@config_files) > 1);

my @matlab_code = &create_all_matlab_commands(@image_folders);

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'pre_birth_diffs');
$opt{error_file} = catfile($opt{error_folder}, 'error.txt');
$opt{runtime} = "1";
if (defined $cfg{job_group}) {
    $opt{job_group} = $cfg{job_group};
}

&FA_job::run_matlab_progam(\@matlab_code,\%opt);

################################################################################
#Functions
################################################################################

sub create_all_matlab_commands {
	my @image_folders = @_;
    my @matlab_code;

    foreach my $image_number (1..scalar(@image_folders)) {
        my $extra_opt = "";

		$matlab_code[0] .= "find_pre_birth_diffs('$config_files[0]',$image_number)\n";
    }

    return @matlab_code;
}

################################################################################
#Documentation
################################################################################
