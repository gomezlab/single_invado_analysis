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
use File::Copy;
use Math::Matlab::Local;
use Getopt::Long;
use Data::Dumper;

use Config::Adhesions qw(ParseConfig);
use Math::Matlab::Extra;
use Emerald;
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l")
  or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

################################################################################
# Main Program
################################################################################
if ($opt{lsf}) {
    my @commands;
	
	#$0 - the name of the program currently running, used to protect against
	#future file name changes
	push @commands, "$0 -cfg $opt{cfg}";
    
    $opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'find_invadopodia');
    if (defined $cfg{job_group}) {
        $opt{job_group} = $cfg{job_group};
    }

    &FA_job::send_general_lsf_program(\@commands,\%opt);

    exit(0);
}

my $data_dir = catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder});

my $output_base = catfile($cfg{exp_results_folder}, $cfg{errors_folder}, 'find_invadopodia');
if (! -e $output_base) {
	mkpath($output_base);
}

my @R_cmds;

my $output_file = catfile($output_base, 'R_out.txt');
push @R_cmds, "R CMD BATCH --vanilla \'--args data_dir=$data_dir\' invado_analysis_lib.R $output_file";

$opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'find_invadopodia');
$opt{error_file} = catfile($opt{error_folder}, 'error.txt');

for (@R_cmds) {
	if ($opt{debug}) {
		print "$_\n";
	} else {
		system($_);
	}
}

################################################################################
#Functions
################################################################################
