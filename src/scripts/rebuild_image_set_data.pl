#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Path;
use File::Spec::Functions;
use Getopt::Long;
use IO::File;

use lib "../lib";
use Config::Adhesions;
use Image::Data::Collection;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|config=s", "debug|d");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Configuration\n" if $opt{debug};

my @needed_vars =
  qw(data_folder results_folder exp_name single_image_folder raw_data_folder general_data_files tracking_files tracking_output_file);
my $ad_conf = new Config::Adhesions(\%opt, \@needed_vars);
my %cfg = $ad_conf->get_cfg_hash;

###############################################################################
# Main Program
###############################################################################

if (-e $cfg{exp_results_folder}) {
    File::Path::rmtree($cfg{exp_results_folder});
}

chdir "../find_cell_features";
system "./collect_mask_set.pl -d -cfg $opt{cfg}";
system "./collect_fa_image_set.pl -d -cfg $opt{cfg}";

chdir "../analyze_cell_features";
system "./track_adhesions.pl -d -cfg $opt{cfg} -o data.stor -i data.stor";
