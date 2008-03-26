#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Temp qw/ tempfile tempdir /;
use File::Path;
use File::Spec::Functions;
use Getopt::Long;
use Data::Dumper;
use Storable;

use lib "../lib";
use Config::Adhesions;
use Image::Data::Collection;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg=s", "debug|d", "input|i=s", "output|o=s");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Configuration\n" if $opt{debug};

my @needed_vars =
  qw(tracking_output_file);
my $ad_conf = new Config::Adhesions(\%opt, \@needed_vars);
my %cfg = $ad_conf->get_cfg_hash;
###############################################################################
# Main Program
###############################################################################

print "\n\nCollecting Tracking Matrix\n" if $opt{debug};
my @tracking_mat;
&read_in_tracking_mat;

print "\n\nGathering Adhesion Lineage Properties\n", if $opt{debug};
my %adh_lineage_props;
&gather_adh_lineage_properties;

print "\n\nFiltering Tracking Matrix\n" if $opt{debug};
&filter_adhesion_lineages;

print "\n\nOutputing Tracking Matrix\n" if $opt{debug};
&output_tracking_mat;

###############################################################################
# Functions
###############################################################################

#######################################
# Tracking Matrix Collection
#######################################
sub read_in_tracking_mat {
    open TRACK_IN, catdir($cfg{results_folder},$cfg{exp_name},$cfg{tracking_output_file}) 
	  or die "Tried to open: ",catdir($cfg{results_folder},$cfg{exp_name},$cfg{tracking_output_file});
    foreach (<TRACK_IN>) {
        chomp($_);
        my @temp = split(",", $_);
        push @tracking_mat, \@temp;
    }
    close TRACK_IN;

    print "Gathered ", scalar(@tracking_mat), " lineages, with ", scalar(@{ $tracking_mat[0] }),
      " timepoints from file ", $cfg{tracking_output_file}, ".\n"
      if $opt{debug};
}

#######################################
# Adhesion Lineage Property Collection
#######################################

sub gather_adh_lineage_properties {
    @{ $adh_lineage_props{longevities} } = &gather_longevities;
}

#######################################
# Adhesion Lineage Filtering
#######################################
sub filter_adhesion_lineages {
	my $min_life = 2;

	my @new_track_mat;
	for my $i (0 .. $#tracking_mat) {
		if (${$adh_lineage_props{longevities}}[$i] >= $min_life) {
			push @new_track_mat, $tracking_mat[$i];
		}
	}
	@tracking_mat = @new_track_mat;
}

#######################################
# Output Tracking Matrix
#######################################
sub output_tracking_mat {
    open OUTPUT, ">filtered_tracking_matrix.csv";

    for my $i (0 .. $#tracking_mat) {
        print OUTPUT join(",", @{ $tracking_mat[$i] }), "\n";
    }

    close OUTPUT;
}
