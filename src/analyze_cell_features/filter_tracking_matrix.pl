#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Path;
use File::Spec::Functions;
use Getopt::Long;
use Data::Dumper;
use Text::CSV::Simple;

use lib "../lib";
use Config::Adhesions;
use Image::Data::Collection;
use Text::CSV::Simple::Extra;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg=s", "debug|d");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Configuration\n" if $opt{debug};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg = $ad_conf->get_cfg_hash;
###############################################################################
# Main Program
###############################################################################

print "\n\nCollecting Tracking Matrix\n" if $opt{debug};
my @tracking_mat = &read_in_tracking_mat(\%cfg, \%opt);

print "\n\nGathering Adhesion Lineage Properties\n", if $opt{debug};
my %lin_props = &collect_ad_props;

print "\n\nFiltering Tracking Matrix\n" if $opt{debug};
my %filtered_matrix_set = &filter_tracking_matrix;

print "\n\nOutputing Tracking Matrix\n" if $opt{debug};
&output_trimmed_matrices;

###############################################################################
# Functions
###############################################################################

sub collect_ad_props {
    my $data_file = catfile($cfg{exp_results_folder},$cfg{adhesion_props_folder},$cfg{lineage_summary_props_file});
    
    my $parser = Text::CSV::Simple->new();
    my @data = $parser->read_file($data_file);
    
    my @header_line =  @{shift @data};
    
    my %lin_props;
    foreach (@data) {
        my @data_line = @{$_};
        foreach my $i (0 .. $#header_line) {
            push @{$lin_props{$header_line[$i]}}, $data_line[$i];
        }
    }
    return %lin_props;
}

sub filter_tracking_matrix {
    my %matrix_set;

    for (2 .. 5) {
        my $required_longevity = $_;
        for my $i (0 .. $#{$lin_props{'longevity'}}) {
            my $this_longev = $lin_props{'longevity'}[$i];
            if ($this_longev >= $required_longevity ||
                $this_longev >= $#tracking_mat) {
                push @{$matrix_set{'longevity'}{$required_longevity}}, $tracking_mat[$i];
            }
        }
    }

    for my $i (0 .. $#tracking_mat) {
        my $in_lin = 0;
        for my $j (0 .. $#{$tracking_mat[$i]}) {
            $in_lin = 1 if ($tracking_mat[$i][$j] >= 0);
            if ($in_lin && $tracking_mat[$i][$j] <= -1) {
                push @{$matrix_set{'dead'}{'dead'}}, $tracking_mat[$i] if ($tracking_mat[$i][$j] == -1);
                $in_lin = 0;
                if ($lin_props{'longevity'}[$i] >= 10) {
                    push @{$matrix_set{'dead'}{'10'}}, $tracking_mat[$i] if ($tracking_mat[$i][$j] == -1);
                }
            }
        }
    }
    
    return %matrix_set;
}

sub output_trimmed_matrices {
    my $base_folder = catdir($cfg{exp_results_folder},$cfg{tracking_folder},'filtered');
    mkpath($base_folder);

    foreach my $type (keys %filtered_matrix_set) {
        my $filtered_path = catdir($base_folder, $type);
        mkpath($filtered_path);
        foreach my $i (keys %{$filtered_matrix_set{$type}}) {
            output_mat_csv(\@{$filtered_matrix_set{$type}{$i}}, catfile($filtered_path, $i . '.csv'));
        }
    }
}
