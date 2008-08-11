#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use lib "../lib";
use lib "../lib/perl";

use strict;
use File::Path;
use File::Basename;
use File::Spec::Functions;
use Getopt::Long;
use Data::Dumper;
use Text::CSV::Simple;

use Config::Adhesions qw(ParseConfig);
use Image::Data::Collection;
use Text::CSV::Simple::Extra;
use Emerald;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg=s", "debug|d", "emerald|e");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Configuration\n" if $opt{debug};

my %cfg = ParseConfig(\%opt);

###############################################################################
# Main Program
###############################################################################
if ($opt{emerald}) {
    my $error_folder = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'track_filter');
    mkpath($error_folder);
    
    my %emerald_opt = ("folder" => $error_folder);
    my @command = "$0 -cfg $opt{cfg}";

    @command = &Emerald::create_general_LSF_commands(\@command,\%emerald_opt);
    &Emerald::send_LSF_commands(\@command);
    exit;
}

print "\n\nCollecting Tracking Matrix\n" if $opt{debug};
my @tracking_mat = &read_in_tracking_mat(\%cfg, \%opt);

print "\n\nGathering Adhesion Lineage Properties\n", if $opt{debug};
my %lin_props = &collect_ad_props;

print "\n\nFiltering Tracking Matrix\n" if $opt{debug};
my %filtered_matrix_set = &filter_tracking_matrix;

print "\n\nOutputing Tracking Matrix\n" if $opt{debug};
&output_filtered_matrices(\%filtered_matrix_set,"");

###############################################################################
# Functions
###############################################################################

sub collect_ad_props {
    my $data_file = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{lineage_summary_props_file});

    my $parser = Text::CSV::Simple->new();
    my @data   = $parser->read_file($data_file);

    my @header_line = @{ shift @data };

    my %lin_props;
    foreach (@data) {
        my @data_line = @{$_};
        foreach my $i (0 .. $#header_line) {
            push @{ $lin_props{ $header_line[$i] } }, $data_line[$i];
        }
    }
    return %lin_props;
}

sub filter_tracking_matrix {
    my %matrix_set;

    for my $required_longevity (5,20) {
        for my $i (0 .. $#{ $lin_props{'longevity'} }) {
            my $this_longev = $lin_props{'longevity'}[$i];
            if (   $this_longev >= $required_longevity
                || $this_longev >= $#tracking_mat) {
                push @{ $matrix_set{'longevity'}{$required_longevity} }, $tracking_mat[$i];
            }
        }
    }

    for my $i (0 .. $#tracking_mat) {
        my $in_lin = 0;
        for my $j (0 .. $#{ $tracking_mat[$i] }) {
            $in_lin = 1 if ($tracking_mat[$i][$j] >= 0);
            if ($in_lin && $tracking_mat[$i][$j] <= -1) {
                push @{ $matrix_set{'no_movie'}{'dead'}{'all'} }, $tracking_mat[$i] if ($tracking_mat[$i][$j] == -1);
                $in_lin = 0;
                if ($lin_props{'longevity'}[$i] >= 10) {
                    push @{ $matrix_set{'no_movie'}{'dead'}{'10'} }, $tracking_mat[$i] if ($tracking_mat[$i][$j] == -1);
                }
            }
        }
    }

    #ad-hoc line to pick out specific lineages
    #@{$matrix_set{'special'}{'high_speed'}} = map $tracking_mat[$_], (146,262,516);
    
    #Filter the tracking matrix if a 'for_vis' folder is available
    my $R_sq_folder = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, 'for_vis');
    if (-d $R_sq_folder) {
        for my $file_name (<$R_sq_folder/*.csv>) {
            $file_name =~ /$R_sq_folder\/(.*).csv/;

            my $parser = Text::CSV::Simple->new;
            eval {
                my @data = $parser->read_file($file_name);
                @{$matrix_set{$1}} = map $tracking_mat[$_->[0] - 1], (@data);
            };
        }
    }

    return %matrix_set;
}

sub output_filtered_matrices {
    my %mat_set = %{$_[0]};
    my $prefix = $_[1];
    
    my $base_folder = catdir($cfg{exp_results_folder}, $cfg{tracking_folder}, 'filtered');
    mkpath(catdir($base_folder,$prefix));

    foreach my $i (keys %mat_set) {
        if (ref($mat_set{$i}) eq "ARRAY") {
            output_mat_csv(\@{ $mat_set{$i} }, catfile($base_folder, $prefix , $i . '.csv'));
        } elsif (ref($mat_set{$i}) eq "HASH") {
            output_filtered_matrices(\%{$mat_set{$i}},catdir($prefix,$i));
        } else {
            die "unexpected data type $prefix";
        }
    }
}

################################################################################
#Documentation
################################################################################

=head1 NAME

filter_tracking_matrix.pl - Filter and output the tracking matrix based on collected FA properties 

=head1 SYNOPSIS

filter_tracking_matrix.pl -cfg FA_config

=head1 Description

After collecting the properties of the FA, this program uses those collected properties to filter the tracking matrix to only those adhesions which meet certain criteria. The criteria can be any property collected in earlier stages.

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=back

Optional parameter(s):

=over 

=item * debug or d: print debuging information during program execution

=item * emerald: setups and runs a job tailored for the LSF job system on emerald

=back

=head1 EXAMPLES

filter_tracking_matrix.pl -cfg FA_config

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 6/6/2008 

=cut
