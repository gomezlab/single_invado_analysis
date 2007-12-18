#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Temp qw/ tempfile tempdir /;
use File::Spec::Functions;
use File::Path;
use Getopt::Long;
use Data::Dumper;
use Storable;
use Statistics::Descriptive;
use Text::CSV;
use IO::File;

#local libraries
use lib "../lib";
use Config::Adhesions;
use Image::Data::Collection;
use Lineage::Properties;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|config=s", "debug|d", "input|i=s", "output|o=s");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Configuration\n" if $opt{debug};

my @needed_vars =
  qw(data_folder results_folder exp_name single_image_folder raw_data_folder general_data_files lineage_analysis_data_files tracking_output_file single_lineage_props_file);
my $ad_conf = new Config::Adhesions(\%opt, \@needed_vars);
my %cfg = $ad_conf->get_cfg_hash;

###############################################################################
# Main Program
###############################################################################

print "\n\nGathering Data Files\n" if $opt{debug};

my @data_files;
push @data_files, split(/\s+/, $cfg{general_data_files});
push @data_files, split(/\s+/, $cfg{lineage_analysis_data_files});
my %data_sets = Image::Data::Collection::gather_data_sets(\%cfg, \%opt, \@data_files);

print "\n\nRemoving Excluded Images\n" if $opt{debug};
%data_sets = Image::Data::Collection::trim_data_sets(\%cfg, \%opt, \%data_sets);

print "\n\nCollecting Tracking Matrix\n" if $opt{debug};
my @tracking_mat = &Image::Data::Collection::read_in_tracking_mat(\%cfg, \%opt);

print "\n\nGathering Adhesion Lineage Properties\n", if $opt{debug};
my %adh_lineage_props = &gather_adh_lineage_properties(\@tracking_mat, \%data_sets);

print "\n\nOutputing Adhesion Lineage Properties\n", if $opt{debug};
&output_adhesion_props(\%cfg, \%adh_lineage_props);

###############################################################################
# Functions
###############################################################################

#######################################
# Adhesion Lineage Property Collection
#######################################

sub gather_adh_lineage_properties {
    my @tracking_mat = @{ $_[0] };
    my %data_sets    = %{ $_[1] };

    my @para_set = (\@tracking_mat, \%data_sets);

    my %adh_lineage_props;
    @{ $adh_lineage_props{longevities} }        = &gather_longevities(@para_set);
    @{ $adh_lineage_props{largest_areas} }      = &gather_largest_areas(@para_set);
    @{ $adh_lineage_props{starting_edge_dist} } = &gather_starting_dist_from_edge(@para_set);
    @{ $adh_lineage_props{speeds} }             = &gather_adhesion_speeds(@para_set);
    @{ $adh_lineage_props{ad_sig} }             = &gather_average_ad_sig(@para_set);

    return %adh_lineage_props;
}

sub gather_largest_areas {
    my @tracking_mat = @{ $_[0] };
    my %data_sets    = %{ $_[1] };

    my @largest_areas;
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    for my $i (0 .. $#tracking_mat) {
        my $largest = 0;
        for my $j (0 .. $#{ $tracking_mat[$i] }) {
            my $i_num = $data_keys[$j];
            if (not defined ${ $data_sets{$i_num}{Area} }[ $tracking_mat[$i][$j] ]) {
                print "$i_num, $tracking_mat[$i][$j]";
                die;
            }
            if (   $tracking_mat[$i][$j] > -1
                && $largest < ${ $data_sets{$i_num}{Area} }[ $tracking_mat[$i][$j] ]) {
                $largest = ${ $data_sets{$i_num}{Area} }[ $tracking_mat[$i][$j] ];
            }
        }
        push @largest_areas, $largest;
    }
    return @largest_areas;
}

sub gather_longevities {
    my @tracking_mat = @{ $_[0] };
    my %data_sets    = %{ $_[1] };

    my @longevities;
    for my $i (0 .. $#tracking_mat) {
        my $count = 0;
        for my $j (0 .. $#{ $tracking_mat[$i] }) {
            $count++ if ($tracking_mat[$i][$j] > -1);
        }
        push @longevities, $count;
    }
    return @longevities;
}

sub gather_starting_dist_from_edge {
    my @tracking_mat = @{ $_[0] };
    my %data_sets    = %{ $_[1] };

    my @starting_dists;
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    for my $i (0 .. $#tracking_mat) {
        for my $j (0 .. $#{ $tracking_mat[$i] }) {
            if ($tracking_mat[$i][$j] > -1) {
                my $adhesion_num = $tracking_mat[$i][$j];
                $starting_dists[$i] = ${ $data_sets{ $data_keys[$j] }{Centroid_dist_from_edge} }[$adhesion_num];
                last;
            }
        }
    }
    return @starting_dists;
}

sub gather_adhesion_speeds {
    my @tracking_mat = @{ $_[0] };
    my %data_sets    = %{ $_[1] };

    my @speeds;
    my @speed_vars;
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    for my $i (0 .. $#tracking_mat) {
        my @temp = ();
        for my $j (0 .. $#{ $tracking_mat[$i] } - 1) {
            if (   $tracking_mat[$i][$j] > -1
                && $tracking_mat[$i][ $j + 1 ] > -1) {
                my $start_ad_num = $tracking_mat[$i][$j];
                my $end_ad_num   = $tracking_mat[$i][ $j + 1 ];

                my $start_x = ${ $data_sets{ $data_keys[$j] }{Centroid_x} }[$start_ad_num];
                my $start_y = ${ $data_sets{ $data_keys[$j] }{Centroid_y} }[$start_ad_num];

                my $end_x = ${ $data_sets{ $data_keys[ $j + 1 ] }{Centroid_x} }[$end_ad_num];
                my $end_y = ${ $data_sets{ $data_keys[ $j + 1 ] }{Centroid_y} }[$end_ad_num];

                push @temp, sqrt(($start_x - $end_x)**2 + ($start_y - $end_y)**2);
            }
        }
        if (@temp) {
            my $stat = Statistics::Descriptive::Full->new();
            $stat->add_data(@temp);

            push @speeds,     $stat->mean();
            push @speed_vars, $stat->variance();
        } else {
            push @speeds,     "NaN";
            push @speed_vars, "NaN";
        }
    }

    return @speeds;
}

sub gather_average_ad_sig {
    my @tracking_mat = @{ $_[0] };
    my %data_sets    = %{ $_[1] };

    my @pax_sig;
    my @pax_sig_vars;
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    for my $i (0 .. $#tracking_mat) {
        my @temp = ();
        for my $j (0 .. $#{ $tracking_mat[$i] }) {
            if ($tracking_mat[$i][$j] > -1) {
                my $ad_num = $tracking_mat[$i][$j];

                push @temp, ${ $data_sets{ $data_keys[$j] }{Average_adhesion_signal} }[$ad_num];
            }
        }
        my $stat = Statistics::Descriptive::Full->new();
        $stat->add_data(@temp);

        push @pax_sig,      $stat->mean();
        push @pax_sig_vars, $stat->variance();
    }

    return @pax_sig;
}

#######################################
# Output Adhesion Lineage Properties
#######################################

sub output_adhesion_props {
    my %cgf               = %{ $_[0] };
    my %adh_lineage_props = %{ $_[1] };

    if (not(-e catdir($cfg{exp_results_folder}, $cfg{lineage_props_folder}))) {
        mkpath(catdir($cfg{exp_results_folder}, $cfg{lineage_props_folder}));
    }

    &output_r_data(\%cfg, \%adh_lineage_props);
}

sub output_r_data {
    my %cgf               = %{ $_[0] };
    my %adh_lineage_props = %{ $_[1] };

    my $output_file = catfile($cfg{exp_results_folder}, $cfg{lineage_props_folder}, $cfg{single_lineage_props_file});

    my @longevities    = @{ $adh_lineage_props{longevities} };
    my @largest_areas  = @{ $adh_lineage_props{largest_areas} };
    my @starting_dists = @{ $adh_lineage_props{starting_edge_dist} };
    my @speeds         = @{ $adh_lineage_props{speeds} };
    my @ad_sig         = @{ $adh_lineage_props{ad_sig} };

    my @all_data =
      map { [ $longevities[$_], $largest_areas[$_], $starting_dists[$_], $speeds[$_], $ad_sig[$_] ] }
      (0 .. $#longevities);

    my $out_hand = new IO::File ">" . $output_file
      or die "Couldn't create single lineage properties file: $output_file";

    my $csv = Text::CSV->new();
    $csv->print($out_hand, [qw(longevity largest_area s_dist_from_edge speed ad_sig)]);
    print $out_hand "\n";
    for my $i (0 .. $#all_data) {
        $csv->print($out_hand, \@{ $all_data[$i] });
        print $out_hand "\n";
    }

    $out_hand->close;
}
