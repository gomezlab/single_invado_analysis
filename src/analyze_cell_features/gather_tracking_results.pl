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

use lib "../lib";
use Config::Adhesions;
use Image::Data::Collection;

#Perl built-in variable that controls buffering print output, 1 turns off 
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|config=s", "debug|d", "input|i=s", "output|o=s");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Configuration\n" if $opt{debug};

my @needed_vars =
  qw(data_folder results_folder exp_name single_image_folder raw_data_folder general_data_files lineage_analysis_data_files tracking_output_file);
my $ad_conf = new Config::Adhesions(\%opt, \@needed_vars);
my %cfg = $ad_conf->get_cfg_hash;

###############################################################################
# Main Program
###############################################################################
my %data_sets;

print "\n\nGathering Data Files\n" if $opt{debug};

my @data_files;
push @data_files, split(/\s+/, $cfg{general_data_files}); 
push @data_files, split(/\s+/, $cfg{lineage_analysis_data_files}); 

%data_sets = Image::Data::Collection::gather_data_sets(\%cfg,\%opt,\@data_files);

print "\n\nRemoving Excluded Images\n" if $opt{debug};
%data_sets = Image::Data::Collection::trim_data_sets(\%cfg,\%opt,\%data_sets);

print "\n\nCollecting Tracking Matrix\n" if $opt{debug};
my @tracking_mat = &read_in_tracking_mat;

print "\n\nGathering Adhesion Lineage Properties\n", if $opt{debug};
my %adh_lineage_props;
&gather_adh_lineage_properties;

print "\n\nOutputing Adhesion Lineage Properties\n", if $opt{debug};
&output_adhesion_props;

###############################################################################
# Functions
###############################################################################

#######################################
# Tracking Matrix Collection
#######################################

sub read_in_tracking_mat {
	my @tracking_mat;
	
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

	 return @tracking_mat; 
}

#######################################
# Adhesion Lineage Property Collection
#######################################

sub gather_adh_lineage_properties {
    #@{ $adh_lineage_props{area_shifts} }        = &gather_area_changes;
    @{ $adh_lineage_props{areas} }              = &gather_areas;
    @{ $adh_lineage_props{longevities} }        = &gather_longevities;
    @{ $adh_lineage_props{largest_areas} }        = &gather_largest_areas;
    @{ $adh_lineage_props{starting_edge_dist} } = &gather_starting_dist_from_edge;
}

sub gather_areas {
    my @ad_areas;
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    for my $j (0 .. $#{ $tracking_mat[0] }) {
        my $i_num       = $data_keys[$j];
        my @these_areas = @{ $data_sets{$i_num}{Area} };
        for my $i (0 .. $#tracking_mat) {
            if ($tracking_mat[$i][$j] > -1) {
                $ad_areas[$i][$j] = $these_areas[ $tracking_mat[$i][$j] ];
            } else {
                $ad_areas[$i][$j] = 0;
            }
        }
    }
    return @ad_areas;
}

sub gather_largest_areas {
    my @largest_areas;
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    for my $i (0 .. $#tracking_mat) {
		my $largest = 0;
		for my $j (0 .. $#{ $tracking_mat[$i] }) {
			my $i_num       = $data_keys[$j];
			#my @these_areas = @{ $data_sets{$i_num}{Area} };
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

sub gather_area_changes {
    my @area_shifts;
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    for (0 .. $#data_keys) {
        next if ($_ == $#data_keys);

        my $i_num      = $data_keys[$_];
        my $next_i_num = $data_keys[ $_ + 1 ];
        @area_shifts = &find_ad_area_changes($i_num, $next_i_num, \@area_shifts);
    }
    return @area_shifts;
}

sub find_ad_area_changes {
    my ($i_num, $next_i_num, $areas_ref) = @_;
    my @area_shifts = @$areas_ref;

    my @cur_areas  = @{ $data_sets{$i_num}{Area} };
    my @next_areas = @{ $data_sets{$next_i_num}{Area} };
    my $lin_num    = $#tracking_mat;

    my $area_shift_pos = $#{ $area_shifts[0] } + 1;

    for my $i (0 .. $lin_num) {
        if (   ($tracking_mat[$i][$area_shift_pos] > -1)
            && ($tracking_mat[$i][ $area_shift_pos + 1 ] > -1)) {
            my $starting_ad = $tracking_mat[$i][$area_shift_pos];
            my $ending_ad   = $tracking_mat[$i][ $area_shift_pos + 1 ];

            $area_shifts[$i][$area_shift_pos] = $cur_areas[$starting_ad] - $next_areas[$ending_ad];
        } else {
            $area_shifts[$i][$area_shift_pos] = "nan";
        }
    }
    return @area_shifts;
}

sub gather_starting_dist_from_edge {
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

#######################################
# Output Adhesion Lineage Properties
#######################################

sub output_adhesion_props {
    if (not(-e $cfg{lineage_props_folder})) {
        mkpath($cfg{lineage_props_folder});
    }
    &output_adhesion_longevity(@{ $adh_lineage_props{longevities} });
    #&output_area_changes;

	&output_all;

    &output_adhesion_area_seqeunces(\@{ $adh_lineage_props{areas} }, \@{ $adh_lineage_props{longevities} });
    open OUTPUT, ">edge_dist.csv";
    print OUTPUT join("\n", @{ $adh_lineage_props{starting_edge_dist} });
    close OUTPUT;
}

sub output_all {
    my @longevities = @{ $adh_lineage_props{longevities} };
    my @largest_areas = @{ $adh_lineage_props{largest_areas} };
    my @starting_dists = @{ $adh_lineage_props{starting_edge_dist} };
	
    open OUTPUT, ">$cfg{lineage_props_folder}/all.txt" or die "$!";
	print OUTPUT join(",",qw(longev la sd)),"\n";
	print OUTPUT join("\n", map {"$longevities[$_],$largest_areas[$_],$starting_dists[$_]"} (0 .. $#longevities));
	close OUTPUT;
}

sub output_adhesion_longevity {
    my @longevities = @_;
    open OUTPUT, ">$cfg{lineage_props_folder}/longevities.txt" or die "$!";
    print OUTPUT join("\n", @longevities);
    close OUTPUT;
}

sub output_area_changes {
    my @area_shifts = &remove_nans(@{ $adh_lineage_props{area_shifts} });
}

sub output_adhesion_area_seqeunces {
    my @areas       = @{ $_[0] };
    my @longevities = @{ $_[1] };

    my %adhesion_seq_lengths;
    for my $i (0 .. $#tracking_mat) {
        push @{ $adhesion_seq_lengths{ $longevities[$i] } }, $i;
    }

    if (not(-e "$cfg{lineage_props_folder}/area_shifts")) {
        mkpath("$cfg{lineage_props_folder}/area_shifts");
    }

    my %area_sequences;
    for my $i (keys %adhesion_seq_lengths) {
        for my $j (@{ $adhesion_seq_lengths{$i} }) {
            my @temp;
            for my $k (0 .. $#{ $tracking_mat[$j] }) {
                if ($tracking_mat[$j][$k] > -1) {
                    push @temp, $areas[$j][$k];
                }
            }
            push @{ $area_sequences{$i} }, \@temp;
        }
    }

    for my $i (keys %area_sequences) {
        open OUTPUT, ">$cfg{lineage_props_folder}/area_shifts/$i.txt";

        #print OUTPUT, join(",", @{$adhesion_seq_lengths{$i}}),"\n";
        for (@{ $area_sequences{$i} }) {
            print OUTPUT join(",", @{$_}), "\n";
        }

        close OUTPUT;
    }
}

sub remove_nans {
    my @temp;
    my @area_shifts = @{ $adh_lineage_props{area_shifts} };
    for my $i (0 .. $#area_shifts) {
        for my $j (0 .. $#{ $area_shifts[0] }) {
            if (not($area_shifts[$i][$j] =~ /nan/)) {
                push @{ $temp[$i] }, $area_shifts[$i][$j];
            }
        }
    }
    return @temp;
}
