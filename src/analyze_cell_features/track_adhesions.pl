#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Path;
use File::Spec::Functions;
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
  qw(data_folder results_folder exp_name single_image_folder raw_data_folder general_data_files tracking_files tracking_output_file);
my $ad_conf = new Config::Adhesions(\%opt, \@needed_vars);
my %cfg = $ad_conf->get_cfg_hash;

###############################################################################
# Main Program
###############################################################################
my %data_sets;
if (not(defined $opt{input}) || not(-e $opt{input})) {
    print "\n\nGathering Data Files\n" if $opt{debug};
	
	my @data_files;
	push @data_files, split(/\s+/, $cfg{general_data_files}); 
	push @data_files, split(/\s+/, $cfg{tracking_files}); 

	%data_sets = Image::Data::Collection::gather_data_sets(\%cfg,\%opt,\@data_files);

    print "\n\nMaking Comparison Matrices\n" if $opt{debug};
    &make_comp_matices;
} else {
    print "\n\nGathering Data Files/Comparison Matrices from Data File\n" if $opt{debug};
    %data_sets = %{ retrieve($opt{input}) };
    %data_sets = Image::Data::Collection::trim_data_sets(\%cfg,\%opt,\%data_sets);
}

if (defined $opt{output}) {
    print "\n\nStoring Data Files/Comparison Matrices\n" if $opt{debug};
    store \%data_sets, $opt{output};
}

print "\n\nDetermining Tracking Matrix\n" if $opt{debug};
my @tracking_mat;
my %tracking_probs;
my %tracking_facts;
&make_tracking_mat;

print "\n\nTracking Results\n" if $opt{debug};
&output_tracking_facts         if $opt{debug};
print "\n"                     if $opt{debug};

print "\n\nOutputing Tracking Problem Data\n" if $opt{debug};
&output_tracking_probs;

print "\n\nOutputing Tracking Matrix\n" if $opt{debug};
&output_tracking_mat;

###############################################################################
# Functions
###############################################################################

#######################################
# Process Data Sets
#######################################
sub make_comp_matices {
    my @data_keys = sort { $a <=> $b } keys %data_sets;

    for (0 .. $#data_keys) {

        #The last image can not be compared to a future image, so we skip
        #calculations on it
        next if ($_ == $#data_keys);

        #These are the keys we will use for all the subsequent matrix creation
        my ($key_1, $key_2) = @data_keys[ $_, $_ + 1 ];

        print "Working on image Number: $key_1 - ";

        #Gather the Centroid distance matrix
        my @x1 = @{ $data_sets{$key_1}{Centroid_x} };
        my @y1 = @{ $data_sets{$key_1}{Centroid_y} };
        my @x2 = @{ $data_sets{$key_2}{Centroid_x} };
        my @y2 = @{ $data_sets{$key_2}{Centroid_y} };
        @{ $data_sets{$key_1}{Cent_dist} } = &make_dist_mat(\@x1, \@y1, \@x2, \@y2);
        print "Cent_dist Collected - ";

        #Gather the Area difference matrix
        my @area1 = @{ $data_sets{$key_1}{Area} };
        my @area2 = @{ $data_sets{$key_2}{Area} };
        @{ $data_sets{$key_1}{Area_diff} } = &make_abs_diff_mat(\@area1, \@area2);
        print "Area_diff Collected - ";

        my @pix_id1 = @{ $data_sets{$key_1}{PixelIdxList} };
        my @pix_id2 = @{ $data_sets{$key_2}{PixelIdxList} };
        @{ $data_sets{$key_1}{Pix_sim} } = &calc_pix_sim(\@pix_id1, \@pix_id2);
        print "Pix_sim Collected - ";
        print "\r";

    }
}

sub make_dist_mat {
    my ($ref_1, $ref_2, $ref_3, $ref_4) = @_;
    my @x1 = @$ref_1;
    my @y1 = @$ref_2;
    my @x2 = @$ref_3;
    my @y2 = @$ref_4;

    my @diff_mat;

    for my $i (0 .. $#x1) {
        for my $j (0 .. $#x2) {
            $diff_mat[$i][$j] = sqrt(($x1[$i] - $x2[$j])**2 + ($y1[$i] - $y2[$j])**2);
        }
    }
    return @diff_mat;
}

sub make_abs_diff_mat {
    my @mat1 = @{ $_[0] };
    my @mat2 = @{ $_[1] };

    my @diff_mat;

    for my $i (0 .. $#mat1) {
        for my $j (0 .. $#mat2) {
            $diff_mat[$i][$j] = abs($mat1[$i] - $mat2[$j]);
        }
    }
    return @diff_mat;
}

sub calc_pix_sim {
    my @pix_id1 = @{ $_[0] };
    my @pix_id2 = @{ $_[1] };

    my @sim_percents;
    for my $i (0 .. $#pix_id1) {
        my @temp_sim;
        my @pix_list = @{ $pix_id1[$i] };
        for my $j (0 .. $#pix_id2) {
            my @matching_list = @{ $pix_id2[$j] };
            my $match_count   = grep {
                my $poss_match = $_;
                my $a = grep $poss_match == $_, @pix_list;
                $a;
            } @matching_list;
            push @temp_sim, $match_count / scalar(@pix_list);
        }
        push @sim_percents, \@temp_sim;
    }
    return @sim_percents;
}

#######################################
# Tracking Matrix Production
#######################################
sub make_tracking_mat {
    &initialize_tracking_mat;
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    for (0 .. $#data_keys) {
        next if ($_ == $#data_keys);

        my $i_num      = $data_keys[$_];
        my $next_i_num = $data_keys[ $_ + 1 ];

        print "Image #: $i_num - " if $opt{debug};

        &track_live_adhesions($i_num, $next_i_num);
        print "# Tracked: $tracking_facts{$i_num}{live_adhesions} - " if $opt{debug};

        &detect_merged_adhesions($i_num, $next_i_num);
        print "# Merged: $tracking_facts{$i_num}{merged_count}/$tracking_facts{$i_num}{merged_prob_count} - "
          if $opt{debug};

        &detect_new_adhesions($i_num, $next_i_num);
        print "# New: $tracking_facts{$i_num}{new_count}" if $opt{debug};

        &check_tracking_mat_integrity;
        print "\r" if $opt{debug};
    }
}

sub initialize_tracking_mat {
    my $first_key = (sort { $a <=> $b } keys %data_sets)[0];
    print "Intializing the tracking matrix on image #: $first_key\n" if $opt{debug};

    my $num_adhesions = $#{ $data_sets{$first_key}{Area} };

    for (0 .. $num_adhesions) {
        push @{ $tracking_mat[$_] }, $_;
    }
}

sub track_live_adhesions {
    my ($i_num, $next_i_num) = @_;
    my $num_ad_lineages = $#tracking_mat;
    my $cur_step        = $#{ $tracking_mat[0] };

    my @time_step_dists = @{ $data_sets{$i_num}{Cent_dist} };
    my @area_diffs      = @{ $data_sets{$i_num}{Area_diff} };
    my @pix_sims        = @{ $data_sets{$i_num}{Pix_sim} };

    my $pix_sim_indeter_percent = 0.8;
    $pix_sim_indeter_percent = $cfg{pix_sim_indeter_percent} if defined $cfg{pix_sim_indeter_percent};

    for my $i (0 .. $num_ad_lineages) {
        if ($tracking_mat[$i][$cur_step] <= -1) {
            push @{ $tracking_mat[$i] }, -1;
            next;
        }

        $tracking_facts{$i_num}{live_adhesions}++;
        my $adhesion_num              = ${ $tracking_mat[$i] }[$cur_step];
        my @dist_to_next_adhesions    = @{ $time_step_dists[$adhesion_num] };
        my @pix_sim_to_next_adhesions = @{ $pix_sims[$adhesion_num] };

        my @sorted_dist_indexes =
          sort { $dist_to_next_adhesions[$a] <=> $dist_to_next_adhesions[$b] } (0 .. $#dist_to_next_adhesions);

        my @sorted_pix_sim_indexes =
          sort { $pix_sim_to_next_adhesions[$b] <=> $pix_sim_to_next_adhesions[$a] } (0 .. $#pix_sim_to_next_adhesions);

        my $tracking_guess;

        my @best_pix_sims = map {
            $pix_sim_to_next_adhesions[ $sorted_pix_sim_indexes[$_] ] -
              $pix_sim_to_next_adhesions[ $sorted_pix_sim_indexes[$#sorted_pix_sim_indexes] ]
        } (0 .. 20);

        if ($best_pix_sims[0] != 0) {
            $tracking_facts{$i_num}{pos_sim_to_next}++;
        }

        if (   $best_pix_sims[0] != 0
            && $best_pix_sims[1] >= $pix_sim_indeter_percent * $best_pix_sims[0]) {
            $tracking_facts{$i_num}{pix_sim_equal}++;
        }

        if ($best_pix_sims[0]) {
            $tracking_guess = $sorted_pix_sim_indexes[0];
        } else {
            $tracking_guess = $sorted_dist_indexes[0];
        }

        push @{ $tracking_mat[$i] }, $tracking_guess;
    }
}

sub detect_merged_adhesions {
    my ($i_num, $next_i_num) = @_;
    my $num_ad_lineages = $#tracking_mat;
    my $cur_step        = $#{ $tracking_mat[0] };
    my @areas           = @{ $data_sets{$i_num}{Area} };
    my @next_step_areas = @{ $data_sets{$next_i_num}{Area} };
    my @dists           = @{ $data_sets{$i_num}{Cent_dist} };
    my @pix_sims        = @{ $data_sets{$i_num}{Pix_sim} };

    my %dest_adhesions;
    for my $i (0 .. $num_ad_lineages) {
        next if $tracking_mat[$i][$cur_step] <= -1;
        push @{ $dest_adhesions{ $tracking_mat[$i][$cur_step] }{lineage_nums} }, $i;
        push @{ $dest_adhesions{ $tracking_mat[$i][$cur_step] }{starting_ad} },  $tracking_mat[$i][ $cur_step - 1 ];
        push @{ $dest_adhesions{ $tracking_mat[$i][$cur_step] }{ending_ad} },    $tracking_mat[$i][$cur_step];
    }

    $tracking_facts{$i_num}{merged_count}      = 0;
    $tracking_facts{$i_num}{merged_prob_count} = 0;
    for my $i (keys %dest_adhesions) {
        next if scalar(@{ $dest_adhesions{$i}{lineage_nums} }) == 1;

        $tracking_facts{$i_num}{merged_count}++;

        my @merged_ad_nums = @{ $dest_adhesions{$i}{starting_ad} };
        my @lineage_nums   = @{ $dest_adhesions{$i}{lineage_nums} };
        my @starting_ads   = @{ $dest_adhesions{$i}{starting_ad} };
        my @ending_ads     = @{ $dest_adhesions{$i}{ending_ad} };

        my @merged_areas = @areas[@merged_ad_nums];

        my @dist_shifts;
        my @pix_sims_set;
        for (0 .. $#lineage_nums) {
            push @dist_shifts,  $dists[ $starting_ads[$_] ][ $ending_ads[$_] ];
            push @pix_sims_set, $pix_sims[ $starting_ads[$_] ][ $ending_ads[$_] ];
        }

        my $final_merged_area = $next_step_areas[ $ending_ads[0] ];

        my ($merge_guess_index, $guess_type) =
          &select_best_merge_decision(\@merged_areas, \@dist_shifts, \@pix_sims_set, $final_merged_area);

        if (!$guess_type) {
            $tracking_facts{$i_num}{merged_prob_count}++;

            push @{ $tracking_probs{merge}{$i_num} },
              {
                image_nums  => [ $i_num, $next_i_num ],
                starting_ad => \@starting_ads,
                ending_ad   => \@ending_ads,
                winning_ad  => [ $merged_ad_nums[$merge_guess_index] ],
              };
        }

        foreach (0 .. $#lineage_nums) {
            if ($_ != $merge_guess_index) {
                $tracking_mat[ $lineage_nums[$_] ][$cur_step] = -1 * $tracking_mat[ $lineage_nums[$_] ][$cur_step] - 1;
            }
        }

    }
}

sub select_best_merge_decision {
    my @merged_areas = @{ $_[0] };
    my @dist_shifts  = @{ $_[1] };
    my @pix_sims_set = @{ $_[2] };
    my $final_area   = $_[3];

    my @sorted_area_indexes = sort { $merged_areas[$b] <=> $merged_areas[$a] } (0 .. $#merged_areas);
    my $biggest_area_index = $sorted_area_indexes[0];

    my @sorted_dist_indexes = sort { $dist_shifts[$a] <=> $dist_shifts[$b] } (0 .. $#dist_shifts);
    my $shortest_dist_index = $sorted_dist_indexes[0];

    my @sorted_pix_sim_indexes = sort { $pix_sims_set[$b] <=> $pix_sims_set[$a] } (0 .. $#pix_sims_set);
    my $best_pix_sim_index = $sorted_pix_sim_indexes[0];

    #There are several cases to deal with here:
    #
    # 1. The adhesion with the biggest starting area has to shift the least
    # amount and has the highest pixel similarity to the final merged adhesion.
    # This is expected when a large adhesion swallows a small one.
    #
    # 2. The differences in area are negligible or equal, but there is a
    # difference in the centroid shifts. In this case, choose the adhesion
    # which is closest to the final adhesion.

    #Case 1
    if (   $biggest_area_index == $shortest_dist_index
        && $biggest_area_index == $best_pix_sim_index
        && $shortest_dist_index == $best_pix_sim_index) {
        return ($biggest_area_index, 1);
    }

    #Case 2
    my $area_shifts_small = 0;
    my $shift_percent     = 0.25;
    $shift_percent = $cfg{merge_shift_percent} if defined($cfg{merge_shift_percent});
    my $second_biggest_area_index = $sorted_area_indexes[ $#sorted_area_indexes - 1 ];

    my $first_area_diff = $merged_areas[$biggest_area_index] - $merged_areas[$second_biggest_area_index];
    if ($first_area_diff / $merged_areas[$biggest_area_index] <= $shift_percent) {
        $area_shifts_small = 1;
    }
    if ($area_shifts_small) {
        return ($shortest_dist_index, 2);
    }

    #Case 3
    my @sorted_pix_sim_diffs = map {
        $pix_sims_set[ $sorted_pix_sim_indexes[$_] ] -
          $pix_sims_set[ $sorted_pix_sim_indexes[$#sorted_pix_sim_indexes] ]
    } (0 .. $#sorted_pix_sim_indexes);

    if ($sorted_pix_sim_diffs[0] != 0 && $sorted_pix_sim_diffs[1] == 0) {
        return ($sorted_pix_sim_indexes[0], 3);
    }

    return ($biggest_area_index, 0);
}

sub detect_new_adhesions {
    my ($i_num, $next_i_num) = @_;

    my $expected_ad_count = $#{ $data_sets{$next_i_num}{Area} };
    my %expected_ad_nums  = map { $_ => 0 } (0 .. $expected_ad_count);
    my $lineage_length    = $#{ $tracking_mat[0] };

    for my $i (0 .. $#tracking_mat) {
        next if ($tracking_mat[$i][$lineage_length] == -1);
        $expected_ad_nums{ $tracking_mat[$i][$lineage_length] } = 1;
    }

    my $new_adhesions_count = 0;
    for (keys %expected_ad_nums) {
        $tracking_facts{$i_num}{new_count}++;
    }

    for my $i (sort { $a <=> $b } keys %expected_ad_nums) {
        if (not($expected_ad_nums{$i})) {
            my @temp;
            for (0 .. $lineage_length - 1) {
                push @temp, -1;
            }
            push @temp,         $i;
            push @tracking_mat, \@temp;
        }
    }
}

sub check_tracking_mat_integrity {
    my $num_ad_lineages = $#tracking_mat;
    my $max_size        = $#{ $tracking_mat[0] };

    for my $i (0 .. $num_ad_lineages) {
        if ($#{ $tracking_mat[$i] } > $max_size) {
            $max_size = $#{ $tracking_mat[$i] };
        }
    }

    my $prob = 0;
    for my $i (0 .. $num_ad_lineages) {
        $prob = 1 if $#{ $tracking_mat[$i] } != $max_size;
    }

    warn "Problem with tracking matrix." if $prob;
}

#######################################
# Output Tracking Facts
#######################################
sub output_tracking_facts {
    print "# of Adhesion Lineages: ", scalar(@tracking_mat), "\n";
    print "# of Live Adhesions Tracked: ",           &get_all_i_num_count("live_adhesions"),  "\n";
    print "# of Positive Pixel Similarities: ",      &get_all_i_num_count("pos_sim_to_next"), "\n";
    print "# of Indeterminate Pixel Similarities: ", &get_all_i_num_count("pix_sim_equal"),   "\n";
    print "# of Merge Operations/# Merge Problems: ",
      &get_all_i_num_count("merged_count"), "/", &get_all_i_num_count("merged_prob_count");
}

sub get_all_i_num_count {
    my $prop  = $_[0];
    my $total = 0;
    for (keys %tracking_facts) {
        next if not defined $tracking_facts{$_}{$prop};
        $total += $tracking_facts{$_}{$prop};
    }
    return $total;
}

#######################################
# Output Tracking Issues
#######################################
sub output_tracking_probs {
    if (not defined $cfg{tracking_probs_folder}) {
        print "Could not find variable tracking_probs_folder in config, skipping output of tracking problems\n"
          if $cfg{debug};
    } else {
        &output_merge_problems;
    }
}

sub output_merge_problems {
    my $full_probs_folder = catdir($cfg{exp_results_folder}, $cfg{tracking_probs_folder});
    if (not(-e catdir($full_probs_folder, "merge"))) {
        mkpath(catdir($full_probs_folder, "merge"));
    }

    for my $i (keys %{ $tracking_probs{merge} }) {
        mkpath(catdir($full_probs_folder, "merge", $i));
        for my $j (0 .. $#{ $tracking_probs{merge}{$i} }) {
            mkpath(catdir($full_probs_folder, "merge", $i, $j));
            my %merge_prob = %{ $tracking_probs{merge}{$i}->[$j] };
            for my $k (keys %merge_prob) {
                open OUTPUT, ">" . catfile($full_probs_folder, "merge", $i, $j, "$k.csv");
                print OUTPUT join(",", @{ $merge_prob{$k} });
                close OUTPUT;
            }
        }
    }
}

#######################################
# Output the Tracking Matrix
#######################################
sub output_tracking_mat {
    open OUTPUT, ">" . catdir($cfg{exp_results_folder}, $cfg{tracking_output_file});

    for my $i (0 .. $#tracking_mat) {
        print OUTPUT join(",", @{ $tracking_mat[$i] }), "\n";
    }

    close OUTPUT;
}
