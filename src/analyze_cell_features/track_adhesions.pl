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
use Storable;
use Text::CSV;
use IO::File;

use Config::Adhesions qw(ParseConfig);
use Image::Data::Collection;
use Text::CSV::Simple::Extra;
use Emerald;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|config=s", "debug|d", "input|i=s", "output|o=s", "emerald|e") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Configuration\n" if $opt{debug};

my %cfg = &ParseConfig(\%opt);

###############################################################################
#Main Program
###############################################################################

if ($opt{emerald}) {
    my $error_folder = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'tracking');
    mkpath($error_folder);
    
    my %emerald_opt = ("folder" => $error_folder);
    my @command = "$0 -cfg $opt{cfg} -input data.stor";
    @command = &Emerald::create_general_LSF_commands(\@command,\%emerald_opt);
    &Emerald::send_LSF_commands(\@command);
    exit;
}

my %data_sets;
if (not(defined $opt{input}) || defined $opt{output}) {
    print "\n\nGathering Data Files\n" if $opt{debug};

    my @data_files;
    push @data_files, @{$cfg{general_data_files}};
    push @data_files, @{$cfg{tracking_files}};

    %data_sets = &Image::Data::Collection::gather_data_sets(\%cfg, \%opt, \@data_files);

    print "\n\nMaking Comparison Matrices\n" if $opt{debug};
    &make_comp_matices;
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
#&output_tracking_probs;

print "\n\nOutputing Tracking Matrix\n" if $opt{debug};
my $tracking_output_file = catfile($cfg{exp_results_folder}, $cfg{tracking_folder}, $cfg{tracking_output_file});
mkpath(dirname($tracking_output_file));
&output_mat_csv(\@tracking_mat, $tracking_output_file);

###############################################################################
#Functions
###############################################################################

#######################################
# Process Data Sets
#######################################
sub make_comp_matices {
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    for (0 .. $#data_keys) {

        #The last image can not be compared to a future image, so we skip
        #calculations on it, but still save the image data if the output
        #option is specified
        if ($_ == $#data_keys) {
            if (defined $opt{output}) {
                my $key_1 = $data_keys[$_];
                store \%{ $data_sets{$key_1} }, catfile($cfg{individual_results_folder}, $key_1, $opt{output});
                delete $data_sets{$key_1};
            }
            next;
        }

        #These are the keys we will use for all the subsequent matrix creation
        my ($key_1, $key_2) = @data_keys[ $_, $_ + 1 ];

        print "Working on image Number: $key_1 - " if $opt{debug};

        #Gather the Centroid distance matrix
        my @x1 = @{ $data_sets{$key_1}{Centroid_x} };
        my @y1 = @{ $data_sets{$key_1}{Centroid_y} };
        my @x2 = @{ $data_sets{$key_2}{Centroid_x} };
        my @y2 = @{ $data_sets{$key_2}{Centroid_y} };
        @{ $data_sets{$key_1}{Cent_dist} } = &make_dist_mat(\@x1, \@y1, \@x2, \@y2);
        print "Cent_dist Collected - " if $opt{debug};

        #Gather the Area difference matrix
        my @area1 = @{ $data_sets{$key_1}{Area} };
        my @area2 = @{ $data_sets{$key_2}{Area} };
        @{ $data_sets{$key_1}{Area_diff} } = &make_abs_diff_mat(\@area1, \@area2);
        print "Area_diff Collected - " if $opt{debug};

        #Gather the Pixel Similarity matrix
        my @pix_id1 = @{ $data_sets{$key_1}{PixelIdxList} };
        my @pix_id2 = @{ $data_sets{$key_2}{PixelIdxList} };
        @{ $data_sets{$key_1}{Pix_sim} } = &calc_pix_sim(\@pix_id1, \@pix_id2);

        print "Pix_sim Collected" if $opt{debug};
        print "\r"                if $opt{debug};

        if (defined $opt{output}) {
            store \%{ $data_sets{$key_1} }, catfile($cfg{individual_results_folder}, $key_1, $opt{output});
            delete $data_sets{$key_1};
        }
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
                if ($a > 1) {
                    die;
                }
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
    my @data_keys;
    if (scalar(keys %data_sets) == 0) {
        @data_keys = Image::Data::Collection::gather_sorted_image_numbers(\%cfg);
    } else {
        @data_keys = sort { $a <=> $b } keys %data_sets;
    }

    #This loop contains all the commands used to create the tracking matrix.
    #During each loop another time point is added to all lineages. There are
    #four discrete steps in this process, each one addressing a specific
    #operation in adhesion model, either birth, merging/death or survival.
    #
    #	Tracking Matrix Building Operations:
    #
    #	1. Initialize the tracking matrix with all the adhesion numbers present
    #	   in the first image number (initial births)
    #	2. Try to find a match for all the live adhesions in the next frame of
    #	   the image set (survival)
    #	3. Resolve the adhesion lineages that appear to merge, either by
    #	   picking a lineage as the best match or declaring the other lineages
    #	   dead (death/merging)
    #	4. Find any adhesions in the next frame not assigned to a lineage and
    #	   form a new lineage for each adhesion (birth)
    #	5. Repeat 2-5 until the data from each image has been used to build the
    #	   tracking matrix

    for (0 .. $#data_keys - 1) {
        my $i_num      = $data_keys[$_];
        my $next_i_num = $data_keys[ $_ + 1 ];
        
        #Collect data seta and comparison matrices if not present
        if (not(defined $data_sets{$i_num})) {
            if (defined $opt{input}) {
                %{ $data_sets{$i_num} } = %{ retrieve(catfile($cfg{individual_results_folder}, $i_num, $opt{input})) };
            } else {
                die "Unable to find the data sets and comparison matrices for image number \"$i_num\".";
            }
        }
        if (not(defined $data_sets{$next_i_num})) {
            if (defined $opt{input}) {
                %{ $data_sets{$next_i_num} } =
                  %{ retrieve catfile($cfg{individual_results_folder}, $next_i_num, $opt{input}) };
            } else {
                die "Unable to find the data sets and comparison matrices for image number \"$next_i_num\".";
            }
        }
        
        #STEP 1
        #Start the tracking matrix, if this is the first time throught the loop
        &initialize_tracking_mat($data_keys[0]) if $_ == 0;
        &check_tracking_mat_integrity if $_ == 0;

        #Begin tracking
        print "Image #: $i_num - " if $opt{debug};

        #STEP 2
        &track_live_adhesions($i_num);
        print "# Tracked: $tracking_facts{$i_num}{live_adhesions} - " if $opt{debug};

        #STEP 3
        &detect_merged_adhesions($i_num, $next_i_num);
        print "# Merged: $tracking_facts{$i_num}{merged_count}/$tracking_facts{$i_num}{merged_prob_count} - "
          if $opt{debug};

        #STEP 4
        &detect_new_adhesions($i_num, $next_i_num);
        print "# New: $tracking_facts{$i_num}{new_count}" if $opt{debug};

        &check_tracking_mat_integrity;
        print "\r" if $opt{debug};

        delete $data_sets{$i_num};
    }
}

sub initialize_tracking_mat {
    my $first_key = (sort { $a <=> $b } keys %data_sets)[0];

    my $num_adhesions = $#{ $data_sets{$first_key}{Area} };

    for (0 .. $num_adhesions) {
        push @{ $tracking_mat[$_] }, $_;
    }
    print "Intialized the tracking matrix on image #: $first_key\n" if $opt{debug};
}

sub track_live_adhesions {
    my ($i_num) = @_;
    my $cur_step = $#{ $tracking_mat[0] };

    my @time_step_dists  = @{ $data_sets{$i_num}{Cent_dist} };
    my @time_step_p_sims = @{ $data_sets{$i_num}{Pix_sim} };

    my $pix_sim_indeter_percent = 0.8;
    $pix_sim_indeter_percent = $cfg{pix_sim_indeter_percent} if defined $cfg{pix_sim_indeter_percent};

    #This function makes a tracking guess for all the adhesion lineages living
    #in the current image. The tracking guess is based on two factors: the
    #distance to the adhesions in the next image and the percentage of pixels
    #that overlap with the each of the adhesions in the next image (pixel
    #similarity). There are several cases we need to deal with in making the
    #tracking guess:
    #
    #    Live Tracking Cases:
    #
    #    1. There is a single pixel similarity value that is much higher than the
    #       other values, select the adhesion with the best similarity
    #
    #    2. There are multiple adhesions with close pixel similarity value to the
    #       current adhesion, select the adhesion with the smallest distance between
    #       the centriods
    #
    #    3. There is no pixel similarity to the next adhesion set, select the
    #       adehsion with the smallest distance between the centriods
    #
    #    Notes:
    #
    #    -how close the similarity measures have to be to trigger using centroid
    #    distance between adhesions is defined in $pix_sim_indeter_percent, which
    #    defaults to 0.8

    for my $i (0 .. $#tracking_mat) {
        
        #The tracking matrix code for dead lineage is any number less than or
        #equal to -1, add another -1 to those lineages to make sure they stay
        #the proper length, then skip to the next lineage
        if ($tracking_mat[$i][$cur_step] <= -1) {
            push @{ $tracking_mat[$i] }, -1;
            next;
        }

        $tracking_facts{$i_num}{live_adhesions}++;
        my $adhesion_num      = ${ $tracking_mat[$i] }[$cur_step];
        my @dist_to_next_ads  = @{ $time_step_dists[$adhesion_num] };
        my @p_sim_to_next_ads = @{ $time_step_p_sims[$adhesion_num] };

        my @sorted_dist_indexes =
          sort { $dist_to_next_ads[$a] <=> $dist_to_next_ads[$b] } (0 .. $#dist_to_next_ads);

        my @sorted_p_sim_indexes =
          sort { $p_sim_to_next_ads[$b] <=> $p_sim_to_next_ads[$a] } (0 .. $#p_sim_to_next_ads);

        my $high_p_sim = $p_sim_to_next_ads[ $sorted_p_sim_indexes[0] ];

        if ($sorted_p_sim_indexes[0] != $sorted_dist_indexes[0] && $high_p_sim > 0) {
            $tracking_facts{$i_num}{dist_p_sim_guess_diff}++;
        }

        my $tracking_guess;
        #Case 1 and 2
        if ($high_p_sim > 0) {
            my @p_sim_close_ad_nums = grep {
                if ($p_sim_to_next_ads[$_] >= $high_p_sim * $pix_sim_indeter_percent) {
                    1;
                } else {
                    0;
                }
            } (0 .. $#p_sim_to_next_ads);
            
            my @close_p_sim_by_dist_indexes = sort { $dist_to_next_ads[$a] <=> $dist_to_next_ads[$b] } @p_sim_close_ad_nums;

            if (scalar(@p_sim_close_ad_nums) > 1) {
                foreach my $ad_num (@close_p_sim_by_dist_indexes[1.. $#close_p_sim_by_dist_indexes]) {
                    if (exists $tracking_facts{$i_num}{split_birth_quant}[$ad_num]) {
                        if ($tracking_facts{$i_num}{split_birth_quant}[$ad_num] < $p_sim_to_next_ads[$ad_num]) {
                            $tracking_facts{$i_num}{split_birth}[$ad_num] = $close_p_sim_by_dist_indexes[0];
                            $tracking_facts{$i_num}{split_birth_quant}[$ad_num] = $p_sim_to_next_ads[$ad_num];
                        }
                    } else {
                        $tracking_facts{$i_num}{split_birth}[$ad_num] = $close_p_sim_by_dist_indexes[0];
                        $tracking_facts{$i_num}{split_birth_quant}[$ad_num] = $p_sim_to_next_ads[$ad_num];
                    }
                }
                $tracking_facts{$i_num}{multiple_good_p_sims}++;
                if ($close_p_sim_by_dist_indexes[0] != $sorted_p_sim_indexes[0]) {
                    $tracking_facts{$i_num}{best_pix_sim_not_selected}++;
                }
            }

            $tracking_guess = $close_p_sim_by_dist_indexes[0];
        } else {

            #Case 3
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
    my @dists           = @{ $data_sets{$i_num}{Cent_dist} };
    my @pix_sims        = @{ $data_sets{$i_num}{Pix_sim} };

    # %dest_adhesions will act as the lookup table for adhesions predicted to
    # merge
    my %dest_adhesions;
    for my $i (0 .. $num_ad_lineages) {
        next if $tracking_mat[$i][$cur_step] <= -1;

        my $this_ad = $tracking_mat[$i][$cur_step];
        push @{ $dest_adhesions{$this_ad}{lineage_nums} }, $i;
        push @{ $dest_adhesions{$this_ad}{starting_ad} },  $tracking_mat[$i][ $cur_step - 1 ];
        push @{ $dest_adhesions{$this_ad}{ending_ad} },    $tracking_mat[$i][$cur_step];
    }

    $tracking_facts{$i_num}{merged_count}      = 0;
    $tracking_facts{$i_num}{merged_prob_count} = 0;
    for my $i (keys %dest_adhesions) {
        next if scalar(@{ $dest_adhesions{$i}{lineage_nums} }) == 1;

        $tracking_facts{$i_num}{merged_count}++;

        my @lineage_nums = @{ $dest_adhesions{$i}{lineage_nums} };
        my @starting_ads = @{ $dest_adhesions{$i}{starting_ad} };
        my @ending_ad    = @{ $dest_adhesions{$i}{ending_ad} };

        my @merged_areas = @areas[@starting_ads];
        my @dist_shifts  = map { $dists[ $starting_ads[$_] ][ $ending_ad[$_] ] } (0 .. $#lineage_nums);
        my @pix_sims_set = map { $pix_sims[ $starting_ads[$_] ][ $ending_ad[$_] ] } (0 .. $#lineage_nums);

        my @merge_decisions = &select_best_merge_decision(\@merged_areas, \@dist_shifts, \@pix_sims_set);

        foreach (0 .. $#lineage_nums) {
            if ($merge_decisions[$_]{dead}) {
                $tracking_mat[ $lineage_nums[$_] ][$cur_step] = -1;
            } elsif (not $merge_decisions[$_]{winner}) {
                $tracking_mat[ $lineage_nums[$_] ][$cur_step] =
                  -1 * ($tracking_mat[ $lineage_nums[$_] ][$cur_step] + 2);
            }
        }
    }
}

sub select_best_merge_decision {
    my @merged_areas = @{ $_[0] };
    my @dist_shifts  = @{ $_[1] };
    my @pix_sims_set = @{ $_[2] };

    my @sorted_area_indexes = sort { $merged_areas[$b] <=> $merged_areas[$a] } (0 .. $#merged_areas);
    my $biggest_area_index = $sorted_area_indexes[0];

    my @sorted_dist_indexes = sort { $dist_shifts[$a] <=> $dist_shifts[$b] } (0 .. $#dist_shifts);
    my $shortest_dist_index = $sorted_dist_indexes[0];

    #There are several cases to deal with in picking the adhesion which will
    #continue:
    #
    # 1. The adhesion with the biggest starting area has to shift the least
    # amount, choose the largest adhesion.
    #
    # 2. The differences in area between the largest adhesions in the merge is
    # negligible, but there is a difference in the centroid shifts.  In this
    # case, choose the adhesion from the large area adhesions which shifts the
    # least distance.

    my @merge_decisions = map { { winner => 0, dead => 0, case => 0 } } (0 .. $#merged_areas);

    #Case 1
    if ($biggest_area_index == $shortest_dist_index) {
        $merge_decisions[$biggest_area_index]{winner} = 1;
        $merge_decisions[$biggest_area_index]{dead}   = 0;
        $merge_decisions[$biggest_area_index]{case}   = 1;
    } else {

        #Case 2
        my $shift_percent = 0.75;
        $shift_percent = $cfg{merge_shift_percent} if defined($cfg{merge_shift_percent});

        my @areas_close_indexes = grep {
            if ($merged_areas[$biggest_area_index] / $merged_areas[$_]) {
                1;
            } else {
                0;
            }
        } (0 .. $#merged_areas);

        #GET COUNTS ON NUMBER OF SINGLE MEMBER

        my @areas_close_dist_sort = sort { $dist_shifts[$a] <=> $dist_shifts[$b] } (@areas_close_indexes);

        $merge_decisions[ $areas_close_dist_sort[0] ]{winner} = 1;
        $merge_decisions[ $areas_close_dist_sort[0] ]{dead}   = 0;
        $merge_decisions[ $areas_close_dist_sort[0] ]{case}   = 2;
    }
    
    #To detect dead adhesions we will examine the pixel similarities counts, if
    #the similarity to the next adhesion is zero, then we know there wasn't any
    #overlap with the merged adhesion and it is unlikely that the event
    #constitutes a merge. Instead, it is more likely that we are observing a
    #death event.
    for (0 .. $#merge_decisions) {
        if ($pix_sims_set[$_] == 0) {
            $merge_decisions[$_]{dead} = 1;
        }
    }

    return @merge_decisions;
}

sub detect_new_adhesions {
    my ($i_num, $next_i_num) = @_;

    my $expected_ad_count = $#{ $data_sets{$next_i_num}{Area} };
    my %expected_ad_nums  = map { $_ => 0 } (0 .. $expected_ad_count);
    my $lineage_length    = $#{ $tracking_mat[0] };

    for my $i (0 .. $#tracking_mat) {
        next if ($tracking_mat[$i][$lineage_length] <= -1);
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
            #Check to see if this unassigned adhesion was part of a split birth,
            #if so replace the last entry in the tracking matrix with the
            #adhesion number that it merged from + 1 times -1
            if (defined $tracking_facts{$i_num}{split_birth}[$i] && $tracking_facts{$i_num}{split_birth}[$i] >= 0) {
                $temp[$#temp] = -1*($tracking_facts{$i_num}{split_birth}[$i] + 2);
            }
            push @temp,         $i;
            push @tracking_mat, \@temp;
        }
    }
}

sub check_tracking_mat_integrity {
    if (grep $#{$tracking_mat[$_]} != $#{$tracking_mat[0]}, (0 .. $#tracking_mat)) {
        die "Tracking matrix entries not at same length on cycle ", $#{$tracking_mat[0]};
    }

    my @assigned_ad_nums = map {
        if ($tracking_mat[$_][-1] >= 0) {
            $tracking_mat[$_][-1];
        } else {
            ();
        }
    } (0 .. $#tracking_mat);
    @assigned_ad_nums = sort {$a <=> $b} @assigned_ad_nums;
    
    for (0 .. $#assigned_ad_nums) {
        if ($_ != $assigned_ad_nums[$_]) {
            print Dumper(@assigned_ad_nums);
            print "$_ => $assigned_ad_nums[$_]\n";
            die "Not all ad numbers filled in on cycle ",  $#{$tracking_mat[0]} + 1;
        }
    }
}

#######################################
# Output Tracking Facts
#######################################
sub output_tracking_facts {
    print "# of Adhesion Lineages: ", scalar(@tracking_mat), "\n";
    print "# of Live Adhesions Tracked: ",      &get_all_i_num_count("live_adhesions"),            "\n";
    print "# of Best Pixel Sim Not Selected: ", &get_all_i_num_count("best_pix_sim_not_selected"), "\n";
    print "# of Multiple Good P Sims: ",        &get_all_i_num_count("multiple_good_p_sims"),      "\n";
    print "# of Dist and P Sim Guesses Diff: ", &get_all_i_num_count("dist_p_sim_guess_diff"),     "\n";
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

    mkpath(catdir($full_probs_folder, "merge"));

    for my $i (keys %{ $tracking_probs{merge} }) {
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

################################################################################
#Documentation
################################################################################

=head1 NAME

track_adhesions.pl - Track the focal adhesions identified in prior steps 

=head1 SYNOPSIS

track_adhesions.pl -cfg FA_config

=head1 Description

This program handles tracking the FA from birth until death. The process is
split into four parts:

=over 

=item 1. Initialize the tracking sequence on the FA present in the first image

=item 2. Find the best match for each FA in the next image

=item 3. Resolve any conflicts where two FA are predicted to become the same FA,
either by merging or death

=item 4. Add any unassigned FA into the list of current FA and execute steps 2-4
until all the image are cycled through

=back

The results of this algorithm is a matrix with assigns each adhesion number in
each image to part of a sequence, which identifies the identified adhesions in
each image that make up the entire adhesion life cycle. Functionally, the matrix
is output as a csv file that is used by other programs to determine the
properties of the adhesions and visualize them.

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=back

Optional parameter(s):

=over 

=item * input or i: Look for files in each inidividual image directory with this
file name, which should be a Storable object created by an earlier program
execution 

=item * output or o: Create files using the spectified file name in each
individual directory with Storable containing FA information; can be loaded
using the -i flag

=item * debug or d: print debuging information during program execution

=item * emerald: setups and runs a job tailored for the LSF job system on emerald

=back

=head1 EXAMPLES

track_adhesions.pl -cfg FA_config

OR

track_adhesions.pl -cfg FA_config -i data.stor -o data.stor

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 6/6/2008 

=cut
