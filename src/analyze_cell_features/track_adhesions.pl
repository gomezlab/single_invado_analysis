#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Temp qw/ tempfile tempdir /;
use File::Path;
use Config::General qw/ ParseConfig /;
use Getopt::Long;
use Data::Dumper;
use Storable;
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg=s", "debug|d", "input|i=s", "output|o=s");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Configuration\n" if $opt{debug};
my %cfg = &get_config;

###############################################################################
# Main Program
###############################################################################
my %data_sets;
if (not(defined $opt{input}) || not(-e $opt{input})) {
    print "\n\nGathering Data Files\n" if $opt{debug};
    &gather_data_sets;

    print "\n\nMaking Comparison Matrices\n" if $opt{debug};
    &make_comp_matices;
} else {
    print "\n\nGathering Data Files/Comparison Matrices from Data File\n";
    %data_sets = %{ retrieve($opt{input}) };
}

if (defined $opt{output}) {
    store \%data_sets, $opt{output};
}

print "\n\nDetermining Tracking Matrix\n" if $opt{debug};
my @tracking_mat;
my %tracking_probs;
&make_tracking_mat;

print "\n\nOutputing Tracking Problem Data\n" if $opt{debug};
&output_tracking_probs;

print "\n\nOutputing Tracking Matrix\n" if $opt{debug};
&output_tracking_mat;

###############################################################################
# Functions
###############################################################################

#######################################
# Config Collection
#######################################
sub get_config {
    my %cfg = ParseConfig(
        -ConfigFile            => $opt{cfg},
        -MergeDuplicateOptions => 1,
    );

    #check to see if a file for the frames that should be excluded from the
    #analysis is included, if it is, collect the data from it, otherwise, set
    #exclude_image_nums to 0
    if (defined $cfg{exclude_file}) {
        open EX_INPUT, $cfg{exclude_file} or die "Can't open the specified exclude file: $cfg{exclude_file}";
        my $temp_line = <EX_INPUT>;
        close EX_INPUT;

        chomp($temp_line);
        @{ $cfg{exclude_image_nums} } = split(",", $temp_line);
    } else {
        $cfg{exclude_image_nums} = 0;
    }
    if ($opt{debug}) {
        print "Image numbers to be excluded:", join(", ", @{ $cfg{exclude_image_nums} }), "\n";
    }
	
	if (not defined $cfg{results_folder}) {
		die "ERROR: The location of the folder that contains the results of the focal ",
		  "adhesion identification must be specified in the config file with the ",
		  "variable \"results_folder\".";
	}

	if (not defined $cfg{single_image_folder}) {
		die "ERROR: The location in the results folder where each of the individual ", 
		  "images and the associated data about those images must be specified in ",
		  "the config file with the variable \"single_image_folder\".";
	}

	if (not defined $cfg{raw_data_folder}) {
		die "ERROR: The folder in each of the individual image folders which contains ",
		  "all the collected properties must be specified in the config file with the ",
		  "variable \"raw_data_folder\".";
	}

    if (not defined $cfg{tracking_files}) {
        die "ERROR: The files that will be used for tracking must be specified in ",
          "the config file using the name tracking_files\n";
    }

    if (not defined $cfg{tracking_output_file}) {
        die "ERROR: The file name for writing the tracking data has not been ",
          "specified, the file name must be specified using the variable name ",
          "\"tracking_output_file\" in the config file.\n";
    }

    if (defined $opt{input} && not(-e $opt{input})) {
        warn "Warning: The data file specified after \'-i\' ($opt{input}) ",
          "does not exist, input data will be gathered from standard data files.\n";
    }

    return %cfg;
}

#######################################
# Data Set Collection
#######################################
sub gather_data_sets {

    my @folders = <$cfg{results_folder}/$cfg{single_image_folder}/*/$cfg{raw_data_folder}>;

    my @data_files;
    push @data_files, split(/\s/, $cfg{general_data_files});
    push @data_files, split(/\s/, $cfg{tracking_files});
	my $image_count = 0;

    foreach my $this_folder (@folders) {
        my $i_num;
        if ($this_folder =~ /$cfg{results_folder}\/$cfg{single_image_folder}\/(.*)\/$cfg{raw_data_folder}/) {
            $i_num = $1;
            if (not($i_num =~ /^\d*$/)) {
                print "ERROR: Problem finding image number in folder: $this_folder, ",
                  "skipping folder in further computation.\n"
                  if $opt{debug};
                next;
            }
        } else {
            print "ERROR: Problem finding image number in folder name: $this_folder, ",
              "skipping folder in further computation.\n"
              if $opt{debug};
            next;
        }
		
		#Skip further processin on images in the excluded list
		next if grep $i_num == $_, @{ $cfg{exclude_image_nums} };

		$image_count++;

        foreach my $file (@data_files) {
            if (-e "$this_folder/$file" && -f "$this_folder/$file") {
                @{ $data_sets{$i_num}{$file} } = &gather_data_from_matlab_file("$this_folder/$file");
                if ($file eq "Centroid") {
                    @{ $data_sets{$i_num}{ $file . "_x" } } = &process_x_centroid_data(@{ $data_sets{$i_num}{$file} });
                    @{ $data_sets{$i_num}{ $file . "_y" } } = &process_y_centroid_data(@{ $data_sets{$i_num}{$file} });
                }
                if ($file eq "Area" || $file eq "Centroid_dist_from_edge") {
                    @{ $data_sets{$i_num}{$file} } = map sprintf("%d", $_), @{ $data_sets{$i_num}{$file} };
                }
            } elsif (-e "$this_folder/$file" && -d "$this_folder/$file") {
                if ($file eq "PixelIdxList") {
                    @{ $data_sets{$i_num}{$file} } = &gather_PixelIdxList_data("$this_folder/$file");
                }
            } else {
                print "ERROR: Problem finding data file ($file) in folder: $this_folder.\n" if $opt{debug};
            }
        }
    }
    if ($opt{debug}) {
        print "Data collected from ", $image_count, " images. ",
          "Data files gathered for each image include: ", join(", ", @data_files), "\n";
    }
}

sub process_x_centroid_data {
    my @centroid = @_;
    my @x;

    for (0 .. $#centroid / 2) {
        push @x, $centroid[ $_ * 2 ];
    }

    return @x;
}

sub process_y_centroid_data {
    my @centroid = @_;
    my @y;

    for (0 .. $#centroid / 2) {
        push @y, $centroid[ $_ * 2 + 1 ];
    }

    return @y;
}

sub gather_data_from_matlab_file {
    my ($file) = @_;

    open INPUT, "$file" or die "Problem opening $file";
    my $this_line = <INPUT>;
    chomp($this_line);
    my @in = split("   ", $this_line);
    close INPUT;

    shift @in if ($in[0] eq "");
    return @in;
}

sub gather_PixelIdxList_data {
	my $folder = $_[0];
	
	my @all_pixid;

	my @ad_folders = <$folder/*>;
	my $total_ad   = scalar(@ad_folders);
	
	for my $ad_num (1..$total_ad) {
		my @temp;
		open INPUT, "$folder/$ad_num" or die "Unable to open PixelIdxList file: $folder/$ad_num";
		for (<INPUT>) {
			chomp($_);
			my $converted_line = sprintf("%d",$_);
			if ($converted_line =~ /(\d+)/) {
				push @temp, $1;
			}
		}
		close INPUT;
		push @all_pixid, \@temp;
	}
	return @all_pixid;
}

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

		my @pix_id1 = @{$data_sets{$key_1}{PixelIdxList}};
		my @pix_id2 = @{$data_sets{$key_2}{PixelIdxList}};
		@{ $data_sets{$key_1}{Pix_sim} } = &calc_pix_sim(\@pix_id1,\@pix_id2);
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
    my @mat1 = @{$_[0]};
    my @mat2 = @{$_[1]};

    my @diff_mat;

    for my $i (0 .. $#mat1) {
        for my $j (0 .. $#mat2) {
            $diff_mat[$i][$j] = abs($mat1[$i] - $mat2[$j]);
        }
    }
    return @diff_mat;
}

sub calc_pix_sim {
	my @pix_id1 = @{$_[0]};
	my @pix_id2 = @{$_[1]};

	my @sim_percents;
	for my $i (0 .. $#pix_id1) {
		my @temp_sim;
		my @pix_list = @{$pix_id1[$i]};
		for my $j (0 .. $#pix_id2) {
			my @matching_list = @{$pix_id2[$j]};
			my $match_count = grep { my $poss_match = $_;
							  		 my $a = grep $poss_match == $_, @pix_list;
							  		 $a; } @matching_list;
			push @temp_sim, $match_count/scalar(@pix_list);
		}
		push @sim_percents, \@temp_sim;
	}
	return @sim_percents;
}
#######################################
# Tracking Matrix Collection
#######################################
sub make_tracking_mat {
    &initialize_tracking_mat;
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    for (0 .. $#data_keys) {
        next if ($_ == $#data_keys);

        my $i_num      = $data_keys[$_];
        my $next_i_num = $data_keys[ $_ + 1 ];

        print "Image #: $i_num - " if $opt{debug};

        my $num_live_adhesions = &track_live_adhesions($i_num, $next_i_num);

        print "# Tracked: $num_live_adhesions - " if $opt{debug};

        my ($num_merged, $num_prob_merged) = &detect_merged_adhesions($i_num, $next_i_num);
        print "# Merged: $num_merged/$num_prob_merged - " if $opt{debug};

        my $num_new = &detect_new_adhesions($next_i_num);
		print "# New: $num_new" if $opt{debug};

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

    my @time_step_dists      = @{ $data_sets{$i_num}{Cent_dist} };
    my @this_time_step_areas = @{ $data_sets{$i_num}{Area} };
    my @next_time_step_areas = @{ $data_sets{$next_i_num}{Area} };
	my @area_diffs 			 = @{ $data_sets{$i_num}{Area_diff} };
	my @pix_sims 			 = @{ $data_sets{$i_num}{Pix_sim} };

    my $live_adhesions = 0;
    for my $i (0 .. $num_ad_lineages) {
        if ($tracking_mat[$i][$cur_step] <= -1) {
            push @{ $tracking_mat[$i] }, $tracking_mat[$i][$cur_step];
            next;
        }

        $live_adhesions++;
        my $adhesion_num              = ${ $tracking_mat[$i] }[$cur_step];
        my @dist_to_next_adhesions    = @{ $time_step_dists[$adhesion_num] };
		my @pix_sim_to_next_adhesions = @{ $pix_sims[$adhesion_num]};

        my @sorted_dist_indexes =
          sort { $dist_to_next_adhesions[$a] <=> $dist_to_next_adhesions[$b] } (0 .. $#dist_to_next_adhesions);
		
		my @sorted_pix_index = sort { $pix_sim_to_next_adhesions[$b] <=> $pix_sim_to_next_adhesions[$a] } (0 .. $#dist_to_next_adhesions);
		
		if ($sorted_pix_index[0] != $sorted_dist_indexes[0]) {
			my @props = ($adhesion_num,
						 $this_time_step_areas[$adhesion_num], 
						 $dist_to_next_adhesions[$sorted_dist_indexes[0]], 
						 $dist_to_next_adhesions[$sorted_pix_index[0]],
						 $pix_sim_to_next_adhesions[$sorted_dist_indexes[0]],
						 $pix_sim_to_next_adhesions[$sorted_pix_index[0]],
						 $pix_sim_to_next_adhesions[$sorted_pix_index[0]] - $pix_sim_to_next_adhesions[$sorted_pix_index[$#sorted_pix_index]],
						 $area_diffs[$adhesion_num][$sorted_dist_indexes[0]],
						 $area_diffs[$adhesion_num][$sorted_pix_index[0]],
						 );
			
			open OUTPUT, ">>tracking_probs.csv";
			print OUTPUT join(",",@props),"\n";
			close OUTPUT;
		}
		
		my $tracking_guess;

		my $best_worst_diff = $pix_sim_to_next_adhesions[$sorted_pix_index[0]] - 
							  $pix_sim_to_next_adhesions[$sorted_pix_index[$#sorted_pix_index]];
		
		if ($best_worst_diff > 0.01) {
			$tracking_guess = $sorted_pix_index[0];
		} else {
			$tracking_guess = $sorted_dist_indexes[0];
		}

        push @{ $tracking_mat[$i] }, $tracking_guess;
    }
    return $live_adhesions;
}

sub detect_merged_adhesions {
    my ($i_num, $next_i_num) = @_;
    my $num_ad_lineages = $#tracking_mat;
    my $cur_step        = $#{ $tracking_mat[0] };
    my @areas           = @{ $data_sets{$i_num}{Area} };
    my @next_step_areas = @{ $data_sets{$next_i_num}{Area} };
    my @dists           = @{ $data_sets{$i_num}{Cent_dist} };

    my %dest_adhesions;
    for my $i (0 .. $num_ad_lineages) {
        next if $tracking_mat[$i][$cur_step] <= -1;
        push @{ $dest_adhesions{ $tracking_mat[$i][$cur_step] }{lineage_nums} }, $i;
        push @{ $dest_adhesions{ $tracking_mat[$i][$cur_step] }{starting_ad} },  $tracking_mat[$i][ $cur_step - 1 ];
        push @{ $dest_adhesions{ $tracking_mat[$i][$cur_step] }{ending_ad} },    $tracking_mat[$i][$cur_step];
    }

    my $merged_count      = 0;
    my $merged_prob_count = 0;
    for my $i (keys %dest_adhesions) {
        if (scalar(@{ $dest_adhesions{$i}{lineage_nums} }) > 1) {
            $merged_count++;

            my @merged_ad_nums = @{ $dest_adhesions{$i}{starting_ad} };
            my @lineage_nums   = @{ $dest_adhesions{$i}{lineage_nums} };
            my @starting_ads   = @{ $dest_adhesions{$i}{starting_ad} };
            my @ending_ads     = @{ $dest_adhesions{$i}{ending_ad} };

            my @merged_areas        = @areas[@merged_ad_nums];

            my @dist_shifts;
            for (0 .. $#lineage_nums) {
                push @dist_shifts, $dists[ $starting_ads[$_] ][ $ending_ads[$_] ];
            }

			my $final_merged_area = $next_step_areas[$ending_ads[0]];
			
			my ($merge_guess_index,$guess_type) = &select_best_merge_decision(\@merged_areas, \@dist_shifts, $final_merged_area);

			if (!$guess_type) {
                $merged_prob_count++;
                push @{ $tracking_probs{merge}{$i_num} }, { 
                    image_nums  => [$i_num, $next_i_num],
                    starting_ad => \@starting_ads,
                    ending_ad   => \@ending_ads,
                    winning_ad  => [$merged_ad_nums[$merge_guess_index]],
                };
            }

            foreach (0 .. $#lineage_nums) {
                if ($_ != $merge_guess_index) {
                    $tracking_mat[ $lineage_nums[$_] ][$cur_step] = -1;
                }
            }
        }
    }
    return ($merged_count, $merged_prob_count);
}

sub select_best_merge_decision {
	my @merged_areas = @{$_[0]};
	my @dist_shifts  = @{$_[1]};
	my $final_area   = $_[2];

	my @sorted_area_indexes = sort { $merged_areas[$a] <=> $merged_areas[$b] } (0 .. $#merged_areas);
	my $biggest_area_index  = $sorted_area_indexes[$#sorted_area_indexes];
	
	my @sorted_dist_indexes = sort { $dist_shifts[$a] <=> $dist_shifts[$b] } (0 .. $#dist_shifts);
	my $shortest_dist_index = $sorted_dist_indexes[0];

	#There are several cases to deal with here:
	# 1. The adhesion with the biggest starting area has to shift the least amount.
	#    This is expected when a large adhesion swallows a small one.
	# 2. The differences in area are negligible or equal, but there is a 
	#	 difference in the centroid shifts. In this case, choose the adhesion 
	#    which is closest to the final adhesion. 
	
	#Case 1
	if ($biggest_area_index == $shortest_dist_index) {
		return ($biggest_area_index,1);
	}

	#Case 2
	my $area_shifts_small = 0;
	my $shift_percent = 0.25;
	$shift_percent = $cfg{merge_shift_percent} if defined($cfg{merge_shift_percent});
	my $second_biggest_area_index = $sorted_area_indexes[$#sorted_area_indexes - 1];

	my $first_area_diff = $merged_areas[$biggest_area_index] - $merged_areas[$second_biggest_area_index];
	if ($first_area_diff/$merged_areas[$biggest_area_index]	<= $shift_percent) {
		$area_shifts_small = 1;
	}
	if ($area_shifts_small) {
		return ($shortest_dist_index, 2);
	}

	return ($biggest_area_index, 0);
}

sub detect_new_adhesions {
    my $next_i_num = $_[0];

    my $expected_ad_count = $#{ $data_sets{$next_i_num}{Area} };
    my %expected_ad_nums  = map { $_ => 0 } (0 .. $expected_ad_count);
    my $lineage_length    = $#{ $tracking_mat[0] };

    for my $i (0 .. $#tracking_mat) {
        next if ($tracking_mat[$i][$lineage_length] == -1);
        $expected_ad_nums{ $tracking_mat[$i][$lineage_length] } = 1;
    }

    my $new_adhesions_count = 0;
    for (keys %expected_ad_nums) {
        $new_adhesions_count++ if (not($expected_ad_nums{$_}));
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
	return $new_adhesions_count;
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
# Output Tracking Issues
#######################################
sub output_tracking_probs {
    if (not defined $cfg{tracking_probs_folder}) {
		print "Could not find variable tracking_probs_folder in config, skipping ",
		  "output of tracking problems\n" if $cfg{debug};
		return 0;
	}
	&output_merge_problems;
}

sub output_merge_problems {
	my $full_probs_folder = $cfg{results_folder} . "/" . $cfg{tracking_probs_folder};
	if (not(-e "$full_probs_folder/merge")) {
        mkpath("$full_probs_folder/merge");
    }

    for my $i (keys %{ $tracking_probs{merge} }) {
        mkpath("$full_probs_folder/merge/$i");
        for my $j (0 .. $#{ $tracking_probs{merge}{$i} }) {
        	mkpath("$full_probs_folder/merge/$i/$j");
			my %merge_prob = %{$tracking_probs{merge}{$i}->[$j]};
			for my $k (keys %merge_prob) {
            	open OUTPUT, ">$full_probs_folder/merge/$i/$j/$k.csv";
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
    open OUTPUT, ">$cfg{results_folder}/$cfg{tracking_output_file}";

    for my $i (0 .. $#tracking_mat) {
        print OUTPUT join(",", @{ $tracking_mat[$i] }), "\n";
    }

    close OUTPUT;
}
