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

    print "\n\nRemoving Excluded Images\n" if $opt{debug};
    &trim_data_sets;

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

    #check to see if files have been specified for the information used for the
    #tracking, if it has not, die with a message to user to specify these files
    if (not(defined $cfg{tracking_files})) {
        die "ERROR: The files that will be used for tracking must be specified in ",
          "the config file using the name tracking_files\n";
    }

    if (not(defined $cfg{tracking_output_file})) {
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
    my @folders = <$cfg{results_folder}/*/$cfg{raw_data_folder}>;

    my %data_gathered_from;
    foreach my $this_folder (@folders) {
        my $i_num;

        if ($this_folder =~ /$cfg{results_folder}\/(.*)\/$cfg{raw_data_folder}/) {
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

        foreach my $file (split(/\s/, $cfg{tracking_files})) {
            if (-e "$this_folder/$file") {
                @{ $data_sets{$i_num}{$file} } = &gather_data_from_matlab_file("$this_folder/$file");
                push @{ $data_gathered_from{$i_num} }, $file;
                if ($file eq "Centroid") {
                    @{ $data_sets{$i_num}{ $file . "_x" } } = &process_x_centroid_data(@{ $data_sets{$i_num}{$file} });
                    @{ $data_sets{$i_num}{ $file . "_y" } } = &process_y_centroid_data(@{ $data_sets{$i_num}{$file} });
                }
            } else {
                print "ERROR: Problem finding data file ($file) in folder: $this_folder.\n" if $opt{debug};
            }
        }
    }

    if ($opt{debug}) {
        print "Data collected from ", scalar(keys %data_gathered_from), " images. ",
          "Data files gathered for each image include: ", join(", ", split(/\s/, $cfg{tracking_files})), "\n";
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

sub trim_data_sets {
    my @excluded_nums;
    for my $ex_num (@{ $cfg{exclude_image_nums} }) {
        my $match_found = 0;
        for my $this_num (sort { $a <=> $b } keys %data_sets) {
            if ($ex_num == $this_num) {
                delete $data_sets{$ex_num};
                push @excluded_nums, $ex_num;
                $match_found = 1;
            }
        }
        if (not($match_found)) {
            print "ERROR: No match found for the excluded image num ($ex_num)\n";
        }
    }
    if ($opt{debug}) {
        print "Image number removed from further consideration: ", join(", ", @excluded_nums), "\n";
    }
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

        #Gather the Centroid distance matrix
        my @x1 = @{ $data_sets{$key_1}{Centroid_x} };
        my @y1 = @{ $data_sets{$key_1}{Centroid_y} };
        my @x2 = @{ $data_sets{$key_2}{Centroid_x} };
        my @y2 = @{ $data_sets{$key_2}{Centroid_y} };
        @{ $data_sets{$key_1}{Cent_dist} } = &make_dist_mat(\@x1, \@y1, \@x2, \@y2);

        #Gather the Area difference matrix
        my @area1 = @{ $data_sets{$key_1}{Area} };
        my @area2 = @{ $data_sets{$key_2}{Area} };
        @{ $data_sets{$key_1}{Area_diff} } = &make_abs_diff_mat(\@area1, \@area2);
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
    my ($ref_1, $ref_2) = @_;
    my @area1 = @$ref_1;
    my @area2 = @$ref_2;

    my @diff_mat;

    for my $i (0 .. $#area1) {
        for my $j (0 .. $#area2) {
            $diff_mat[$i][$j] = abs($area1[$i] - $area2[$j]);
        }
    }
    return @diff_mat;
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

        print "Working on Image #: $i_num - " if $opt{debug};

        my $num_live_adhesions = &track_live_adhesions($i_num, $next_i_num);

        print "# of Adhesions Being Tracked: $num_live_adhesions - " if $opt{debug};

        my ($num_merged, $num_prob_merged) = &detect_merged_adhesions($i_num, $next_i_num);
        print "# of Merged Ad: $num_merged/$num_prob_merged" if $opt{debug};

        &detect_new_adhesions($next_i_num);
        &fill_in_tracking_mat_holes;
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

    my $live_adhesions = 0;
    for my $i (0 .. $num_ad_lineages) {
        if ($tracking_mat[$i][$cur_step] <= -1) {
            push @{ $tracking_mat[$i] }, $tracking_mat[$i][$cur_step];
            next;
        }

        $live_adhesions++;
        my $adhesion_num           = ${ $tracking_mat[$i] }[$cur_step];
        my @dist_to_next_adhesions = @{ $time_step_dists[$adhesion_num] };

        my @sorted_dist_indexes =
          sort { $dist_to_next_adhesions[$a] <=> $dist_to_next_adhesions[$b] } (0 .. $#dist_to_next_adhesions);

        my @dist_sorted_adhesions_areas = @next_time_step_areas[@sorted_dist_indexes];

        push @{ $tracking_mat[$i] }, $sorted_dist_indexes[0];
    }
    return $live_adhesions;
}

sub detect_merged_adhesions {
    my ($i_num, $next_i_num) = @_;
    my $num_ad_lineages = $#tracking_mat;
    my $cur_step        = $#{ $tracking_mat[0] };
    my @areas           = @{ $data_sets{$i_num}{Area} };
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

            my @merged_areas        = @areas[@merged_ad_nums];
            my @sorted_area_indexes = sort { $merged_areas[$a] <=> $merged_areas[$b] } (0 .. $#merged_areas);
            my $biggest_area_index  = $sorted_area_indexes[$#sorted_area_indexes];

            my @starting_ad = @{ $dest_adhesions{$i}{starting_ad} };
            my @ending_ad   = @{ $dest_adhesions{$i}{ending_ad} };

            my @dist_shifts;
            for (0 .. $#lineage_nums) {
                push @dist_shifts, $dists[ $starting_ad[$_] ][ $ending_ad[$_] ];
            }
            my @sorted_dist_indexes = sort { $dist_shifts[$b] <=> $dist_shifts[$a] } (0 .. $#dist_shifts);
            my $shortest_dist_index = $sorted_dist_indexes[$#sorted_dist_indexes];

            if ($biggest_area_index != $shortest_dist_index) {
                $merged_prob_count++;
                @{ $tracking_probs{merge}{$i_num}{image_nums} }  = ($i_num, $next_i_num);
                @{ $tracking_probs{merge}{$i_num}{starting_ad} } = @starting_ad;
                @{ $tracking_probs{merge}{$i_num}{ending_ad} }   = @ending_ad;
            }

            foreach (0 .. $#lineage_nums) {
                if ($_ != $biggest_area_index) {
                    $tracking_mat[ $lineage_nums[$_] ][$cur_step] *= -1;
                }
            }
        }
    }
    return ($merged_count, $merged_prob_count);
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
}

sub fill_in_tracking_mat_holes {
    my $num_ad_lineages = $#tracking_mat;
    my $max_size        = $#{ $tracking_mat[0] };

    for my $i (0 .. $num_ad_lineages) {
        if ($#{ $tracking_mat[$i] } > $max_size) {
            $max_size = $#{ $tracking_mat[$i] };
        }
    }

    for my $i (0 .. $num_ad_lineages) {
        while ($#{ $tracking_mat[$i] } != $max_size) {
            push @{ $tracking_mat[$i] }, -1;
        }
    }
}

sub min_index {
    my $min_i = 0;
    for (0 .. $#_) {
        if ($_[$_] < $_[$min_i]) {
            $min_i = $_;
        }
    }
    return $min_i;
}

sub max_index {
    my $max_i = 0;
    for (0 .. $#_) {
        if ($_[$_] > $_[$max_i]) {
            $max_i = $_;
        }
    }
    return $max_i;
}

#######################################
# Output Tracking Issues
#######################################
sub output_tracking_probs {
    &output_merge_problems;
}

sub output_merge_problems {
    if (not(-e "$cfg{tracking_probs_folder}/merge")) {
        mkpath("$cfg{tracking_probs_folder}/merge");
    }

    for my $i (keys %{ $tracking_probs{merge} }) {
        mkpath("$cfg{tracking_probs_folder}/merge/$i");
        for my $j (keys %{ $tracking_probs{merge}{$i} }) {
            open OUTPUT, ">$cfg{tracking_probs_folder}/merge/$i/$j.csv";
            print OUTPUT join(",", @{ $tracking_probs{merge}{$i}{$j} });
            close OUTPUT;
        }
    }
}

#######################################
# Output the Tracking Matrix
#######################################
sub output_tracking_mat {
    open OUTPUT, ">$cfg{tracking_output_file}";

    for my $i (0 .. $#tracking_mat) {
        print OUTPUT join(",", @{ $tracking_mat[$i] }), "\n";
    }

    close OUTPUT;
}
