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

print "\n\nGathering Data Files\n" if $opt{debug};
&gather_data_sets;

print "\n\nRemoving Excluded Images\n" if $opt{debug};
&trim_data_sets;

print "\n\nCollecting Tracking Matrix\n" if $opt{debug};
my @tracking_mat;
&gather_tracking_mat;

print "\n\nGathering Adhesion Lineage Properties\n", if $opt{debug};
my %adh_lineage_props;
&gather_adh_lineage_properties;

print "\n\nOutputing Adhesion Lineage Properties\n", if $opt{debug};
&output_adhesion_props;

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
        open EX_INPUT, $cfg{exclude_file}
          or die "Can't open the specified exclude file: $cfg{exclude_file}";
        my $temp_line = <EX_INPUT>;
        close EX_INPUT;

        chomp($temp_line);
        @{ $cfg{exclude_image_nums} } = split(/\s/, $temp_line);

        shift @{ $cfg{exclude_image_nums} } if (${ $cfg{exclude_image_nums} }[0] =~ //);
    } else {
        @{ $cfg{exclude_image_nums} } = (0);
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

    if (not(defined $cfg{lineage_props_folder})) {
        die "ERROR: A folder must be specified in the config file which specifies ",
          "where the lineage properties will be output, the variable name is ", "
		  \"lineage_props_folder\".";
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

                if ($file eq "Area" || $file eq "Centroid_dist_from_edge") {
                    @{ $data_sets{$i_num}{$file} } = map sprintf("%d", $_), @{ $data_sets{$i_num}{$file} };
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
# Tracking Matrix Collection
#######################################

sub gather_tracking_mat {
    open TRACK_IN, "$cfg{tracking_output_file}" or die;
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
    @{ $adh_lineage_props{area_shifts} }        = &gather_area_changes;
    @{ $adh_lineage_props{areas} }              = &gather_areas;
    @{ $adh_lineage_props{longevities} }        = &gather_longevities;
    @{ $adh_lineage_props{starting_edge_dist} } = &gather_starting_dist_from_edge;
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
    &output_area_changes;
    &output_adhesion_area_seqeunces(\@{ $adh_lineage_props{areas} }, \@{ $adh_lineage_props{longevities} });
    open OUTPUT, ">edge_dist.csv";
    print OUTPUT join("\n", @{ $adh_lineage_props{starting_edge_dist} });
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
