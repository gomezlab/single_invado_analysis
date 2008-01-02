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
use threads;
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

} else {
    print "\n\nGathering Data Files/Comparison Matrices from Data File\n" if $opt{debug};
    %data_sets = %{ retrieve($opt{input}) };
    &trim_data_sets;
}

if (defined $opt{output}) {
    print "\n\nStoring Data Files/Comparison Matrices\n" if $opt{debug};
    store \%data_sets, $opt{output};
}

print "\n\nReading Tracking Matrix\n" if $opt{debug};
my @tracking_mat;
&gather_tracking_matrix;

print "\n\nWriting SQL Statements\n" if $opt{debug};
&write_sql_statements;


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

    #Check for the set of variables in the config file that must be specified
    if (not defined $cfg{results_folder}) {
        die "ERROR: The location of the folder that contains the results of the focal ",
          "adhesion identification must be specified in the config file with the ", "variable \"results_folder\".";
    }

    if (not defined $cfg{exp_name}) {
        die "ERROR: The name of the experiment must be specified in the config ",
          "file with the variable \"exp_name\".";
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

    #Compute a few config variables from the provided values:
    $cfg{results_data_folder} = "$cfg{results_folder}/$cfg{exp_name}/$cfg{single_image_folder}/";
    $cfg{exp_result_folder}   = "$cfg{results_folder}/$cfg{exp_name}";

    return %cfg;
}

#######################################
# Data Set Collection
#######################################
sub gather_data_sets {
    my @folders = <$cfg{results_data_folder}/*/$cfg{raw_data_folder}>;

    my @data_files;
    push @data_files, split(/\s/, $cfg{general_data_files});
    push @data_files, split(/\s/, $cfg{lineage_analysis_data_files});
    my $image_count = 0;

    foreach my $this_folder (@folders) {
        my $i_num;
        if ($this_folder =~ /$cfg{results_data_folder}\/(.*)\/$cfg{raw_data_folder}/) {
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
        print "Data collected from ", $image_count, " images. ", "Data files gathered for each image include: ",
          join(", ", @data_files), "\n";
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

    for my $ad_num (1 .. $total_ad) {
        my @temp;
        open INPUT, "$folder/$ad_num" or die "Unable to open PixelIdxList file: $folder/$ad_num";
        for (<INPUT>) {
            chomp($_);
            my $converted_line = sprintf("%d", $_);
            if ($converted_line =~ /(\d+)/) {
                push @temp, $1;
            }
        }
        close INPUT;
        push @all_pixid, \@temp;
    }
    return @all_pixid;
}

sub trim_data_sets {
    my @excluded_nums;
    for my $ex_num (@{ $cfg{exclude_image_nums} }) {
        for my $this_num (sort { $a <=> $b } keys %data_sets) {
            if ($ex_num == $this_num) {
                delete $data_sets{$this_num};
                push @excluded_nums, $ex_num;
            }
        }
    }
    if ($opt{debug} && scalar(@excluded_nums) != 0) {
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
# Process Data Sets
#######################################
sub gather_tracking_matrix {
	open TRAC_MAT, "$cfg{exp_result_folder}/$cfg{tracking_output_file}" or die "$!";
	while (<TRAC_MAT>) {
		my @temp = split(",",$_);
		push @tracking_mat, \@temp;
	}
	close TRAC_MAT;
}


#######################################
# Process Data Sets
#######################################
sub write_sql_statements {
	&write_image_statements;
	&write_adhesion_statements;
	&write_lineage_statements;
}

sub write_image_statements {
	open STATE, ">/Users/mbergins/Documents/School\ Documents/2007/Fall/BCB Modules/Databases/final_project/database/image_insert.sql";
	for (keys %data_sets) {
		print STATE "INSERT INTO image SET exp_number = 1, image_number = $_;\n";
	}
	print STATE "INSERT INTO image SET exp_number = 1, image_number = 167, blurry = TRUE;\n";
	print STATE "INSERT INTO image SET exp_number = 1, image_number = 168, blurry = TRUE;";
	close STATE;
}

sub write_adhesion_statements {
	open STATE, ">/Users/mbergins/Documents/School\ Documents/2007/Fall/BCB Modules/Databases/final_project/database/adhesion_insert.sql";
	my $ad_count = 0;
	for my $i_num (keys %data_sets) {
		my @areas     = @{ $data_sets{$i_num}{Area} };
		my @cent_x    = @{ $data_sets{$i_num}{Centroid_x} };
		my @cent_y    = @{ $data_sets{$i_num}{Centroid_y} };
		my @cent_dist = @{ $data_sets{$i_num}{Centroid_dist_from_edge} };

		for my $i (0 .. $#areas) {
			print STATE "INSERT INTO adhesion SET image_number = $i_num, image_adhesion_number = $i, overall_adhesion_number = $ad_count, area = $areas[$i], centroid_x = $cent_x[$i], centroid_y = $cent_y[$i], centroid_dist_edge = $cent_dist[$i];\n";
			${ $data_sets{$i_num}{$i} } = $ad_count;
			$ad_count++;
		}
	}
	close STATE;
}

sub write_lineage_statements {
	open STATE, ">/Users/mbergins/Documents/School\ Documents/2007/Fall/BCB Modules/Databases/final_project/database/lineage_insert.sql";
	for my $i (0 .. $#tracking_mat) {
		my @overall_ad_list;
		my $hit = 0;
		for my $j (0 .. $#{$tracking_mat[$i]} - 1) {
			my $temp = $j + 1;

			$temp = $temp + 2 if $j > 165;
			
			my $padded_i_num = sprintf("%03d",$temp);
			
			if ($tracking_mat[$i][$j] > -1) {
				$hit = 1;
				if (not exists $data_sets{$padded_i_num}{$tracking_mat[$i][$j]}) {
				}
				else {
					push @overall_ad_list, ${$data_sets{$padded_i_num}{$tracking_mat[$i][$j]}};
				}
			}
		}
		#print "PROB\n\n" if not $hit;
		print STATE "INSERT INTO lineage SET lineage_number = $i, adhesion_number_list = \"",join(",",@overall_ad_list),"\", speed = ",rand(0.7), ";\n";
	}
	close STATE;
}
