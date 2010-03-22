#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use lib "../lib";
use lib "../lib/perl";

use strict;
use File::Temp qw/ tempfile tempdir /;
use File::Spec::Functions;
use File::Path;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use Storable;
use Statistics::Descriptive;
use Statistics::Distributions;
use POSIX;

use Config::Adhesions qw(ParseConfig);
use Image::Data::Collection;
use Text::CSV::Simple::Extra;
use Emerald;
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|config=s", "debug|d", "lsf|l") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my %cfg = ParseConfig(\%opt);

###############################################################################
# Main Program
###############################################################################

if ($opt{lsf}) {
    #$0 - the name of the program currently running, used to protect against
    #future file name changes
    my @command = "$0 -cfg $opt{cfg}";
    $command[0] .= " -skip_lin_regions" if $opt{skip_lin_regions};
    $opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'track_analysis');
    $opt{resource} = "mem32";
    $opt{queue} = "week";
    if (defined $cfg{job_group}) {
        $opt{job_group} = $cfg{job_group};
    }

    &FA_job::send_general_lsf_program(\@command,\%opt);
    
    exit;
}

print "Gathering/Converting Data Files\n" if $opt{debug};
my @data_files;
push @data_files, @{ $cfg{general_data_files} };
push @data_files, @{ $cfg{lineage_analysis_data_files} };
my %data_sets = &Image::Data::Collection::gather_data_sets(\%cfg, \%opt, \@data_files);
%data_sets = &convert_data_to_units(\%data_sets, \%cfg);
my @available_data_types = &gather_data_types;
if ($opt{debug}) {
	print "Data types detected:\n";
	print join("\n", @available_data_types);
	print "\n";
}

print "\n\nCollecting Tracking Matrix\n" if $opt{debug};
my @tracking_mat = &Image::Data::Collection::read_in_tracking_mat(\%cfg, \%opt);

print "\n\nCreating/Outputing Overall Cell Property Files\n" if $opt{debug};
&gather_and_output_overall_cell_properties;

print "\n\nCreating/Outputing Individual Adhesion Property Files\n" if $opt{debug};
my @single_ad_props = &gather_single_ad_props(\%cfg, \%opt);
&output_single_adhesion_props(@single_ad_props);
@single_ad_props = ();
undef @single_ad_props;

print "\n\nCreating/Outputing Adhesion Lineage Property Files\n", if $opt{debug};
&gather_and_output_lineage_properties;

###############################################################################
# Functions
###############################################################################

#######################################
#Raw Data Gathering/Converting
#######################################
sub convert_data_to_units {
    my %data_sets = %{ $_[0] };
    my %cfg       = %{ $_[1] };

    my $lin_conv_factor = $cfg{pixel_size};
    my $sq_conv_factor  = $lin_conv_factor**2;
    my @no_conversion =
	  qw(Class Centroid_x Centroid_y Eccentricity Solidity
	  Background_corrected_signal Angle_to_center Orientation
	  Shrunk_corrected_signal Cell_mean_intensity Outside_mean_intensity
	  Cell_not_ad_mean_intensity Adhesion_mean_intensity CB_corrected_signal
	  Global_gel_diff Local_gel_diff End_local_gel_diff Local_gel_diff_corr 
	  First_local_gel_diff);

    for my $time (sort keys %data_sets) {
        for my $data_type (keys %{ $data_sets{$time} }) {
            if (grep $data_type eq $_, qw(Centroid_dist_from_edge Centroid_dist_from_center MajorAxisLength 
                                          MinorAxisLength)) {
                @{ $data_sets{$time}{$data_type} } = map $lin_conv_factor * $_, @{ $data_sets{$time}{$data_type} };
            } elsif (grep $data_type eq $_, qw(Area Cell_size)) {
                @{ $data_sets{$time}{$data_type} } = map $sq_conv_factor * $_, @{ $data_sets{$time}{$data_type} };
            } elsif ($data_type eq "Edge_speed" || $data_type eq "Edge_projection") {
                foreach (0 .. $#{ $data_sets{$time}{$data_type} }) {
                    @{ $data_sets{$time}{$data_type}[$_] } = map $lin_conv_factor * $_,
                      @{ $data_sets{$time}{$data_type}[$_] };
                }
            } elsif ((grep $data_type eq $_, @no_conversion)
                || ($data_type =~ /adhesion_signal/)) {

                #This is the arbitrary units place, don't do any unit
                #conversion
            } else {
                die "Unable to determine the appropriate conversion factor (if any) for $data_type.";
            }
        }
    }
    return %data_sets;
}

sub gather_data_types {
    my $first_key = (sort keys %data_sets)[0];
    return (keys %{ $data_sets{$first_key} });
}

#######################################
#Single Adhesion Props
#######################################
sub gather_single_ad_props {
    my @data;

    my @possible_data_types =
      qw(Area Average_adhesion_signal Eccentricity Solidity MajorAxisLength 
         MinorAxisLength Centroid_dist_from_edge Centroid_dist_from_center 
         Variance_adhesion_signal);

    my @single_ad_data_types = map {
        my $type = $_;
        grep $type eq $_, @available_data_types;
    } @possible_data_types;

    my @i_num_list = sort keys %data_sets;
    my %i_num_to_col = map { $i_num_list[$_] => $_ } (0 .. $#i_num_list);

    my @first_line = qw(I_num ad_num);
    push @first_line, @single_ad_data_types;
    push @data,       \@first_line;

    foreach my $i_num (sort keys %data_sets) {
        my $col                = $i_num_to_col{$i_num};
        my @ad_nums_to_include = map {
            if ($tracking_mat[$_][$col] >= 0) {
                $tracking_mat[$_][$col];
            } else {
                ();
            }
        } (0 .. $#tracking_mat);

        @ad_nums_to_include = sort { $a <=> $b } @ad_nums_to_include;

        foreach my $ad_num (@ad_nums_to_include) {
            my @line = ($i_num, $ad_num);
            for my $i (0 .. $#single_ad_data_types) {
                push @line, $data_sets{$i_num}{ $single_ad_data_types[$i] }[$ad_num];
            }
            push @data, \@line;
        }
    }
    return @data;
}

sub output_single_adhesion_props {
    my @single_ad_props = @_;
    mkpath(catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}));

    my $output_file =
      catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{individual_adhesions_props_file});
    &output_mat_csv(\@single_ad_props, $output_file);
}

#######################################
#Overall Props
#######################################
sub gather_and_output_overall_cell_properties {
	my @single_props = qw(Cell_size Cell_mean_intensity Outside_mean_intensity
	Cell_not_ad_mean_intensity Adhesion_mean_intensity);

	foreach my $data_type (@single_props) {
        next if (not(grep $data_type eq $_, @available_data_types));
		
        my @prop_sequence = &gather_single_number_time_series("$data_type");
        my $output_file = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder},"single_props", $data_type . ".csv");
        &output_mat_csv(\@prop_sequence, $output_file);
	}
    
    my @total_ad_size = &gather_total_adhesion_size_time_series;
    my $output_file = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, "single_props", "Ad_size.csv");
    &output_mat_csv(\@total_ad_size, $output_file);
}

sub gather_single_number_time_series {
    my $prop = shift @_;

    my @sequence;
    for my $i_num (sort keys %data_sets) {
        die "Greater than one entry in the Cell_size field in image number $i_num" 
          if scalar(@{$data_sets{$i_num}{$prop}}) > 1;
        
        push @sequence, $data_sets{$i_num}{$prop}[0];
    }
    return @sequence;
}

sub gather_total_adhesion_size_time_series {
    my @total_ad_size;
    for my $i_num (sort keys %data_sets) {
        my $this_total_size = 0;
        for my $this_size (@{$data_sets{$i_num}{Area}}) {
            $this_total_size += $this_size;
        }
        push @total_ad_size, $this_total_size;
    }
    return @total_ad_size;
}

#######################################
#Adhesion Lineage Property Collection
#######################################
sub gather_and_output_lineage_properties {
    mkpath(catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}));
    mkpath(catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{lineage_ts_folder}));

    my %props;
    
    if (grep $_ eq "Edge_speed", @available_data_types) {
        my %edge_data = &gather_edge_velo_data("Edge_speed");
        my $base_dir = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, "Edge_speed");
        &output_mat_csv($edge_data{pre_birth}, catfile($base_dir, "pre_birth.csv"));
        &output_mat_csv($edge_data{post_birth}, catfile($base_dir, "post_birth.csv"));
        &output_mat_csv($edge_data{pre_death}, catfile($base_dir, "pre_death.csv"));
        &output_mat_csv($edge_data{post_death}, catfile($base_dir, "post_death.csv"));
        &output_mat_csv($edge_data{null_data}, catfile($base_dir, "null.csv"));
    }
    if (grep $_ eq "Edge_projection", @available_data_types) {
        my %edge_data = &gather_edge_velo_data("Edge_projection");
        my $base_dir = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, "Edge_projection");
        &output_mat_csv($edge_data{pre_birth}, catfile($base_dir, "pre_birth.csv"));
        &output_mat_csv($edge_data{post_birth}, catfile($base_dir, "post_birth.csv"));
        &output_mat_csv($edge_data{pre_death}, catfile($base_dir, "pre_death.csv"));
        &output_mat_csv($edge_data{post_death}, catfile($base_dir, "post_death.csv"));
        &output_mat_csv($edge_data{null_data}, catfile($base_dir, "null.csv"));
    }

    #Pure Time Series Props
	my @ts_props = qw(Area Global_gel_diff Local_gel_diff Local_gel_diff_corr
		First_local_gel_diff);
    foreach my $data_type (@ts_props) {
        next if (not(grep $data_type eq $_, @available_data_types));

        my $this_data = &gather_prop_seq($data_type);
        &output_prop_time_series($this_data, $data_type);
    }

    $props{longevity}               = &gather_longevities;
    $props{uncertain_longevity}     = &gather_uncertain_longevities;
    $props{merge_count}             = &gather_merge_count;
    $props{death_status}            = &gather_death_status;
    $props{split_birth_status}      = &gather_split_birth_status;
    $props{average_eccentricity}    = &gather_average_eccentricity;
    $props{Average_adhesion_signal} = &gather_prop_seq("Average_adhesion_signal");
    &output_prop_time_series($props{Average_adhesion_signal}, "Average_adhesion_signal");
    $props{ad_sig} = &gather_average_value($props{Average_adhesion_signal});
    undef $props{Average_adhesion_signal};
    	
    $props{Centroid_x} = &gather_prop_seq("Centroid_x");
    $props{start_x} = &gather_first_entry($props{Centroid_x});
    $props{end_x} = &gather_last_entry($props{Centroid_x});
    $props{birth_i_num} = &gather_birth_i_num($props{Centroid_x});
    $props{death_i_num} = &gather_death_i_num($props{Centroid_x});
    undef $props{Centroid_x};

    $props{Centroid_y} = &gather_prop_seq("Centroid_y");
    $props{start_y} = &gather_first_entry($props{Centroid_y});
    $props{end_y} = &gather_last_entry($props{Centroid_y});
    undef $props{Centroid_y};
	
    ($props{speeds}{All}, $props{velocity}) = &gather_adhesion_speeds;
    ($props{average_speeds}, $props{variance_speeds}, $props{max_speeds}) = 
		&gather_speed_props($props{speeds}{All});

    if (grep "Area" eq $_, @available_data_types) {
        $props{Area} = &gather_prop_seq("Area");
        $props{largest_area} = &gather_largest_entry($props{Area});
        $props{largest_area_inum} = &gather_largest_entry_inum($props{Area});
        $props{mean_area} = &gather_average_value($props{Area});
        undef $props{Area};
    }
	
	$props{End_local_gel_diff} = &gather_prop_seq("End_local_gel_diff");
	$props{last_local_gel_diff} = &gather_entry_by_inum($props{End_local_gel_diff}, $props{largest_area_inum});
	undef $props{End_local_gel_diff};

    if (grep "Centroid_dist_from_center" eq $_, @available_data_types) {
        $props{Centroid_dist_from_center} = &gather_prop_seq("Centroid_dist_from_center");
        $props{starting_center_dist} = &gather_first_entry($props{Centroid_dist_from_center});
        $props{ending_center_dist}   = &gather_last_entry($props{Centroid_dist_from_center});
        undef $props{Centroid_dist_from_center};
    }

    if (grep "Centroid_dist_from_edge" eq $_, @available_data_types) {
        $props{Centroid_dist_from_edge} = &gather_prop_seq("Centroid_dist_from_edge");
        $props{starting_edge_dist} = &gather_first_entry($props{Centroid_dist_from_edge});
        $props{ending_edge_dist}   = &gather_last_entry($props{Centroid_dist_from_edge});
        undef $props{Centroid_dist_from_edge};
    }

    my @lin_summary_data = &gather_lineage_summary_data(\%props);
    my $output_file = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{lineage_summary_props_file});
    &output_mat_csv(\@lin_summary_data, $output_file);
    %props = ();

	print "\n\n" if $opt{debug};
}

sub gather_edge_velo_data {
    print "\r", " " x 80, "\rGathering Edge Velocities Summaries" if $opt{debug};
    
    my $edge_type = "Edge_speed";
    if (scalar(@_) > 0) {
        $edge_type = $_[0];
    }

    my $default_val = "NA";
    my $min_living_time = 10;

    my @pre_birth_data;
    my @post_birth_data;
    my @pre_death_data;
    my @post_death_data;
    my @null_data;
    my @data_keys = sort keys %data_sets;
    
    #Cycle through each row of the tracking matrix, extracting all the relavent
    #edge speed data points
    for my $i (0 .. $#tracking_mat) {
        my $first_data_index = (grep $tracking_mat[$i][$_] >= 0, (0 .. $#{ $tracking_mat[$i] }))[0];
        my $last_data_index  = (grep $tracking_mat[$i][$_] >= 0, (0 .. $#{ $tracking_mat[$i] }))[-1];
        
        #Data structure notes:
        #   -the $edge_type variables hold scalar projections between the
        #   current position of the adhesion in any given time points and the
        #   closest edge tracking point in EVERY image, so to get all the points
        #   prior to birth we need only analyze the $edge_type data from the
        #   birth image
        
        #Collect the null data from this row
        for my $this_index ($first_data_index .. $last_data_index) {
            my $this_i_num = $data_keys[$this_index];
            my $this_i_num_index = grep $data_keys[$_] == $this_i_num, (0 .. $#data_keys);
            my $this_ad_num = $tracking_mat[$i][$this_index];
            my @velo_data = @{ $data_sets{$this_i_num}{$edge_type} };
            die if not exists $velo_data[$this_ad_num][$this_i_num_index];
            push @{$null_data[$i]}, $velo_data[$this_ad_num][$this_i_num_index];
        }

        #Determine upto birth data points
        my $birth_ad_num = $tracking_mat[$i][$first_data_index];
        die "First pre-birth adhesion number is not valid ($birth_ad_num)." if $birth_ad_num < 0;
        my $birth_i_num     = $data_keys[$first_data_index];
        my @birth_velo_data = @{ $data_sets{$birth_i_num}{$edge_type} };
        my @this_pre_birth  = @{ $birth_velo_data[$birth_ad_num] }[ 0 .. ($first_data_index - 1) ];
        push @pre_birth_data, \@this_pre_birth;
        
        #Determine post birth and pre death data points
        if ($last_data_index - $first_data_index >= $min_living_time) {
            my @tracking_mat_indexes = ($first_data_index .. $last_data_index);
            die if scalar(@tracking_mat_indexes) < $min_living_time;
            
            for my $this_index (0 .. $#tracking_mat_indexes) {
                my $this_data_index = $tracking_mat_indexes[$this_index];
                my $this_i_num = $data_keys[$this_data_index];
                my $this_i_num_index = grep $data_keys[$_] == $this_i_num, (0 .. $#data_keys);
                my $this_ad_num = $tracking_mat[$i][$this_data_index];
                my @velo_data = @{ $data_sets{$this_i_num}{$edge_type} };

                if ($this_index <= floor($#tracking_mat_indexes/2)) {
                    push @{$post_birth_data[$i]}, $velo_data[$this_ad_num][$this_i_num_index];
                } else {
                    push @{$pre_death_data[$i]}, $velo_data[$this_ad_num][$this_i_num_index];
                }
            }
        } else {
            push @post_birth_data, [];
            push @pre_death_data, [];
        }

        #Determine post death data points
        my $death_ad_num = $tracking_mat[$i][$last_data_index];
        die "Death adhesion number is not valid ($death_ad_num)." if $death_ad_num < 0;
        my $death_i_num     = $data_keys[$last_data_index];
        my @death_velo_data = @{ $data_sets{$death_i_num}{$edge_type} };
        my @this_post_death = @{ $death_velo_data[$death_ad_num] }[ $last_data_index + 1 .. $#{ $death_velo_data[$death_ad_num] } ];
        push @post_death_data, \@this_post_death;
    }

    @null_data = &pad_arrays_to_longest(\@null_data, "push");
    @pre_birth_data = &pad_arrays_to_longest(\@pre_birth_data, "unshift");
    @post_birth_data = &pad_arrays_to_longest(\@post_birth_data, "push");
    @pre_death_data = &pad_arrays_to_longest(\@pre_death_data, "unshift");
    @post_death_data = &pad_arrays_to_longest(\@post_death_data, "push");
    
    return ( pre_birth => \@pre_birth_data, 
             post_birth => \@post_birth_data,
             pre_death => \@pre_death_data,
             post_death => \@post_death_data,
             null_data => \@null_data,
           );
}

sub find_longest_row {
    my @data = @_;

    my @seq_lengths = map scalar(@{ $data[$_] }), (0 .. $#data);
    my $max_length = 0;
    my $row;
    foreach (0 .. $#seq_lengths) {
        my $length = $seq_lengths[$_];
        $row = $_ if $length > $max_length;
        $max_length = $length if $length > $max_length;
    }
    return $max_length;
}

sub pad_arrays_to_longest {
    my @data = @{$_[0]};
    my $add_style = $_[1];
    my $default_val = "NA";
    if (scalar(@_) >= 3) {
        $default_val = ${$_[2]};
    }

    if (!($add_style eq "unshift" || $add_style eq "push")) {
        die "Second parameter to pad_arrays_to_longest must be either unshift or push, got $add_style.\n";
    }

    my $max_length = &find_longest_row(@data);
    for my $i (0 .. $#data) {
        for (1 .. ($max_length - scalar(@{ $data[$i] }))) {
            if ($add_style eq "unshift") {
                unshift @{ $data[$i] }, $default_val;
            } else {
                push @{ $data[$i] }, $default_val;
            }
        }
        die "Padding to longest length failed in post-death edge velocity, got ",
          scalar(@{ $data[$i] }), " expected $max_length" if scalar(@{ $data[$i] }) != $max_length;
    }

    return @data;
}

sub output_prop_time_series {
    my @data      = @{ $_[0] };
    my $file_name = $_[1];
    my $output_file =
      catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{lineage_ts_folder}, $file_name . ".csv");
    &output_mat_csv(\@data, $output_file);
}

sub gather_longevities {
    my @longevities;
    my $default_val = "NA";
    $default_val = $_[0] if (scalar(@_) > 0);

    print "\r", " " x 80, "\rGathering Longevity" if $opt{debug};
    my $default_count = 0;
    for my $i (0 .. $#tracking_mat) {
        if ($tracking_mat[$i][0] > -1 || $tracking_mat[$i][-1] > -1) {
            $default_count++;
            push @longevities, $default_val;
        } else {
            my $count = scalar(grep $tracking_mat[$i][$_] > -1, 0..$#{$tracking_mat[$i]});
            push @longevities, $count;
        }
    }
    return \@longevities;
}

sub gather_uncertain_longevities {
    my @longevities;
    my $default_val = "NA";
    $default_val = $_[0] if (scalar(@_) > 0);

    print "\r", " " x 80, "\rGathering Uncertain Longevity" if $opt{debug};
    my $default_count = 0;
    for my $i (0 .. $#tracking_mat) {
		my $count = scalar(grep $tracking_mat[$i][$_] > -1, 0..$#{$tracking_mat[$i]});
		push @longevities, $count;
    }
    return \@longevities;
}

sub gather_prop_seq {
    my $prop        = $_[0];
    my $default_val = "NaN";
    $default_val = $_[1] if (scalar(@_) > 1);

    print "\r", " " x 80, "\rGathering Prop Seq: $prop" if $opt{debug};
    my @prop_vals;
    my @data_keys = sort keys %data_sets;
    for my $i (0 .. $#tracking_mat) {
        for my $j (0 .. $#{ $tracking_mat[$i] }) {
            my $ad_num = $tracking_mat[$i][$j];

            if ($ad_num <= -1) {
                push @{ $prop_vals[$i] }, $default_val;
                next;
            }

            my $i_num = $data_keys[$j];
            if (not defined ${ $data_sets{$i_num}{$prop} }[$ad_num]) {
                die "Can't find value for property ($prop) in image number ($i_num), ad number ($ad_num)\n"
                  . "Data keys in this I number:\n"
                  . join(" ", keys %{ $data_sets{$i_num} });
            }

            push @{ $prop_vals[$i] }, ${ $data_sets{$i_num}{$prop} }[$ad_num];
        }
    }
    return \@prop_vals;
}

sub gather_largest_entry {
    my @data = @{ $_[0] };

    print "\r", " " x 80, "\rGathering Largest Entries" if $opt{debug};
    my @largest_data;
    for my $i (0 .. $#data) {
        my $largest = 0;
        for my $j (0 .. $#{ $data[$i] }) {
            next if ($data[$i][$j] eq "NaN");
            $largest = $data[$i][$j] if ($largest < $data[$i][$j]);
        }
        push @largest_data, $largest;
    }
    return \@largest_data;
}

sub gather_largest_entry_inum {
	my @data = @{ $_[0] };

	print "\r", " " x 80, "\rGathering Largest Entry Image Numbers" if $opt{debug};
	my @largest_inum;
	for my $i (0 .. $#data) {
		my $largest = 0;
		my $largest_inum = -1;
		for my $j (0 .. $#{ $data[$i] }) {
			next if ($data[$i][$j] eq "NaN");
			if ($largest < $data[$i][$j]) {
					$largest = $data[$i][$j];
					$largest_inum = $j;
			}
		}
		die if ($largest_inum == -1);
		push @largest_inum, $largest_inum;
	}
	return \@largest_inum;
}

sub gather_entry_by_inum {
	my @data = @{ $_[0] };
	my @inum_data = @{ $_[1] };

	my @entries_by_inum;
	for my $i (0 .. $#data) {
		push @entries_by_inum, $data[$i][$inum_data[$i]];
	}
	return \@entries_by_inum;
}

sub gather_birth_i_num {
    my @data = @{ $_[0] };

    print "\r", " " x 80, "\rGathering Birth Image Number" if $opt{debug};
    my @starting_data;
    for my $i (0 .. $#data) {
        my $first_data_index = (grep $data[$i][$_] ne "NaN", (0 .. $#{ $data[$i] }))[0];

        if ($first_data_index == 0) {
            $starting_data[$i] = "NA";
        } else {
            $starting_data[$i] = $first_data_index;
        }
    }
    return \@starting_data;
}

sub gather_death_i_num {
    my @data = @{ $_[0] };

    print "\r", " " x 80, "\rGathering Death Image Number" if $opt{debug};
    my @ending_data;

    my $max_index = $#{$data[0]};
    for my $i (0 .. $#data) {
        my $last_data_index = (grep $data[$i][$_] ne "NaN", (0 .. $#{ $data[$i] }))[-1];
        if ($last_data_index == $max_index) {
            $ending_data[$i] = "NA";
        } else {
            $ending_data[$i] = $last_data_index;
        }
    }
    return \@ending_data;
}

sub gather_first_entry {
    my @data = @{ $_[0] };

    print "\r", " " x 80, "\rGathering First Entries" if $opt{debug};
    my @starting_data;
    for my $i (0 .. $#data) {
        my $first_data_index = (grep $data[$i][$_] ne "NaN", (0 .. $#{ $data[$i] }))[0];
        $starting_data[$i] = $data[$i][$first_data_index];
    }
    return \@starting_data;
}

sub gather_last_entry {
    my @data = @{ $_[0] };

    print "\r", " " x 80, "\rGathering Last Entries" if $opt{debug};
    my @last_data;
    for my $i (0 .. $#data) {
        my $last_data_index = (grep $data[$i][$_] ne "NaN", (0 .. $#{ $data[$i] }))[-1];
        $last_data[$i] = $data[$i][$last_data_index];
    }
    return \@last_data;
}

sub gather_average_value {
    my @signal = @{ $_[0] };

    print "\r", " " x 80, "\rGathering Average Values" if $opt{debug};
    my @average_sig;
    for my $i (0 .. $#signal) {
        my $stat = Statistics::Descriptive::Full->new();
        for my $j (0 .. $#{ $signal[$i] }) {
            $stat->add_data($signal[$i][$j]) if ($signal[$i][$j] ne "NaN" && $signal[$i][$j] ne "nan");
        }

        push @average_sig, $stat->mean();
    }
    return \@average_sig;
}

sub gather_adhesion_speeds {
    my @speed;
    my @velocity;
    my @data_keys = sort keys %data_sets;
    print "\r", " " x 80, "\rGathering Adhesion Speeds" if $opt{debug};
    for my $i (0 .. $#tracking_mat) {
        for my $j (0 .. $#{ $tracking_mat[$i] } - 1) {
            my $start_i_num = $data_keys[$j];
            my $end_i_num   = $data_keys[ $j + 1 ];

            my $start_ad_num = $tracking_mat[$i][$j];
            my $end_ad_num   = $tracking_mat[$i][ $j + 1 ];

            if ($start_ad_num > -1 && $end_ad_num > -1) {
                my $start_x = ${ $data_sets{$start_i_num}{Centroid_x} }[$start_ad_num];
                my $start_y = ${ $data_sets{$start_i_num}{Centroid_y} }[$start_ad_num];

                my $end_x = ${ $data_sets{$end_i_num}{Centroid_x} }[$end_ad_num];
                my $end_y = ${ $data_sets{$end_i_num}{Centroid_y} }[$end_ad_num];

                my $speed = sqrt(($start_x - $end_x)**2 + ($start_y - $end_y)**2);

                push @{ $speed[$i] }, $speed;
                push @{ $velocity[$i] }, [ ($start_x - $end_x, $start_y - $end_y) ];
            } else {
                push @{ $speed[$i] },    "NaN";
                push @{ $velocity[$i] }, "NaN";
            }
        }
    }
    return (\@speed, \@velocity);
}

sub gather_speed_props {
    my @speed = @{ $_[0] };

    my @av_speeds;
    my @var_speeds;
    my @max_speeds;

    print "\r", " " x 80, "\rGathering Adhesion Speed Properties" if $opt{debug};
    for my $i (0 .. $#speed) {
        my $stat = Statistics::Descriptive::Full->new();

        my @ad_speeds = grep $_ ne "NaN", @{ $speed[$i] };

        $stat->add_data(@ad_speeds) if (scalar(@ad_speeds) != 0);

        if ($stat->count() > 0) {
            push @av_speeds,  $stat->mean();
            push @var_speeds, $stat->variance();
            push @max_speeds, $stat->max();
        } else {
            push @av_speeds,  "NaN";
            push @var_speeds, "NaN";
            push @max_speeds, "NaN";
        }
    }
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(grep $_ ne "NaN", @av_speeds);

    return \@av_speeds, \@var_speeds, \@max_speeds;
}

sub gather_merge_count {
    my @merge_count = map 0, (0 .. $#tracking_mat);
    my @data_keys = sort keys %data_sets;
    print "\r", " " x 80, "\rGathering Merge Counts" if $opt{debug};
    for my $i (0 .. $#tracking_mat) {
        for my $j (0 .. $#{ $tracking_mat[$i] }) {
            next if ($tracking_mat[$i][$j] <= -1);
            for my $k (0 .. $#{ $tracking_mat[$i] }) {
                $merge_count[$i]++ if ($tracking_mat[$i][$k] == ($tracking_mat[$i][$j] + 2) * -1);
            }
        }
    }

    return \@merge_count;
}

sub gather_death_status {
    my @death_status = map 0, (0 .. $#tracking_mat);
    print "\r", " " x 80, "\rGathering Death Status" if $opt{debug};
    for my $i (0 .. $#tracking_mat) {
        my $pre_birth_num =
          (grep $tracking_mat[$i][$_] >= 0 && $tracking_mat[$i][ $_ + 1 ] == -1, (0 .. $#{ $tracking_mat[$i] } - 1))[0];

        if (defined $pre_birth_num) {
            $death_status[$i] = 1;
        }
    }

    return \@death_status;
}

sub gather_average_eccentricity {
    print "\r", " " x 80, "\rGathering Average Eccentricity" if $opt{debug};
	my @major_axis = @{&gather_prop_seq("MajorAxisLength")};
	my @minor_axis = @{&gather_prop_seq("MinorAxisLength")};
	
	my @eccentricity;

    for my $i (0 .. $#major_axis) {
        for my $j (0 .. $#{ $major_axis[$i] }) {
			die "Minor Axis value of zero at position: $i $j, value of $minor_axis[$i][$j]" 
				if ($minor_axis[$i][$j] == 0);
			$eccentricity[$i][$j] = $major_axis[$i][$j]/$minor_axis[$i][$j];
        }
    }
	
	my @average_eccentricity = @{&gather_average_value(\@eccentricity)};

	return \@average_eccentricity;
}

sub gather_split_birth_status {
    my @sb_status = map 0, (0 .. $#tracking_mat);
    print "\r", " " x 80, "\rGathering Spilt Birth Status" if $opt{debug};
    for my $i (0 .. $#tracking_mat) {
        my $pre_birth_num =
          (grep $tracking_mat[$i][$_] >= 0 && $tracking_mat[$i][ $_ - 1 ] <= -2, (1 .. $#{ $tracking_mat[$i] }))[0];

        if (defined $pre_birth_num) {
            $sb_status[$i] = 1;
        }
    }

    return \@sb_status;
}

sub gather_lineage_summary_data {
    my %props          = %{ $_[0] };
	my @possible_props = qw(longevity uncertain_longevity largest_area mean_area
		largest_area_inum starting_edge_dist ending_edge_dist
		starting_center_dist ending_center_dist merge_count death_status
		split_birth_status average_speeds max_speeds ad_sig birth_i_num start_x
		start_y death_i_num end_x end_y average_eccentricity last_local_gel_diff);

    my @lin_summary_data;
    for (@possible_props) {
        my $this_prop = $_;
        next if (not(grep $this_prop eq $_, keys %props));

        push @{ $lin_summary_data[0] }, $this_prop;

        for my $i (0 .. $#{ $props{$this_prop} }) {
            push @{ $lin_summary_data[ $i + 1 ] }, $props{$this_prop}[$i];
        }
    }

    return @lin_summary_data;
}

################################################################################
#Documentation
################################################################################

=head1 NAME

gather_tracking_results.pl - Gather the properties of the FA found in prior
steps

=head1 SYNOPSIS

gather_tracking_results.pl -cfg FA_config

=head1 Description

This program gathers the properties of the detected FA. The sets of properties
can be roughly divided into two types: static and dynamic. The static properties
are those that only concern each image FA individually. Such properties include
the size, area and intensity of fluorescence signal. The dynamic properties
include each of the static properties tracked through time. Most of the results
are saved to csv files. 

Required parameter(s):

=over 

=item * cfg or c: the focal adhesion analysis config file

=back

Optional parameter(s):

=over 

=item * tracking_mat: specifies a file name to use instead of the file specified
in the config file for the tracking matrix


=item * skip_lin_regions: debugging command that skips the time consuming
production of the pixel properties results

=item * debug or d: print debuging information during program execution

=item * emerald: setups and runs a job tailored for the LSF job system on emerald

=back

=head1 EXAMPLES

gather_tracking_results.pl -cfg FA_config

=head1 AUTHORS

Matthew Berginski (mbergins@unc.edu)

Documentation last updated: 6/6/2008 

=cut
