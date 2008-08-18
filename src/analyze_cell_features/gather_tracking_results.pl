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
use Getopt::Long;
use Data::Dumper;
use Storable;
use Statistics::Descriptive;
use Statistics::Distributions;

use Config::Adhesions qw(ParseConfig);
use Image::Data::Collection;
use Text::CSV::Simple::Extra;
use Math::R;
use Emerald;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|config=s", "debug|d", "tracking_mat=s", "skip_lin_regions", "emerald|e") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my %cfg = ParseConfig(\%opt);

if (exists $opt{tracking_mat} && not -e catfile($cfg{exp_results_folder}, $cfg{tracking_folder}, $opt{tracking_mat})) {
    die "Could not find the tracking mat specified on the command line, looked in: "
      ,catfile($cfg{exp_results_folder}, $cfg{tracking_folder}, $opt{tracking_mat});
} else {
    my $all_but_type_tracking = '';
    if ($opt{tracking_mat} =~ m/(.*)\.(.*?)/) {
        $all_but_type_tracking = $1;
    } else {
        die "Could not find the proper end of the tracking matrix without file type";
    }

    $cfg{adhesion_props_folder} = catdir($cfg{adhesion_props_folder}, $all_but_type_tracking);
    $cfg{tracking_output_file} = $opt{tracking_mat};
}

###############################################################################
# Main Program
###############################################################################

if ($opt{emerald}) {
    my $error_folder = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'track_analysis');
    mkpath($error_folder);
    
    my %emerald_opt = ("folder" => $error_folder);
    my @command = "$0 -cfg $opt{cfg}";
    $command[0] .= " -tracking_mat $opt{tracking_mat}" if (exists $opt{tracking_mat});

    @command = &Emerald::create_general_LSF_commands(\@command,\%emerald_opt);
    &Emerald::send_LSF_commands(\@command);
    exit;
}

print "Gathering/Converting Data Files\n" if $opt{debug};

my @data_files;
push @data_files, @{ $cfg{general_data_files} };
push @data_files, @{ $cfg{lineage_analysis_data_files} };
my %data_sets = &Image::Data::Collection::gather_data_sets(\%cfg, \%opt, \@data_files);
%data_sets = &Image::Data::Collection::trim_data_sets(\%cfg, \%opt, \%data_sets);
%data_sets = &convert_data_to_units(\%data_sets, \%cfg);
my @available_data_types = &gather_data_types;

print "\n\nCollecting Tracking Matrix\n" if $opt{debug};
my @tracking_mat = &Image::Data::Collection::read_in_tracking_mat(\%cfg, \%opt);

print "\n\nCreating Individual Adhesion Properties Plots\n" if $opt{debug};
my @single_ad_props = &gather_single_ad_props(\%cfg, \%opt);
&output_single_adhesion_props;
#&build_single_ad_plots;

print "\n\nCreating Adhesion Lineage Property Plots\n", if $opt{debug};
my %ad_lineage_props = &gather_ad_lineage_properties;
&output_adhesion_lineage_props;
#&build_lineage_plots;

print "\n\nBuilding R Model Files\n", if $opt{debug};
if (not($opt{skip_lin_regions})) {
    &run_R_linear_region_code;
}

###############################################################################
# Functions
###############################################################################
sub convert_data_to_units {
    my %data_sets = %{ $_[0] };
    my %cfg       = %{ $_[1] };

    my $lin_conv_factor = $cfg{pixel_size} / $cfg{target_unit_size};
    my $sq_conv_factor  = $lin_conv_factor**2;

    for my $time (keys %data_sets) {
        for my $data_type (keys %{ $data_sets{$time} }) {
            if (grep $data_type eq $_, qw(Centroid_x Centroid_y Centroid_dist_from_edge Centroid_dist_from_center)) {
                @{ $data_sets{$time}{$data_type} } = map $lin_conv_factor * $_, @{ $data_sets{$time}{$data_type} };
            } elsif (grep $data_type eq $_, qw(Area Cell_size)) {
                @{ $data_sets{$time}{$data_type} } = map $sq_conv_factor * $_, @{ $data_sets{$time}{$data_type} };
            } elsif ((grep $data_type eq $_, qw(Class Eccentricity Solidity Background_corrected_signal))
                || ($data_type =~ /adhesion_signal/)) {

                #This is the arbitrary units place, don't do any unit
                #conversion
            } else {
                die "Unable to determine the appropriate conversion factor for $data_type.";
            }
        }
    }
    return %data_sets;
}

sub gather_data_types {
    my $first_key = (keys %data_sets)[0];
    return (keys %{ $data_sets{$first_key} });
}

#######################################
#Single Adhesion Props
#######################################
sub gather_single_ad_props {
    my @data;

    my @possible_data_types =
      qw(Area Average_adhesion_signal Eccentricity Solidity Centroid_dist_from_edge Centroid_dist_from_center Variance_adhesion_signal);

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
    mkpath(catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}));

    my $output_file =
      catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{individual_adhesions_props_file});
    &output_mat_csv(\@single_ad_props, $output_file);
}

sub build_single_ad_plots {
    my @r_code;

    my $data_dir = catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder});
    my $plot_dir = catdir($data_dir, $cfg{plot_folder});

    my $xy_default = "pch=19,cex=0.4";

    my $pdf_default       = "height=12, width=12, pointsize=24";
    my @possible_xy_plots = (
        {
            xy        => "adhesions\$Area,adhesions\$Average_adhesion_signal",
            main      => "",
            xlab      => "expression(paste('Area (', mu, m^2, ')'))",
            ylab      => "'Paxillin Concentration (AU)'",
            file_name => "area_vs_pax.pdf",
            plot_para => $xy_default,
            pdf_para  => $pdf_default,
        },
        {
            xy        => "adhesions\$Centroid_dist_from_edge,adhesions\$Centroid_dist_from_center",
            main      => "",
            xlab      => "expression(paste('Distance from Edge (', mu, 'm)'))",
            ylab      => "expression(paste('Distance from Center (', mu, 'm)'))",
            file_name => "center_dist_vs_edge_dist.pdf",
            plot_para => $xy_default,
            pdf_para  => $pdf_default,
        },
        {
            xy        => "adhesions\$Area,adhesions\$Centroid_dist_from_edge",
            main      => "",
            xlab      => "expression(paste('Area (', mu, m^2, ')'))",
            ylab      => "expression(paste('Distance from Edge (', mu, 'm)'))",
            file_name => "area_vs_dist.pdf",
            plot_para => $xy_default,
            pdf_para  => $pdf_default,
        },
        {
            xy        => "adhesions\$Average_adhesion_signal,adhesions\$Centroid_dist_from_edge",
            main      => "",
            xlab      => "'Paxillin Concentration (AU)'",
            ylab      => "expression(paste('Distance from Edge (', mu, 'm)'))",
            file_name => "sig_vs_dist.pdf",
            plot_para => $xy_default,
            pdf_para  => $pdf_default,
        },
        {
            xy        => "adhesions\$Centroid_dist_from_edge,adhesions\$Ecc",
            main      => "",
            xlab      => "expression(paste('Distance from Edge (', mu, 'm)'))",
            ylab      => "Eccentricity",
            file_name => "dist_vs_eccen.pdf",
            plot_para => $xy_default,
            pdf_para  => $pdf_default,
        },
    );

    my @prop_types = @{ $single_ad_props[0] };

    my @xy_plots;
    for (@possible_xy_plots) {
        my %para_set = %{$_};
        $para_set{xy} =~ /\$(.*),.*\$(.*)/ or die $para_set{xy};

        if ((grep $1 eq $_, @prop_types) && (grep $2 eq $_, @prop_types)) {
            push @xy_plots, \%para_set;
        }
    }

    mkpath($plot_dir);
    mkpath(catdir($plot_dir, 'png'));

    my $png_convert_calls;

    #Read in data
    push @r_code, "adhesions = read.table('$data_dir/$cfg{individual_adhesions_props_file}',header=T,sep=',');\n";

    #Build the plots
    foreach (@xy_plots) {
        my %parameters = %{$_};
        my $output_file = catfile($plot_dir, $parameters{file_name});

        my $png_file = $parameters{file_name};
        $png_file =~ s/\.pdf/\.png/;
        my $output_file_png = catfile($plot_dir, 'png', $png_file);

        push @r_code, "pdf('$output_file',$parameters{pdf_para})\n";
        push @r_code, "par(mar=c(4,4,0.5,0.5),bty='n')\n";
        push @r_code, "plot($parameters{xy},xlab=$parameters{xlab},ylab=$parameters{ylab},$parameters{plot_para})\n";
        push @r_code, "dev.off();\n";
        $png_convert_calls .= "convert $output_file $output_file_png\n";
    }
    &Math::R::execute_commands(\@r_code);
    system($png_convert_calls);
}

#######################################
#Adhesion Lineage Property Collection
#######################################
sub gather_ad_lineage_properties {
    my %props;

    $props{longevity}                     = &gather_longevities;
    $props{merge_count}                   = &gather_merge_count;
    $props{death_status}                  = &gather_death_status;
    $props{Average_adhesion_signal}       = &gather_prop_seq("Average_adhesion_signal");
    $props{Background_corrected_signal}       = &gather_prop_seq("Background_corrected_signal");
    $props{Eccentricity}                  = &gather_prop_seq("Eccentricity");
    $props{Solidity}                      = &gather_prop_seq("Solidity");
    $props{ad_sig}                        = &gather_average_ad_sig($props{Average_adhesion_signal});
    $props{Max_adhesion_signal}           = &gather_prop_seq("Max_adhesion_signal");
    ($props{All_speeds}, $props{velocity}) = &gather_adhesion_speeds;

    ($props{average_speeds}, $props{variance_speeds}, $props{max_speeds}) = &gather_speed_props($props{All_speeds});

    if (grep "Area" eq $_, @available_data_types) {
        $props{Area}         = &gather_prop_seq("Area");
        $props{largest_area} = &gather_largest_entry($props{Area});
    }

    if (grep "Centroid_dist_from_center" eq $_, @available_data_types) {
        $props{Centroid_dist_from_center} = &gather_prop_seq("Centroid_dist_from_center");
        $props{starting_center_dist}      = &gather_first_entry($props{Centroid_dist_from_center});
        $props{ending_center_dist}        = &gather_last_entry($props{Centroid_dist_from_center});
    }

    if (grep "Centroid_dist_from_edge" eq $_, @available_data_types) {
        $props{Centroid_dist_from_edge} = &gather_prop_seq("Centroid_dist_from_edge");
        $props{starting_edge_dist}      = &gather_first_entry($props{Centroid_dist_from_edge});
        $props{ending_edge_dist}        = &gather_last_entry($props{Centroid_dist_from_edge});
    }

    return %props;
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
    return \@longevities;
}

sub gather_prop_seq {
    my $prop        = $_[0];
    my $default_val = "NaN";
    $default_val = $_[1] if (scalar(@_) > 1);

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
                die "Can't find value for property ($prop) in image number ($i_num), ad number ($ad_num)\n" . 
                    "Data keys in this I number:\n" . join(" ",keys %{$data_sets{$i_num}});
            }

            push @{ $prop_vals[$i] }, ${ $data_sets{$i_num}{$prop} }[$ad_num];
        }
    }
    return \@prop_vals;
}

sub gather_largest_entry {
    my @data = @{ $_[0] };

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

sub gather_first_entry {
    my @data = @{ $_[0] };

    my @starting_data;
    for my $i (0 .. $#data) {
        my $first_data_index = (grep $data[$i][$_] ne "NaN", (0 .. $#{ $data[$i] }))[0];
        $starting_data[$i] = $data[$i][$first_data_index];
    }
    return \@starting_data;
}

sub gather_last_entry {
    my @data = @{ $_[0] };

    my @last_data;
    for my $i (0 .. $#data) {
        my $last_data_index = (grep $data[$i][$_] ne "NaN", (0 .. $#{ $data[$i] }))[-1];
        $last_data[$i] = $data[$i][$last_data_index];
    }
    return \@last_data;
}

sub gather_average_ad_sig {
    my @signal = @{ $_[0] };

    my @pax_sig;
    for my $i (0 .. $#signal) {
        my $stat = Statistics::Descriptive::Full->new();
        for my $j (0 .. $#{ $signal[$i] }) {
            $stat->add_data($signal[$i][$j]) if ($signal[$i][$j] ne "NaN");
        }

        push @pax_sig, $stat->mean();
    }

    return \@pax_sig;
}

sub gather_adhesion_speeds {
    my @speed;
    my @velocity;
    my @data_keys = sort keys %data_sets;
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

    #print "$cfg{exp_name}:", $stat->mean()*60,"\n";
    return \@av_speeds, \@var_speeds, \@max_speeds;
}

sub gather_merge_count {
    my @merge_count = map 0, (0 .. $#tracking_mat);
    my @data_keys = sort keys %data_sets;
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
    my @data_keys = sort keys %data_sets;
    for my $i (0 .. $#tracking_mat) {
        my $in_seq = 0;
        for my $j (0 .. $#{ $tracking_mat[$i] }) {
            $in_seq = 1 if ($tracking_mat[$i][$j] >= 0);
            $death_status[$i] = 1 if ($tracking_mat[$i][$j] == -1 && $in_seq);
            $in_seq = 0 if ($tracking_mat[$i][$j] <= -1 && $in_seq);

        }
    }

    return \@death_status;
}

sub output_adhesion_lineage_props {
    mkpath(catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}));
    mkpath(catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{lineage_ts_folder}));

    my @lin_summary_data = &gather_lineage_summary_data;
    my $output_file = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{lineage_summary_props_file});
    &output_mat_csv(\@lin_summary_data, $output_file);

    my @ts_props = qw(All_speeds Area Centroid_dist_from_edge
      Centroid_dist_from_center Average_adhesion_signal 
      Max_adhesion_signal Eccentricity Solidity Background_corrected_signal);
    foreach (@ts_props) {
        my $this_result = $_;
        next if (not(grep $this_result eq $_, keys %ad_lineage_props));

        my $output_file =
          catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{lineage_ts_folder},
            $this_result . ".csv");
        &output_mat_csv($ad_lineage_props{$_}, $output_file);
    }
}

sub gather_lineage_summary_data {
    my @possible_props = qw(longevity largest_area starting_edge_dist ending_edge_dist
      starting_center_dist ending_center_dist merge_count death_status average_speeds max_speeds ad_sig);

    my @lin_summary_data;
    for (@possible_props) {
        my $this_prop = $_;
        next if (not(grep $this_prop eq $_, keys %ad_lineage_props));

        push @{ $lin_summary_data[0] }, $this_prop;

        for my $i (0 .. $#{ $ad_lineage_props{$this_prop} }) {
            push @{ $lin_summary_data[ $i + 1 ] }, $ad_lineage_props{$this_prop}[$i];
        }
    }

    return @lin_summary_data;
}

sub build_lineage_plots {
    my @r_code;

    my $data_dir = catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder});
    my $plot_dir = catdir($data_dir, $cfg{plot_folder});

    my $xy_default = "pch=19,cex=0.4";

    my $pdf_default       = "height=12, width=12, pointsize=24";
    my @possible_xy_plots = (
        {
            xy        => "lineages\$longevity,lineages\$starting_edge_dist",
            xlab      => "'Longevity (min)'",
            ylab      => "expression(paste('Starting Distance from Edge (', mu, 'm)'))",
            file_name => "longev_vs_s_dist.pdf",
            plot_para => $xy_default,
            pdf_para  => $pdf_default,
        },
        {
            xy        => "lineages\$longevity,lineages\$largest_area",
            xlab      => "'Longevity (min)'",
            ylab      => "expression(paste('Largest Area Attained (', mu, m^2, ')'))",
            file_name => "longev_vs_largest_area.pdf",
            plot_para => $xy_default,
            pdf_para  => $pdf_default,
        },
        {
            xy        => "lineages\$longevity,lineages\$ad_sig",
            xlab      => "'Longevity (min)'",
            ylab      => "'Paxillin Concentration (AU)'",
            file_name => "longev_vs_pax.pdf",
            plot_para => $xy_default,
            pdf_para  => $pdf_default,
        },
    );

    my @summary_data = &gather_lineage_summary_data;
    my @prop_types   = @{ $summary_data[0] };

    my @xy_plots;
    for (@possible_xy_plots) {
        my %para_set = %{$_};
        $para_set{xy} =~ /\$(.*),.*\$(.*)/ or die $para_set{xy};

        if ((grep $1 eq $_, @prop_types) && (grep $2 eq $_, @prop_types)) {
            push @xy_plots, \%para_set;
        }
    }

    mkpath($plot_dir);
    mkpath(catdir($plot_dir, 'png'));

    my $png_convert_calls;

    #Read in data
    push @r_code, "lineages = read.table('$data_dir/$cfg{lineage_summary_props_file}',header=T,sep=',');\n";

    #Build the plots
    foreach (@xy_plots) {
        my %parameters = %{$_};
        my $output_file = catfile($plot_dir, $parameters{file_name});

        my $png_file = $parameters{file_name};
        $png_file =~ s/\.pdf/\.png/;
        my $output_file_png = catfile($plot_dir, 'png', $png_file);

        push @r_code, "pdf('$output_file',$parameters{pdf_para})\n";
        push @r_code, "par(mar=c(4,4,0.5,0.5),bty='n')\n";
        push @r_code, "plot($parameters{xy},xlab=$parameters{xlab},ylab=$parameters{ylab},$parameters{plot_para})\n";
        push @r_code, "dev.off();\n";
        $png_convert_calls .= "convert $output_file $output_file_png\n";
    }
    &Math::R::execute_commands(\@r_code);
    system($png_convert_calls);
}

#######################################
# Linear Region Finding Code
#######################################
sub run_R_linear_region_code {
    my $data_dir = catdir($cfg{exp_results_folder},$cfg{adhesion_props_folder},$cfg{lineage_ts_folder});
    my $R_cmd = "R CMD BATCH --vanilla --args -$data_dir linear_regions.R";
    system($R_cmd);
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
