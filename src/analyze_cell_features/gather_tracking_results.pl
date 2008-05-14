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
use Imager;
use Statistics::Descriptive;
use Statistics::Distributions;

#local libraries
use lib "../lib";
use Config::Adhesions;
use Image::Data::Collection;
use Text::CSV::Simple::Extra;
use Math::R;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|config=s", "debug|d", "input|i=s", "output|o=s", "skip_pix_props", "tracking_mat=s");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg     = $ad_conf->get_cfg_hash;

if (exists $opt{tracking_mat}) {
    $cfg{tracking_output_file} = $opt{tracking_mat};
}
if (exists $opt{tracking_mat}) {
    $cfg{adhesion_props_folder} = catdir($cfg{adhesion_props_folder}, dirname($opt{tracking_mat}));
}

###############################################################################
# Main Program
###############################################################################
print "Gathering/Converting Data Files\n" if $opt{debug};

my @data_files;
push @data_files, split(/\s+/, $cfg{general_data_files});
push @data_files, split(/\s+/, $cfg{lineage_analysis_data_files});
my %data_sets = &Image::Data::Collection::gather_data_sets(\%cfg, \%opt, \@data_files);
%data_sets = &Image::Data::Collection::trim_data_sets(\%cfg, \%opt, \%data_sets);
%data_sets = &convert_data_to_units(\%data_sets, \%cfg);
my @available_data_types = &gather_data_types;

print "\n\nCollecting Tracking Matrix\n" if $opt{debug};
my @tracking_mat = &Image::Data::Collection::read_in_tracking_mat(\%cfg, \%opt);

print "\n\nCreating Pixel Properties Plots\n" if $opt{debug};
my @image_values;
if (not($opt{skip_pix_props})) {
    @image_values = &gather_image_value_props(\%cfg, \%opt);
    &output_image_props;
}
&build_image_props_plots;

print "\n\nCreating Individual Adhesion Properties Plots\n" if $opt{debug};
my @single_ad_props = &gather_single_ad_props(\%cfg, \%opt);
&output_single_adhesion_props;
&build_single_ad_plots;

print "\n\nCreating Adhesion Lineage Property Plots\n", if $opt{debug};
my %ad_lineage_props = &gather_ad_lineage_properties;
&output_adhesion_lineage_props;
&build_lineage_plots;

print "\n\nGathering Adhesion Property Sequences\n", if $opt{debug};
my %ad_lineage_prop_seqs = &gather_property_sequences(\@tracking_mat, \%data_sets);
&output_adhesion_prop_seqs;

#print "\n\nCreating Vis Project Files\n", if $opt{debug};
#create_vis_project_files;

if (exists $cfg{treatment_time}) {
    print "\n\nGathering Treatment Results\n", if $opt{debug};
    my %treatment_results = &gather_treatment_results;
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
            } elsif (   (grep $data_type eq $_, qw(Class))
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

######################################
#Pixel Value Props
######################################
sub gather_image_value_props {
    my %cfg = %{ $_[0] };
    my %opt = %{ $_[1] };

    my @focal_image_files = sort <$cfg{individual_results_folder}/*/$cfg{adhesion_image_file}>;
    my @cell_mask_files   = sort <$cfg{individual_results_folder}/*/cell_mask.png>;

    my @overall_stats = [qw(ImageNum CellSize AdhesionTotalSize Minimum Maximum Average)];

    print "Working on image #: " if $opt{debug};

    foreach my $i (0 .. $#focal_image_files) {
        my $image_num;
        if ($focal_image_files[$i] =~ /$cfg{individual_results_folder}\/(\d+)\/$cfg{adhesion_image_file}/) {
            $image_num = $1;
            next if grep $1 == $_, @{ $cfg{exclude_image_nums} };
        }
        print "$image_num " if $opt{debug};

        my $focal_img = Imager->new;
        $focal_img->read(file => "$focal_image_files[$i]") or die;

        my $cell_mask_img = Imager->new;
        if (scalar(@cell_mask_files) > 0) {
            $cell_mask_img->read(file => "$cell_mask_files[$i]") or die;
        }

        my @vals;
        my @intra_cellular;
        for my $y (0 .. $focal_img->getheight - 1) {
            my @gray = $focal_img->getsamples(y => $y, type => 'float');
            push @vals, @gray;

            my @in_cell_indexes;
            if (scalar(@cell_mask_files) > 0) {
                my @cell_mask = $cell_mask_img->getsamples(y => $y);
                my @in_cell_indexes = grep $cell_mask[$_], (0 .. $#cell_mask);
                push @intra_cellular, @gray[@in_cell_indexes];
            } else {
                push @intra_cellular, @gray;
            }

        }

        my $total_stat = Statistics::Descriptive::Full->new();
        $total_stat->add_data(@vals);

        my $in_cell_stat = Statistics::Descriptive::Full->new();
        $in_cell_stat->add_data(@intra_cellular);

        my $all_ad_size = 0;
        $all_ad_size += $_ foreach (@{ $data_sets{$image_num}{Area} });

        push @overall_stats,
          [
            $image_num,       $data_sets{$image_num}{"Cell_size"}[0], $all_ad_size,
            $total_stat->min, $total_stat->max,                       $in_cell_stat->mean
          ];
    }

    return @overall_stats;
}

sub output_image_props {
    if (not(-e catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}))) {
        mkpath(catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}));
    }

    my $output_file = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{image_props_file});
    &output_mat_csv(\@image_values, $output_file);
}

sub build_image_props_plots {
    my $pdf_default = "width=12, height=12, pointsize=24";

    my @plot_vars = (
        {
            pdf_para  => $pdf_default,
            xy        => "image_values\$ImageNum,image_values\$Maximum",
            file_name => "pix_max.pdf",
            xlab      => "\"Image Number\"",
            ylab      => "\"Maximum Normalized Fluorescence (AU)\"",
            plot_opt  => "type=\"l\",ylim=c(0,1)",
        },
        {
            pdf_para  => $pdf_default,
            xy        => "image_values\$ImageNum,image_values\$Average",
            file_name => "pix_average.pdf",
            xlab      => "\"Image Number\"",
            ylab      => "\"Average Normalized Fluorescence (AU)\"",
            plot_opt  => "type=\"l\"",
        },
    );

    my $png_convert_calls;

    my $data_dir = catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder});
    my $plot_dir = catdir($data_dir, $cfg{plot_folder});

    mkpath($plot_dir);
    mkpath(catdir($plot_dir, 'png'));

    my @r_code;

    foreach (@plot_vars) {
        my %para = %{$_};
        my $output_file = catfile($plot_dir, $para{file_name});

        my $png_file = $para{file_name};
        $png_file =~ s/\.pdf/\.png/;
        my $output_file_png = catfile($plot_dir, 'png', $png_file);

        #Read in data
        push @r_code, "image_values = read.table('$data_dir/$cfg{image_props_file}',header=T,sep=',');\n";

        #Build the plots
        push @r_code, "pdf('$output_file',$para{pdf_para})\n";
        push @r_code, "par(mar=c(4,4,0.5,0.5),bty='n')\n";
        push @r_code, "plot($para{xy},xlab=$para{xlab},ylab=$para{ylab},$para{plot_opt})\n";
        push @r_code, "dev.off();\n";

        $png_convert_calls .= "convert $output_file $output_file_png\n";
    }
    &Math::R::execute_commands(\@r_code);
    system($png_convert_calls);
}

#######################################
#Single Adhesion Props
#######################################
sub gather_single_ad_props {
    my @data;

    my @possible_data_types =
      qw(Area Average_adhesion_signal Centroid_dist_from_edge Centroid_dist_from_center Variance_adhesion_signal);

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
    if (not(-e catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}))) {
        mkpath(catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}));
    }

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

    $props{longevity}               = &gather_longevities;
    $props{merge_count}             = &gather_merge_count;
    $props{Average_adhesion_signal} = &gather_prop_seq("Average_adhesion_signal");
    $props{ad_sig}                  = &gather_average_ad_sig($props{Average_adhesion_signal});
    $props{Average_pixel_adhesion_signal} = &gather_prop_seq("Average_pixel_adhesion_signal");
    $props{Max_adhesion_signal} = &gather_prop_seq("Max_adhesion_signal");
    ($props{All_speeds},$props{velocity})              = &gather_adhesion_speeds;

    ($props{average_speeds}, $props{variance_speeds}, $props{max_speeds}) = &gather_speed_props($props{All_speeds});

    if (grep "Area" eq $_, @available_data_types) {
        $props{Area}         = &gather_prop_seq("Area");
        $props{largest_area} = &gather_largest_areas($props{Area});
    }

    if (grep "Centroid_dist_from_center" eq $_, @available_data_types) {
        $props{Centroid_dist_from_center} = &gather_prop_seq("Centroid_dist_from_center");
        $props{starting_center_dist}      = &gather_starting_dist($props{Centroid_dist_from_center});
    }

    if (grep "Centroid_dist_from_edge" eq $_, @available_data_types) {
        $props{Centroid_dist_from_edge} = &gather_prop_seq("Centroid_dist_from_edge");
        $props{starting_edge_dist}      = &gather_starting_dist($props{Centroid_dist_from_edge});
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
                print "$i_num, $ad_num, $prop";
                die;
            }

            push @{ $prop_vals[$i] }, ${ $data_sets{$i_num}{$prop} }[$ad_num];
        }
    }
    return \@prop_vals;
}

sub gather_largest_areas {
    my @areas = @{ $_[0] };

    my @largest_areas;
    for my $i (0 .. $#areas) {
        my $largest = 0;
        for my $j (0 .. $#{ $areas[$i] }) {
            next if ($areas[$i][$j] eq "NaN");
            $largest = $areas[$i][$j] if ($largest < $areas[$i][$j]);
        }
        push @largest_areas, $largest;
    }
    return \@largest_areas;
}

sub gather_starting_dist {
    my @dists = @{ $_[0] };

    my @starting_dists;
    for my $i (0 .. $#dists) {
        for my $j (0 .. $#{ $dists[$i] }) {
            if ($dists[$i][$j] ne "NaN") {
                $starting_dists[$i] = $dists[$i][$j];
                last;
            }
        }
    }
    return \@starting_dists;
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
                push @{ $velocity[$i] }, [($start_x - $end_x,$start_y - $end_y)];
            } else {
                push @{ $speed[$i] }, "NaN";
                push @{ $velocity[$i] }, "NaN";
            }
        }
    }

    return (\@speed,\@velocity);
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

sub output_adhesion_lineage_props {
    if (not(-e catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}))) {
        mkpath(catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}));
    }
    if (not(-e catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{lineage_ts_folder}))) {
        mkpath(catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{lineage_ts_folder}));
    }

    my @lin_summary_data = &gather_lineage_summary_data;
    my $output_file = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{lineage_summary_props_file});
    &output_mat_csv(\@lin_summary_data, $output_file);

    my @ts_props = qw(All_speeds Area Centroid_dist_from_edge
      Centroid_dist_from_center Average_adhesion_signal Average_pixel_adhesion_signal 
      Max_adhesion_signal);
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
    my @possible_props = qw(longevity largest_area starting_edge_dist
      starting_center_dist starting_center_dist merge_count average_speeds
      max_speeds ad_sig);

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
# Vis Project
######################################

sub create_vis_project_files {
    my $base = catdir($cfg{results_folder}, 'vis_project', $cfg{exp_name});
    mkpath($base);
    
    my @local_tracking_mat;
    for (2) {
        my $required_longevity = $_;
        for my $i (0 .. $#{ $ad_lineage_props{'longevity'} }) {
            my $this_longev = $ad_lineage_props{'longevity'}[$i];
            if (   $this_longev >= $required_longevity
                || $this_longev >= $#tracking_mat) {
                push @local_tracking_mat, $tracking_mat[$i];
            }
        }
    }
    
    my @data_keys = sort keys %data_sets;
    my %i_num_to_col = map {$data_keys[$_] => $_} (0 .. $#data_keys);
    
    #Remove the last data key since it lacks several of the properties
    pop @data_keys;

    for my $i_num (@data_keys) {
        my @header = ("# vtk DataFile Version 1.0\n",
                $cfg{exp_name} . ", Image Num: $i_num Data Set\n",
                "ASCII\n\n",
                );

        my $col_num = $i_num_to_col{$i_num};
            
        my @ad_nums = map {
            if($local_tracking_mat[$_][$col_num] >= 0) {
                $local_tracking_mat[$_][$col_num];
            } else {
                ();
            }
        } (0 .. $#local_tracking_mat);
        my $point_count = scalar(@ad_nums);
        
        my @grid_data = &gather_grid_data($i_num, \@ad_nums, $point_count);
        my @cell_data = &gather_cell_data($i_num, \@ad_nums, $point_count);
        my @cell_type_data = &gather_cell_type_data($i_num, \@ad_nums, $point_count);
        my @size_data = &gather_size_data($i_num, \@ad_nums, $point_count);
        my @pax_data = &gather_pax_data($i_num, \@ad_nums, $point_count);
        my @velo_data = &gather_velo_data($i_num, \@ad_nums, $point_count, $col_num);
        my @total_lifetime_data = &gather_total_lifetime_data($i_num, \@ad_nums, \@local_tracking_mat, $col_num);
        my @lifetime_reminaing_data = &gather_lifetime_remaining_data($i_num, \@ad_nums, \@local_tracking_mat, $col_num);
        my @min_size_data = &gather_min_size_data($i_num, \@ad_nums, \@local_tracking_mat, $col_num, $point_count);
        my @max_size_data = &gather_lifetime_remaining_data($i_num, \@ad_nums, \@local_tracking_mat, $col_num);
        

        mkpath(catdir($base, 'G1'));
        mkpath(catdir($base, 'G2'));
        mkpath(catdir($base, 'G3'));
        open OUTPUT, '>' . catfile($base, 'G1', $i_num .'.vtk');
        print OUTPUT @header, @grid_data, @cell_data, @cell_type_data, 
                     @size_data, @lifetime_reminaing_data, @velo_data;
        close OUTPUT;
        
        open OUTPUT, '>' . catfile($base,'G2', $i_num .'.vtk');
        print OUTPUT @header, @grid_data, @cell_data, @cell_type_data, 
                     @size_data, @pax_data, @velo_data;
        close OUTPUT;
        
        open OUTPUT, '>' . catfile($base,'G3', $i_num .'_min.vtk');
        print OUTPUT @header, @grid_data, @cell_data, @cell_type_data, @min_size_data;
        close OUTPUT;
        
        open OUTPUT, '>' . catfile($base,'G3', $i_num .'_max.vtk');
        print OUTPUT @header, @grid_data, @cell_data, @cell_type_data, @max_size_data;
        close OUTPUT;
        
        open OUTPUT, '>' . catfile($base,'G3', $i_num .'_cur.vtk');
        print OUTPUT @header, @grid_data, @cell_data, @cell_type_data, @size_data;
        close OUTPUT;
     }
}

sub gather_grid_data {
    my $i_num = $_[0];
    my @ad_nums = @{$_[1]};
    my $point_count = $_[2];
    
    my @lines = ("DATASET UNSTRUCTURED_GRID\n",
                 "POINTS ". $point_count. " float\n",
                ); 
    
    for (@ad_nums) {
        #die "$i_num, $_" if not(defined($data_sets{$i_num}{Centroid_x}[$_]));
        push @lines, sprintf("%.1f",$data_sets{$i_num}{Centroid_x}[$_]). " " . sprintf("%.1f",$data_sets{$i_num}{Centroid_y}[$_]) . " 0.0\n";
    }
    push @lines, "\n";

    return @lines;
}

sub gather_cell_data {
    my $i_num = $_[0];
    my @ad_nums = @{$_[1]};
    my $point_count = $_[2];
    
    my @lines = ("CELLS $point_count 1\n"); 
    
    for (0 .. $#ad_nums) {
        my $num = $_ + 1;
        #die "$i_num, $_" if not(defined($data_sets{$i_num}{Centroid_x}[$_]));
        push @lines, "1 " . $num . "\n";
    }
    push @lines, "\n";

    return @lines;
}

sub gather_cell_type_data {
    my $i_num = $_[0];
    my @ad_nums = @{$_[1]};
    my $point_count = $_[2];
    
    my @lines = ("CELL_TYPES $point_count\n"); 
    
    for (0 .. $#ad_nums) {
        #die "$i_num, $_" if not(defined($data_sets{$i_num}{Centroid_x}[$_]));
        push @lines, "1\n";
    }
    push @lines, "\n";

    return @lines;
}

sub gather_size_data {
    my $i_num = $_[0];
    my @ad_nums = @{$_[1]};
    my $point_count = $_[2];
    
    my @lines = ("POINT_DATA $point_count\n",
                 "SCALARS size float\n",
                 "LOOKUP_TABLE default\n",
                 ); 
    
    for (@ad_nums) {
        #die "$i_num, $_" if not(defined($data_sets{$i_num}{Centroid_x}[$_]));
        push @lines, $data_sets{$i_num}{Area}[$_] . "\n";
    }
    push @lines, "\n";

    return @lines;
}

sub gather_min_max_size_mat {
    my $i_num = $_[0];
    my @ad_nums = @{$_[1]};
    my @local_tracking_mat = @{$_[2]};
    my $col_num = $_[3];
    
    my @data_keys = sort keys %data_sets;
    my %i_num_to_col = map {$data_keys[$_] => $_} (0 .. $#data_keys);
    
    #Remove the last data key since it lacks several of the properties
    pop @data_keys;
    
    my @min_max; 
        
    for my $ad_num (@ad_nums) {
        my $row_num = (grep $local_tracking_mat[$_][$col_num] == $ad_num, (0 .. $#local_tracking_mat))[0];
        
        my $cur_size = $data_sets{$i_num}{Area}[$ad_num];

        my $increasing = 1;
        my $offset = 1;
        my @this_mm = ($cur_size);
        while ($increasing) {
            my $bottom_offset = $offset;
            my $top_offset = $offset;

            $top_offset = $#{$local_tracking_mat[$row_num]} - $col_num - 1 if ($col_num + $offset >= $#{$local_tracking_mat[$row_num]} - 1);
            $bottom_offset = 0 if ($col_num - $offset >= 0);

            last if ($local_tracking_mat[$row_num][$bottom_offset] < 0);
            next if ($local_tracking_mat[$row_num][$bottom_offset] < 0);

            last if ($local_tracking_mat[$row_num][$top_offset] < 0);
            next if ($local_tracking_mat[$row_num][$top_offset] < 0);

            my $early_ad = $local_tracking_mat[$row_num][$col_num - $bottom_offset];
            my $early_size =  $data_sets{$data_keys[$col_num - $bottom_offset]}{Area}[$early_ad];
            if (not($early_size < $this_mm[0])) {
                $increasing = 0;
            } else {
                unshift(@this_mm, $early_size);
            }

            my $late_ad = $local_tracking_mat[$row_num][$col_num + $top_offset];
            if (not(defined($data_sets{$data_keys[$col_num + $top_offset]}{Area}[$late_ad]))) {
                die $data_keys[$col_num + $offset];
            }
            my $late_size = $data_sets{$data_keys[$col_num + $top_offset]}{Area}[$late_ad];
            if (not($late_size > $this_mm[-1])) {
                $increasing = 0;
            } else {
                push(@this_mm, $late_size);
            }

            $offset++;
        }

        if (scalar(@this_mm) < 3) {
            my $decreasing = 1;
            $offset = 1;
            @this_mm = ($cur_size);
            while ($decreasing) {
                my $bottom_offset = $offset;
                my $top_offset = $offset;

                $top_offset = $#{$local_tracking_mat[$row_num]} - $col_num - 1 if ($col_num + $offset >= $#{$local_tracking_mat[$row_num]} - 1);
                $bottom_offset = 0 if ($col_num - $offset >= 0);

                last if ($local_tracking_mat[$row_num][$bottom_offset] < 0);
                next if ($local_tracking_mat[$row_num][$bottom_offset] < 0);

                last if ($local_tracking_mat[$row_num][$top_offset] < 0);
                next if ($local_tracking_mat[$row_num][$top_offset] < 0);

                my $early_ad = $local_tracking_mat[$row_num][$col_num - $bottom_offset];
                my $early_size =  $data_sets{$data_keys[$col_num - $bottom_offset]}{Area}[$early_ad];
                if (not($early_size > $this_mm[0])) {
                    $decreasing = 0;
                } else {
                    unshift(@this_mm, $early_size);
                }

                my $late_ad = $local_tracking_mat[$row_num][$col_num + $top_offset];
                if (not(defined($data_sets{$data_keys[$col_num + $top_offset]}{Area}[$late_ad]))) {
                    die $data_keys[$col_num + $offset];
                }
                my $late_size = $data_sets{$data_keys[$col_num + $top_offset]}{Area}[$late_ad];
                if (not($late_size < $this_mm[-1])) {
                    $decreasing = 0;
                } else {
                    push(@this_mm, $late_size);
                }

                $offset++;
            }
        }
        
        if (scalar(@this_mm) > 3) {
            if ($this_mm[0] > $this_mm[-1]) {
                push @min_max, [$this_mm[-1],$this_mm[0]];
            } elsif ($this_mm[-1] > $this_mm[0]) {
                push @min_max, [$this_mm[0],$this_mm[-1]];
            } else {
                push @min_max, [0,0];
            }
        } else {
            push @min_max, [0,0];
        }
    }

    return @min_max;
}

sub gather_min_size_data {
    my $i_num = $_[0];
    my @ad_nums = @{$_[1]};
    my @local_tracking_mat = @{$_[2]};
    my $col_num = $_[3];
    my $point_count = $_[4];

    my @lines = ("POINT_DATA $point_count\n",
                 "SCALARS size float\n",
                 "LOOKUP_TABLE default\n",
                 ); 

    my @min_max = &gather_min_max_size_mat(@_);
    
    for my $ref (@min_max) {
        push @lines, ${$ref}[0] . "\n";
    }
    
    return @lines;
}

sub gather_max_size_data {
    my $i_num = $_[0];
    my @ad_nums = @{$_[1]};
    my @local_tracking_mat = @{$_[2]};
    my $col_num = $_[3];
    my $point_count = $_[4];

    my @lines = ("POINT_DATA $point_count\n",
                 "SCALARS size float\n",
                 "LOOKUP_TABLE default\n",
                 ); 

    my @min_max = &gather_min_max_size_mat(@_);
    
    for my $ref (@min_max) {
        push @lines, ${$ref}[1] . "\n";
    }
    
    return @lines;
}

sub gather_pax_data {
    my $i_num = $_[0];
    my @ad_nums = @{$_[1]};
    my $point_count = $_[2];
    
    my @lines = ("SCALARS pax float\n",
                 "LOOKUP_TABLE default\n",
                 ); 
    
    for (@ad_nums) {
        #die "$i_num, $_" if not(defined($data_sets{$i_num}{Centroid_x}[$_]));
        push @lines, $data_sets{$i_num}{Average_adhesion_signal}[$_] . "\n";
    }
    push @lines, "\n";

    return @lines;
}

sub gather_total_lifetime_data {
    my $i_num = $_[0];
    my @ad_nums = @{$_[1]};
    my @local_tracking_mat = @{$_[2]};
    my $col_num = $_[3];
    
    my @lines = ("SCALARS totallife float\n",
                 "LOOKUP_TABLE default\n",
                 ); 

    for my $ad_num (@ad_nums) {
        my @row_num = grep $local_tracking_mat[$_][$col_num] == $ad_num, (0 .. $#local_tracking_mat);
        die "Too many row matches in total lifetime data gathering" if (scalar(@row_num) > 1);
        push @lines, $ad_lineage_props{longevity}[$row_num[0]] . "\n";
    }
    push @lines, "\n";

    return @lines;
}

sub gather_lifetime_remaining_data {
    my $i_num = $_[0];
    my @ad_nums = @{$_[1]};
    my @local_tracking_mat = @{$_[2]};
    my $col_num = $_[3];
    
    my @lines = ("SCALARS liferemaining float\n",
                 "LOOKUP_TABLE default\n",
                 ); 

    for my $ad_num (@ad_nums) {
        my @row_num = grep $local_tracking_mat[$_][$col_num] == $ad_num, (0 .. $#local_tracking_mat);
        die "Too many row matches in lifetime remaining data gathering" if (scalar(@row_num) > 1);
        
        my @time_points_left = grep $local_tracking_mat[$row_num[0]][$_] >= 0, ($col_num + 1 .. $#{$local_tracking_mat[$row_num[0]]}); 
        
        push @lines, scalar(@time_points_left) . "\n";
    }
    push @lines, "\n";

    return @lines;
}

sub gather_velo_data {
    my $i_num = $_[0];
    my @ad_nums = @{$_[1]};
    my $point_count = $_[2];
    my $col_num = $_[3];
    
    my @lines = ("VECTORS velocity float\n",
                 ); 

    for (@ad_nums) {
        if (ref($ad_lineage_props{velocity}[$_][$col_num]) eq "ARRAY") {
            my @velo = map sprintf('%.5f',$_), @{$ad_lineage_props{velocity}[$_][$col_num]};
            push @lines, join(" ",(@velo, 0.0)) . "\n";
        } else {
            push @lines, "0.0 0.0 0.0\n";
        }
    }
    push @lines, "\n";

    return @lines;
}

#######################################
# Adhesion Lineage Property Sequence
######################################
sub gather_property_sequences {
    my @tracking_mat = @{ $_[0] };
    my %data_sets    = %{ $_[1] };

    my %seqs;
    $seqs{area}{increasing} = &gather_patterned_property_seqs("Area", \@tracking_mat, \%data_sets, \&greater_or_equal);
    $seqs{area}{decreasing} = &gather_patterned_property_seqs("Area", \@tracking_mat, \%data_sets, \&less_or_equal);

    $seqs{ad_sig}{increasing} =
      &gather_patterned_property_seqs("Average_adhesion_signal", \@tracking_mat, \%data_sets, \&greater_or_equal);
    $seqs{ad_sig}{decreasing} =
      &gather_patterned_property_seqs("Average_adhesion_signal", \@tracking_mat, \%data_sets, \&less_or_equal);
    return %seqs;
}

sub gather_patterned_property_seqs {
    my $property     = $_[0];
    my @tracking_mat = @{ $_[1] };
    my %data_sets    = %{ $_[2] };
    my $func_ref     = $_[3];

    my @image_nums = sort keys %data_sets;

    my @image_num_seqs;

    for my $i (0 .. $#tracking_mat) {
        my @data_seq;
        my @i_num_seq;
        for my $j (0 .. $#{ $tracking_mat[$i] }) {
            next if ($tracking_mat[$i][$j] <= -1);
            push @i_num_seq, $image_nums[$j];
            push @data_seq,  $data_sets{ $image_nums[$j] }{$property}[ $tracking_mat[$i][$j] ];
        }

        my @indexes = &gather_indexes_seqs_by_func(\@data_seq, $func_ref);
        for my $j (0 .. $#indexes) {
            push @{ $image_num_seqs[$i] }, [ @i_num_seq[ @{ $indexes[$j] } ] ];
        }
        if (not(@indexes)) {
            push @{ $image_num_seqs[$i] }, ();
        }
    }
    return \@image_num_seqs;
}

sub gather_indexes_seqs_by_func {
    my @set           = @{ $_[0] };
    my $deciding_func = $_[1];

    my @deciding_func_true = map { $deciding_func->($set[$_], $set[ $_ - 1 ]) } (1 .. $#set);
    unshift @deciding_func_true, 0;

    my @index_seqs;
    @{ $index_seqs[0] } = ();

    my $cur_index = 0;
    while ($cur_index <= $#deciding_func_true) {
        while ($cur_index <= $#deciding_func_true && $deciding_func_true[$cur_index]) {
            push @{ $index_seqs[$#index_seqs] }, $cur_index;
            $cur_index++;
        }
        while ($cur_index <= $#deciding_func_true && !$deciding_func_true[$cur_index]) {
            if ($#{ $index_seqs[$#index_seqs] } > -1) {
                push @index_seqs, [];
            }
            $cur_index++;
        }
    }

    for (0 .. $#index_seqs) {
        my @this_seq = @{ $index_seqs[$_] };
        if ($#this_seq > -1) {
            unshift @this_seq, $this_seq[0] - 1;
            @{ $index_seqs[$_] } = @this_seq;
        } else {
            pop @index_seqs;
        }
    }

    return @index_seqs;
}

sub greater_or_equal {
    if ($_[0] >= $_[1]) {
        return 1;
    } else {
        return 0;
    }
}

sub less_or_equal {
    if ($_[0] <= $_[1]) {
        return 1;
    } else {
        return 0;
    }
}

sub output_adhesion_prop_seqs {
    if (not(-e catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}))) {
        mkpath(catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder}));
    }

    &output_sequence_trimmed_mat(\@{ $ad_lineage_prop_seqs{area}{increasing} }, "_increasing_area");
}

sub output_sequence_trimmed_mat {
    my @i_num_lists = @{ $_[0] };
    my $file_suffix = $_[1];

    my @i_nums = sort keys %data_sets;
    my %i_num_to_index = map { $i_nums[$_] => $_ } (0 .. $#i_nums);

    my @trimmed_tracking_mat;

    for my $i (0 .. $#i_num_lists) {
        for my $j (0 .. $#{ $i_num_lists[$i] }) {
            my @indexes_to_include = map { $i_num_to_index{$_} } @{ $i_num_lists[$i][$j] };
            next if (not(@indexes_to_include));
            next if (scalar(@indexes_to_include) < 5);

            push @trimmed_tracking_mat, [ map { -1 } (0 .. $#{ $tracking_mat[$i] }) ];
            @{ $trimmed_tracking_mat[$#trimmed_tracking_mat] }[@indexes_to_include] =
              @{ $tracking_mat[$i] }[@indexes_to_include];
        }
    }

    my $output_filename;
    if ($cfg{tracking_output_file} =~ m/(.*)(\..*)/) {
        $output_filename = $1 . $file_suffix . $2;
    } else {
        $output_filename = "seq_prop.csv";
    }

    my $output_file = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $output_filename);
    &output_mat_csv(\@trimmed_tracking_mat, $output_file);
}

#######################################
# Treatment Analysis
######################################
sub gather_treatment_results {
    my %treat_ad;
    #($treat_ad{'all_summary'}, $treat_ad{'all_error'}) = &gather_treatment_summary(\@tracking_mat,"Max_adhesion_signal");
    #($treat_ad{'all_summary_aver'}, $treat_ad{'all_error_aver'}) = &gather_treatment_summary;
    #($treat_ad{'across_summary'}, $treat_ad{'across_error'}) = &gather_across_treatment_summary("Max_adhesion_signal");
    ($treat_ad{'across_summary_aver'}, $treat_ad{'across_error_aver'}) = &gather_across_treatment_summary;
    my @output_mat = map [$treat_ad{'all_summary'}[$_], 
                          $treat_ad{'all_error'}[$_],
                          $treat_ad{'all_summary_aver'}[$_], 
                          $treat_ad{'all_error_aver'}[$_],
                          $treat_ad{'across_summary'}[$_],
                          $treat_ad{'across_error'}[$_],
                          $treat_ad{'across_summary_aver'}[$_],
                          $treat_ad{'across_error_aver'}[$_],
                          ], (0 .. $#{$treat_ad{'all_summary'}});
    
    &output_mat_csv(\@output_mat, catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, "treatment.csv"));
    
    return %treat_ad;
}

sub gather_across_treatment_summary {
    my @cross_indexes = &gather_cross_indexes;
    
    my @treatment_ads = grep {
        my @this_seq = @{$_};
        if (   $this_seq[$cross_indexes[0]] >= 0
            && $this_seq[$cross_indexes[1]] >= 0) {
            1;
        }
    } @tracking_mat;
    
    &output_mat_csv(\@treatment_ads, catfile($cfg{exp_results_folder}, $cfg{tracking_folder}, "treatment.csv"));

    return &gather_treatment_summary(\@treatment_ads,@_);
}

sub gather_cross_indexes {
    my @mat_index_to_i_num = sort keys %data_sets; 
    @mat_index_to_i_num = map sprintf("%0" . length($mat_index_to_i_num[-1]) . "d",$_), @mat_index_to_i_num;
    
    my @cross_indexes;

    for my $i (reverse(0 .. $#mat_index_to_i_num)) {
        if ($mat_index_to_i_num[$i] < $cfg{treatment_time}) {
            $cross_indexes[0] = $i;
            last;
        }
    }
    for my $i (0 .. $#mat_index_to_i_num) {
        if ($mat_index_to_i_num[$i] > $cfg{treatment_time}) {
            $cross_indexes[1] = $i;
            last;
        }
    }

    return @cross_indexes;
}

sub gather_treatment_summary {
    my @local_tracking_mat = @tracking_mat;
    if (scalar(@_) > 0) {
        @local_tracking_mat = @{$_[0]};
    }
    my $property = "Average_adhesion_signal";
    if (scalar(@_) > 1) {
        $property = $_[1];
    }

    my @mat_index_to_i_num = sort keys %data_sets; 
    @mat_index_to_i_num = map sprintf("%0" . length($mat_index_to_i_num[-1]) . "d",$_), @mat_index_to_i_num;
    
    my @cross_indexes = &gather_cross_indexes;
    
    #Determine the properties before treatment
    my $pix_count = 0;
    my $ad_sig_total = 0;
    my $stat = new Statistics::Descriptive::Full;
    my @per_ad_sig;
    for my $i (0 .. $cross_indexes[0]) {
        for my $j (0 .. $#local_tracking_mat) {
            my $ad_num = $local_tracking_mat[$j][$i];
            next if ($ad_num < 0); 
            
            my $area = $data_sets{$mat_index_to_i_num[$i]}{Area}[$ad_num];
            my $ad_sig = $data_sets{$mat_index_to_i_num[$i]}{$property}[$ad_num];
            
            $pix_count += $area;
            $ad_sig_total += $ad_sig;
            $stat->add_data($ad_sig);
            
            $per_ad_sig[$j][0] += 1; 
            $per_ad_sig[$j][1] += $ad_sig; 
        }
    }
    for (0 .. $#per_ad_sig) {
        $per_ad_sig[$_] = $per_ad_sig[$_][1]/$per_ad_sig[$_][0];
    }
    
    my $base_concentration = ($ad_sig_total/$stat->count);
    #print "Base Error: ", sqrt($stat->variance/$stat->count), "\n";

    #Determine the properties after treatment
    my @summary_ts;
    my @error_ts;
    my $per_ad = 0;
    my $base_ad = 0;
    my @values;
    mkpath(catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder},'raw_treatment'));
    for my $i ($cross_indexes[1] .. $#{$local_tracking_mat[0]}) {
        my $pix_count = 0;
        my $ad_sig_fold_total = 0;
        
        my $stat = new Statistics::Descriptive::Full;
        for my $j (0 .. $#local_tracking_mat) {
            my $ad_num = $local_tracking_mat[$j][$i];
            
            push @{$values[$i]}, 0 if ($ad_num < 0);
            next if ($ad_num < 0); 

            my $area = $data_sets{$mat_index_to_i_num[$i]}{Area}[$ad_num];
            my $ad_sig = $data_sets{$mat_index_to_i_num[$i]}{$property}[$ad_num];
            
            $pix_count += $area;
            
            if (defined($per_ad_sig[$j])) {
                $per_ad++;
                $ad_sig_fold_total += $ad_sig/$per_ad_sig[$j];
                push @{$values[$i]}, $ad_sig/$per_ad_sig[$j];
                $stat->add_data($ad_sig/$per_ad_sig[$j]);
            } else {
                $base_ad++;
                $ad_sig_fold_total += $ad_sig/$base_concentration;
                push @{$values[$i]}, $ad_sig/$base_concentration;
                $stat->add_data($ad_sig/$base_concentration);
            }
        }
        &output_mat_csv([$values[$i]], catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, 'raw_treatment', $i . '.csv'));
        push @summary_ts, $ad_sig_fold_total/$stat->count;
        push @error_ts, sqrt($stat->variance/$stat->count);
    }
    &output_mat_csv(\@values, catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, 'raw_treatment', 'all.csv'));
    #print "Base Ad: $base_ad, Per: $per_ad\n\n";
    
    return \@summary_ts, \@error_ts;
}

sub build_treatment_plots {
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
