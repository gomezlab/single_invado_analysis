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
&output_single_adhesion_props_density;
&build_single_ad_plots;

print "\n\nCreating Adhesion Lineage Property Plots\n", if $opt{debug};
my %ad_lineage_props = &gather_ad_lineage_properties;
&output_adhesion_lineage_props;
&build_lineage_plots;

print "\n\nGathering Adhesion Property Sequences\n", if $opt{debug};
my %ad_lineage_prop_seqs = &gather_property_sequences(\@tracking_mat, \%data_sets);

#&output_adhesion_prop_seqs;

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
            } elsif (grep $data_type eq $_, qw(Average_adhesion_signal Variance_adhesion_signal PixelIdxList)) {

                #This is the arbitrary units place, don't do any unit
                #conversion
            } else {
                die "Unable to determine the appropriate conversion factor for $data_type.";
            }
        }
    }
    return %data_sets;
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
        $cell_mask_img->read(file => "$cell_mask_files[$i]") or die;

        my @vals;
        my @intra_cellular;
        for my $y (0 .. $focal_img->getheight - 1) {
            my @gray = $focal_img->getsamples(y => $y, type => 'float');
            push @vals, @gray;

            my @cell_mask = $cell_mask_img->getsamples(y => $y);
            my @in_cell_indexes = grep $cell_mask[$_], (0 .. $#cell_mask);
            push @intra_cellular, @gray[@in_cell_indexes];
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

    my @data_types =
      qw(Area Average_adhesion_signal Centroid_dist_from_edge Centroid_dist_from_center Variance_adhesion_signal);

    my @first_line = qw(I_num ad_num);
    push @first_line, @data_types;
    push @data,       \@first_line;

    my @i_num_list = sort keys %data_sets;
    my %i_num_to_col = map { $i_num_list[$_] => $_ } (0 .. $#i_num_list);

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
            for my $i (0 .. $#data_types) {
                push @line, $data_sets{$i_num}{ $data_types[$i] }[$ad_num];
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

sub output_single_adhesion_props_density {
    my $output_file = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, 'density.csv');

    &collect_density;
}

sub collect_density {
    my $box_num = 100;

    my $prop_1 = "Area";
    my $prop_2 = "Centroid_dist_from_center";

    my $col_num_1 = (grep $single_ad_props[0][$_] eq $prop_1, (0 .. $#{ $single_ad_props[0] }))[0];
    my $col_num_2 = (grep $single_ad_props[0][$_] eq $prop_2, (0 .. $#{ $single_ad_props[0] }))[0];

    my @prop_1 = map $single_ad_props[$_][$col_num_1], (1 .. $#single_ad_props);
    my @prop_2 = map $single_ad_props[$_][$col_num_2], (1 .. $#single_ad_props);

    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@prop_1);
}

sub build_single_ad_plots {
    my @r_code;

    my $data_dir = catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder});
    my $plot_dir = catdir($data_dir, $cfg{plot_folder});

    my $xy_default = "pch=19,cex=0.4";

    #my $xy_default = "pch=19";
    my $pdf_default = "height=12, width=12, pointsize=24";
    my @xy_plots    = (
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
    $props{longevities}               = &gather_longevities;
    $props{Area}                      = &gather_prop_seq("Area");
    $props{largest_areas}             = &gather_largest_areas($props{Area});
    $props{Centroid_dist_from_center} = &gather_prop_seq("Centroid_dist_from_center");
    $props{starting_center_dist}      = &gather_starting_dist($props{Centroid_dist_from_center});
    $props{Centroid_dist_from_edge}   = &gather_prop_seq("Centroid_dist_from_edge");
    $props{starting_edge_dist}        = &gather_starting_dist($props{Centroid_dist_from_edge});
    $props{Average_adhesion_signal}   = &gather_prop_seq("Average_adhesion_signal");
    $props{ad_sig}                    = &gather_average_ad_sig($props{Average_adhesion_signal});
    $props{merge}                     = &gather_merge_count;
    $props{All_speeds}                = &gather_adhesion_speeds;

    ($props{average_speeds}, $props{variance_speeds}, $props{max_speeds}) 
      = &gather_speed_props($props{All_speeds});

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

            push @{ $prop_vals[$i] }, $default_val if ($ad_num <= -1);
            next if ($ad_num <= -1);

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
    my @data_keys = sort keys %data_sets;
    for my $i (0 .. $#tracking_mat) {
        for my $j (0 .. $#{ $tracking_mat[$i] } - 1) {
            my $start_i_num = $data_keys[$j];
            my $end_i_num   = $data_keys[ $j + 1 ];

            my $start_ad_num = $tracking_mat[$i][$j];
            my $end_ad_num   = $tracking_mat[$i][ $j + 1 ];

            if ($start_ad_num > -1 && $end_ad_num > -1) {
                my $start_x = ${ $data_sets{ $data_keys[$j] }{Centroid_x} }[$start_ad_num];
                my $start_y = ${ $data_sets{ $data_keys[$j] }{Centroid_y} }[$start_ad_num];

                my $end_x = ${ $data_sets{ $data_keys[ $j + 1 ] }{Centroid_x} }[$end_ad_num];
                my $end_y = ${ $data_sets{ $data_keys[ $j + 1 ] }{Centroid_y} }[$end_ad_num];
                
                my $speed = sqrt(($start_x - $end_x)**2 + ($start_y - $end_y)**2);
                
                if ($speed > 10) {
                    #print $i, ",";
                }

                push @{ $speed[$i] }, $speed;
            } else {
                push @{ $speed[$i] }, "NaN";
            }
        }
    }

    return \@speed;
}

sub gather_speed_props {
    my @speed = @{ $_[0] };

    my @av_speeds;
    my @var_speeds;
    my @max_speeds;

    for my $i (0 .. $#speed) {
        my $stat = Statistics::Descriptive::Full->new();
        
        my @these_speeds = grep $_ ne "NaN", @{$speed[$i]};

        $stat->add_data(@these_speeds) if (scalar(@these_speeds) != 0);
        
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
    #print $stat->max();
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

    my @longevities           = @{ $ad_lineage_props{longevities} };
    my @largest_areas         = @{ $ad_lineage_props{largest_areas} };
    my @starting_edge_dists   = @{ $ad_lineage_props{starting_edge_dist} };
    my @starting_center_dists = @{ $ad_lineage_props{starting_center_dist} };
    my @merge_counts          = @{ $ad_lineage_props{merge} };
    my @average_speeds        = @{ $ad_lineage_props{average_speeds} };
    my @max_speeds            = @{ $ad_lineage_props{max_speeds} };
    my @ad_sig                = @{ $ad_lineage_props{ad_sig} };

    my @all_data =
      map {
        [
            $longevities[$_],           $largest_areas[$_],  $starting_edge_dists[$_],
            $starting_center_dists[$_], $average_speeds[$_], $max_speeds[$_],
            $ad_sig[$_]
        ]
      } (0 .. $#longevities);

    unshift @all_data, [qw(longevity largest_area s_dist_from_edge s_dist_from_center speed max_speed ad_sig)];

    my $output_file = catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{lineage_summary_props_file});
    &output_mat_csv(\@all_data, $output_file);

    my @single_output_props =
      qw(All_speeds Area Centroid_dist_from_edge Centroid_dist_from_center Average_adhesion_signal);
    foreach (@single_output_props) {
        my $output_file =
          catfile($cfg{exp_results_folder}, $cfg{adhesion_props_folder}, $cfg{lineage_ts_folder}, $_ . ".csv");
        &output_mat_csv($ad_lineage_props{$_}, $output_file);
    }
}

sub build_lineage_plots {
    my @r_code;

    my $data_dir = catdir($cfg{exp_results_folder}, $cfg{adhesion_props_folder});
    my $plot_dir = catdir($data_dir, $cfg{plot_folder});

    my $xy_default = "pch=19,cex=0.4";

    #my $xy_default = "pch=19";
    my $pdf_default = "height=12, width=12, pointsize=24";
    my @xy_plots    = (
        {
            xy        => "lineages\$longevity,lineages\$s_dist_from_edge",
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
