#!/usr/bin/perl -w

###############################################################################
#
# Global Variables and Modules
#
###############################################################################
use strict;
use File::Temp qw/ tempfile tempdir /;
use File::Spec::Functions;
use File::Path;
use Getopt::Long;
use Data::Dumper;
use Storable;
use Statistics::Descriptive;

#local libraries
use lib "../lib";
use Config::Adhesions;
use Image::Data::Collection;
use Image::Data::Writing;
use Math::R;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|config=s", "debug|d", "input|i=s", "output|o=s");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Configuration\n" if $opt{debug};

my @needed_vars = qw(data_folder results_folder exp_name single_image_folder 
    raw_data_folder general_data_files lineage_analysis_data_files 
    tracking_output_file single_lineage_props_file);
my $ad_conf = new Config::Adhesions(\%opt, \@needed_vars);
my %cfg = $ad_conf->get_cfg_hash;

###############################################################################
#
# Main Program
#
###############################################################################
print "\n\nGathering Data Files\n" if $opt{debug};

my @data_files;
push @data_files, split(/\s+/, $cfg{general_data_files});
push @data_files, split(/\s+/, $cfg{lineage_analysis_data_files});
my %data_sets = &Image::Data::Collection::gather_data_sets(\%cfg, \%opt, \@data_files);

print "\n\nRemoving Excluded Images\n" if $opt{debug};
%data_sets = &Image::Data::Collection::trim_data_sets(\%cfg, \%opt, \%data_sets);

%data_sets = &convert_data_to_units(\%data_sets,\%cfg);

print "\n\nCollecting Tracking Matrix\n" if $opt{debug};
my @tracking_mat = &Image::Data::Collection::read_in_tracking_mat(\%cfg, \%opt);

print "\n\nCollecting Individual Adhesion Properties\n" if $opt{debug};
my @single_ad_props = &gather_single_adh_props(\%cfg, \%opt);

print "\n\nGathering Adhesion Lineage Properties\n", if $opt{debug};
my %adh_lineage_props = &gather_adh_lineage_properties(\@tracking_mat,\%data_sets);

print "\n\nGathering Adhesion Property Sequences\n", if $opt{debug};
my %adh_lineage_prop_seqs = &gather_property_sequences(\@tracking_mat, \%data_sets);

print "\n\nOutputing Single Adhesion Properties\n", if $opt{debug};
&output_single_adhesion_props;

print "\n\nOutputing Adhesion Lineage Properties\n", if $opt{debug};
&output_adhesion_props;

print "\n\nOutputing Adhesion Lineage Sequence Properties\n", if $opt{debug};
&output_adhesion_prop_seqs;

print "\n\nBuilding Plots\n", if $opt{debug};
&build_r_plots;

###############################################################################
#
# Functions
#
###############################################################################

sub convert_data_to_units {
    my %data_sets = %{$_[0]};
    my %cfg = %{$_[1]};
    
    my $lin_conv_factor = $cfg{pixel_size}/$cfg{target_unit_size};
    my $sq_conv_factor = $lin_conv_factor**2;

    for my $time (keys %data_sets) {
        for my $data_type (keys %{$data_sets{$time}}) {
            if ($data_type eq "Centroid_x" || "Centroid_y" || "Centroid_dist_from_edge") {
                @{$data_sets{$time}{$data_type}} = map $lin_conv_factor*$_, @{$data_sets{$time}{$data_type}};
            }
            if ($data_type eq "Area") {
                @{$data_sets{$time}{$data_type}} = map $sq_conv_factor*$_, @{$data_sets{$time}{$data_type}};
            }
        }
    }
    return %data_sets;
}

####################################### 
#
#Single Lineage Props Collection
#
#######################################
sub gather_single_adh_props {
    my @data;
    
    my @data_types = qw(Area Average_adhesion_signal Centroid_dist_from_edge Variance_adhesion_signal);
    
    my @first_line = qw(I_num ad_num);
    push @first_line, @data_types;
    push @data, \@first_line;
    
    for my $i_num (sort keys %data_sets) {
        for my $ad_num (0 .. $#{$data_sets{$i_num}{$data_types[0]}}) {
            my @line = ($i_num,$ad_num);
            for my $i (0 .. $#data_types) {
                push @line, $data_sets{$i_num}{$data_types[$i]}[$ad_num];
            }
            push @data, \@line;
        }
    }
    return @data;
}

#######################################
#
#Adhesion Lineage Property Collection
#
#######################################
sub gather_adh_lineage_properties {
    my @tracking_mat = @{ $_[0] };
    my %data_sets    = %{ $_[1] };

    my @para_set = (\@tracking_mat, \%data_sets);

    my %props;
    $props{longevities}        = &gather_longevities(@para_set);
    $props{largest_areas}      = &gather_largest_areas(@para_set);
    $props{starting_edge_dist} = &gather_starting_dist_from_edge(@para_set);
    $props{ad_sig}             = &gather_average_ad_sig(@para_set);

    ($props{speeds}, $props{max_speeds}) = &gather_adhesion_speeds(@para_set);
    
    return %props;
}

sub gather_longevities {
    my @tracking_mat = @{ $_[0] };
    my %data_sets    = %{ $_[1] };

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

sub gather_largest_areas {
    my @tracking_mat = @{ $_[0] };
    my %data_sets    = %{ $_[1] };

    my @largest_areas;
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    for my $i (0 .. $#tracking_mat) {
        my $largest = 0;
        for my $j (0 .. $#{ $tracking_mat[$i] }) {
            my $i_num = $data_keys[$j];
            next if ($tracking_mat[$i][$j] <= -1);

            if (not defined ${ $data_sets{$i_num}{Area} }[ $tracking_mat[$i][$j] ]) {
                print "$i_num, $tracking_mat[$i][$j]";
                die;
            }
            if (   $tracking_mat[$i][$j] > -1
                && $largest < ${ $data_sets{$i_num}{Area} }[ $tracking_mat[$i][$j] ]) {
                $largest = ${ $data_sets{$i_num}{Area} }[ $tracking_mat[$i][$j] ];
            }
        }
        push @largest_areas, $largest;
    }
    return \@largest_areas;
}

sub gather_starting_dist_from_edge {
    my @tracking_mat = @{ $_[0] };
    my %data_sets    = %{ $_[1] };

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
    return \@starting_dists;
}

sub gather_adhesion_speeds {
    my @tracking_mat = @{ $_[0] };
    my %data_sets    = %{ $_[1] };

    my @speeds;
    my @max_speeds;
    my @speed_vars;
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    for my $i (0 .. $#tracking_mat) {
        my @temp = ();
        for my $j (0 .. $#{ $tracking_mat[$i] } - 1) {
            if (   $tracking_mat[$i][$j] > -1
                && $tracking_mat[$i][ $j + 1 ] > -1) {
                my $start_ad_num = $tracking_mat[$i][$j];
                my $end_ad_num   = $tracking_mat[$i][ $j + 1 ];

                my $start_x = ${ $data_sets{ $data_keys[$j] }{Centroid_x} }[$start_ad_num];
                my $start_y = ${ $data_sets{ $data_keys[$j] }{Centroid_y} }[$start_ad_num];

                my $end_x = ${ $data_sets{ $data_keys[ $j + 1 ] }{Centroid_x} }[$end_ad_num];
                my $end_y = ${ $data_sets{ $data_keys[ $j + 1 ] }{Centroid_y} }[$end_ad_num];

                push @temp, sqrt(($start_x - $end_x)**2 + ($start_y - $end_y)**2);
            }
        }
        if (@temp) {
            my $stat = Statistics::Descriptive::Full->new();
            $stat->add_data(@temp);

            push @speeds,     $stat->mean();
            push @speed_vars, $stat->variance();
            push @max_speeds, $stat->max();
        } else {
            push @speeds,     "NaN";
            push @speed_vars, "NaN";
            push @max_speeds, "NaN";
        }
    }

    return \@speeds, \@max_speeds;
}

sub gather_average_ad_sig {
    my @tracking_mat = @{ $_[0] };
    my %data_sets    = %{ $_[1] };

    my @pax_sig;
    my @pax_sig_vars;
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    for my $i (0 .. $#tracking_mat) {
        my @temp = ();
        for my $j (0 .. $#{ $tracking_mat[$i] }) {
            if ($tracking_mat[$i][$j] > -1) {
                my $ad_num = $tracking_mat[$i][$j];
                push @temp, ${ $data_sets{ $data_keys[$j] }{Average_adhesion_signal} }[$ad_num];
            }
        }
        my $stat = Statistics::Descriptive::Full->new();
        $stat->add_data(@temp);

        push @pax_sig,      $stat->mean();
        push @pax_sig_vars, $stat->variance();
    }

    return \@pax_sig;
}

####################################### 
#
# Adhesion Lineage Property Sequence 
#
######################################
sub gather_property_sequences {
    my @tracking_mat = @{ $_[0] };
    my %data_sets    = %{ $_[1] };

    my %seqs;
    $seqs{area}{increasing} = &gather_patterned_property_seqs("Area",\@tracking_mat,\%data_sets,\&greater_or_equal);
    $seqs{area}{decreasing} = &gather_patterned_property_seqs("Area",\@tracking_mat,\%data_sets,\&less_or_equal);

    $seqs{ad_sig}{increasing} = &gather_patterned_property_seqs("Average_adhesion_signal",\@tracking_mat,\%data_sets,\&greater_or_equal);
    $seqs{ad_sig}{decreasing} = &gather_patterned_property_seqs("Average_adhesion_signal",\@tracking_mat,\%data_sets,\&less_or_equal);
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
        for my $j (0 .. $#{$tracking_mat[$i]}) {
            next if ($tracking_mat[$i][$j] <= -1);
            push @i_num_seq, $image_nums[$j];
            push @data_seq, $data_sets{$image_nums[$j]}{$property}[$tracking_mat[$i][$j]];
        }

        my @indexes = &gather_indexes_seqs_by_func(\@data_seq,$func_ref);
        for my $j (0 .. $#indexes) {
            push @{$image_num_seqs[$i]}, [@i_num_seq[@{$indexes[$j]}]];
        }
        if (not(@indexes)) {
            push @{$image_num_seqs[$i]}, ();
        }
    }
    return \@image_num_seqs;
}

sub gather_indexes_seqs_by_func {
    my @set = @{ $_[0] };
    my $deciding_func = $_[1];

    my @deciding_func_true = map { $deciding_func->($set[$_], $set[ $_ - 1 ]) } (1 .. $#set);
    unshift @deciding_func_true, 0;

    my @index_seqs;
    @{$index_seqs[0]} = ();
    
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
        my @this_seq = @{$index_seqs[$_]};
        if ($#this_seq > -1) {
            unshift @this_seq, $this_seq[0] - 1;
            @{$index_seqs[$_]} = @this_seq;
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

####################################### 
#
#Output Adhesion Lineage Properties 
#
#######################################
sub output_single_adhesion_props {
    if (not(-e catdir($cfg{exp_results_folder}, $cfg{lineage_props_folder}))) {
        mkpath(catdir($cfg{exp_results_folder}, $cfg{lineage_props_folder}));
    }

    &output_all_single_ad_data;
}

sub output_all_single_ad_data {
    my $output_file = catfile($cfg{exp_results_folder}, $cfg{lineage_props_folder}, $cfg{individual_adhesions_props_file});
    &Image::Data::Writing::output_mat_csv(\@single_ad_props,$output_file);
}

####################################### 
#
#Output Adhesion Lineage Properties 
#
#######################################
sub output_adhesion_props {
    if (not(-e catdir($cfg{exp_results_folder}, $cfg{lineage_props_folder}))) {
        mkpath(catdir($cfg{exp_results_folder}, $cfg{lineage_props_folder}));
    }

    &output_all_prop_data;
}

sub output_all_prop_data {
    my @longevities    = @{ $adh_lineage_props{longevities} };
    my @largest_areas  = @{ $adh_lineage_props{largest_areas} };
    my @starting_dists = @{ $adh_lineage_props{starting_edge_dist} };
    my @speeds         = @{ $adh_lineage_props{speeds} };
    my @max_speeds     = @{ $adh_lineage_props{max_speeds} };
    my @ad_sig         = @{ $adh_lineage_props{ad_sig} };

    my @all_data =
      map { [ $longevities[$_], $largest_areas[$_], $starting_dists[$_], $speeds[$_], $max_speeds[$_], $ad_sig[$_] ] }
      (0 .. $#longevities);
    
    unshift @all_data, [qw(longevity largest_area s_dist_from_edge speed max_speed ad_sig)];
    
    my $output_file = catfile($cfg{exp_results_folder}, $cfg{lineage_props_folder}, $cfg{single_lineage_props_file});
    &Image::Data::Writing::output_mat_csv(\@all_data,$output_file);
}

####################################### 
#
#Output Adhesion Lineage Prop Sequences 
#
#######################################
sub output_adhesion_prop_seqs {
    if (not(-e catdir($cfg{exp_results_folder}, $cfg{lineage_props_folder}))) {
        mkpath(catdir($cfg{exp_results_folder}, $cfg{lineage_props_folder}));
    }
    
    &output_sequence_trimmed_mat(\@{$adh_lineage_prop_seqs{area}{increasing}},"_increasing_area");
}

sub output_sequence_trimmed_mat {
    my @i_num_lists = @{$_[0]};
    my $file_suffix = $_[1];

    my @i_nums = sort keys %data_sets;
    my %i_num_to_index = map { $i_nums[$_] => $_ } (0 .. $#i_nums);

    my @trimmed_tracking_mat;

    for my $i (0 .. $#i_num_lists) {
        for my $j (0 .. $#{$i_num_lists[$i]}) {
            my @indexes_to_include = map { $i_num_to_index{$_} } @{$i_num_lists[$i][$j]};
            next if (not(@indexes_to_include));
            next if (scalar(@indexes_to_include) < 5);

            push @trimmed_tracking_mat, [map {-1} (0 .. $#{$tracking_mat[$i]})];
            @{$trimmed_tracking_mat[$#trimmed_tracking_mat]}[@indexes_to_include] = @{$tracking_mat[$i]}[@indexes_to_include];
        }
    }
    
    my $output_filename;
    if ($cfg{tracking_output_file} =~ m/(.*)(\..*)/) {
        $output_filename = $1 . $file_suffix . $2;
    } else { 
        $output_filename = "seq_prop.csv";
    }

    my $output_file = catfile($cfg{exp_results_folder}, $cfg{lineage_props_folder}, $output_filename);
    &Image::Data::Writing::output_mat_csv(\@trimmed_tracking_mat,$output_file);
}


####################################### 
#
#Plotting 
#
#######################################

sub build_r_plots {
    my @r_code;
    
    my $data_dir = catdir($cfg{exp_results_folder},$cfg{lineage_props_folder});
    my $plot_dir = catdir($data_dir,$cfg{plot_folder});
    
    my $xy_default = "pch=19,cex=0.1";
    my $pdf_default = "pointsize=19";
    my @xy_plots = ({xy => "lineages\$longevity,lineages\$s_dist_from_edge",
                     main => "",
                     xlab => "'Longevity (min)'",
                     ylab => "expression(paste('Starting Distance from Edge (', mu, 'm)'))",
                     file_name => "longev_vs_s_dist.pdf",
                     plot_para => $xy_default,
                     pdf_para => $pdf_default,},
                     {xy => "lineages\$longevity,lineages\$ad_sig",
                     main => "",
                     xlab => "'Longevity (min)'",
                     ylab => "'Paxillin Concentration (AU)'",
                     file_name => "longev_vs_pax.pdf",
                     plot_para => $xy_default,
                     pdf_para => $pdf_default,},
                     {xy => "adhesions\$Area,adhesions\$Average_adhesion_signal",
                     main => "",
                     xlab => "expression(paste('Area (', mu, m^2, ')'))",
                     ylab => "'Paxillin Concentration (AU)'",
                     file_name => "area_vs_pax.pdf",
                     plot_para => $xy_default,
                     pdf_para => $pdf_default,},
                     {xy => "adhesions\$Area,adhesions\$Centroid_dist_from_edge",
                     main => "",
                     xlab => "expression(paste('Area (', mu, m^2, ')'))",
                     ylab => "expression(paste('Distance from Edge (', mu, 'm)'))",
                     file_name => "area_vs_dist.pdf",
                     plot_para => $xy_default,
                     pdf_para => $pdf_default,},
                     {xy => "adhesions\$Average_adhesion_signal,adhesions\$Centroid_dist_from_edge",
                     main => "",
                     xlab => "expression(paste('Area (', mu, m^2, ')'))",
                     ylab => "expression(paste('Distance from Edge (', mu, 'm)'))",
                     file_name => "sig_vs_dist.pdf",
                     plot_para => $xy_default,
                     pdf_para => $pdf_default,},
                   );

    mkpath($plot_dir);

    #Read in data
    push @r_code, "lineages = read.table('$data_dir/$cfg{single_lineage_props_file}',header=T,sep=',');\n";
    push @r_code, "adhesions = read.table('$data_dir/$cfg{individual_adhesions_props_file}',header=T,sep=',');\n";
    
    #Build the plots
    foreach (@xy_plots) {
        my %parameters = %{$_};
        my $output_file = catfile($plot_dir,$parameters{file_name});
        
        push @r_code, "pdf('$output_file',$parameters{pdf_para})\n";
        #push @r_code, "par(oma=c(0,0,0,0))\n"; push @r_code,
        #"par(omi=c(0,0,0,0))\n";
        push @r_code, "par(mar=c(4,4,0.5,0.5),bty='n')\n";
        push @r_code, "plot($parameters{xy},xlab=$parameters{xlab},ylab=$parameters{ylab},$parameters{plot_para})\n";
        push @r_code, "dev.off();\n";
    }
    &Math::R::execute_commands(\@r_code);
}
