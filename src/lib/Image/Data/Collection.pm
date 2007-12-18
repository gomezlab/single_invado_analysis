#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use warnings;
use File::Spec;
use Text::CSV::Simple;

package Image::Data::Collection;
###############################################################################
# Functions
###############################################################################

#######################################
# Collection/Verification
#######################################
sub gather_data_sets {
    my %cfg        = %{ $_[0] };
    my %opt        = %{ $_[1] };
    my @data_files = @{ $_[2] };

    my %data_sets;

    my @folders     = <$cfg{individual_results_folder}/*/$cfg{raw_data_folder}>;
    my $image_count = 0;

    foreach my $this_folder (@folders) {
        my $i_num;
        if ($this_folder =~ /$cfg{individual_results_folder}\/(.*)\/$cfg{raw_data_folder}/) {
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
                    delete $data_sets{$i_num}{$file};
                }
                if ($file eq "Area" || $file eq "Centroid_dist_from_edge") {
                    @{ $data_sets{$i_num}{$file} } = map sprintf("%f", $_), @{ $data_sets{$i_num}{$file} };
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

    &check_data_set_lengths(\%data_sets);

    return %data_sets;
}

sub gather_data_from_matlab_file {
    my ($file) = @_;

    my $parser = Text::CSV::Simple->new;
    my @data   = $parser->read_file($file);

    die "Found two lines of data in $file" if scalar(@data) > 1;

    return @{ $data[0] };
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

sub check_data_set_lengths {
    my %data_sets = %{ $_[0] };

    my %data_sets_length;
    for my $key (keys %data_sets) {
        for my $data_type (keys %{ $data_sets{$key} }) {
            $data_sets_length{$key}{$data_type} = scalar(@{ $data_sets{$key}{$data_type} });
        }
    }

    for my $key (sort { $a <=> $b } keys %data_sets_length) {
        my $all_same       = 1;
        my @data_type_keys = keys %{ $data_sets_length{$key} };
        my $length         = $data_sets_length{$key}{ $data_type_keys[0] };
        for my $data_type (keys %{ $data_sets_length{$key} }) {
            if ($data_sets_length{$key}{$data_type} != $length) {
                $all_same = 0;
            }
        }
        if (not $all_same) {
            warn "Data set lengths do not match in image number $key:\n",
              map { "    $_ - $data_sets_length{$key}{$_}\n" } keys %{ $data_sets_length{$key} };
        }
    }
}

########################################
# Other
#######################################
sub trim_data_sets {
    my %cfg       = %{ $_[0] };
    my %opt       = %{ $_[1] };
    my %data_sets = %{ $_[2] };

    my @excluded_nums;
    for my $ex_num (@{ $cfg{exclude_image_nums} }) {
        for my $this_num (sort { $a <=> $b } keys %data_sets) {
            if ($ex_num == $this_num) {
                delete $data_sets{$this_num};
                push @excluded_nums, $ex_num;
            }
        }
    }
    if ($opt{debug}) {
        if (scalar(@excluded_nums) != 0) {
            print "Image number removed from further Dataset: ", join(", ", @excluded_nums), "\n";
        } else {
            print "Image number removed from further Dataset: none\n";
        }
    }

    return %data_sets;
}

sub read_in_tracking_mat {
    my %cfg = %{ $_[0] };
    my %opt = %{ $_[1] };
    my $file = File::Spec->catdir($cfg{results_folder}, $cfg{exp_name}, $cfg{tracking_output_file});

    my $parser = Text::CSV::Simple->new;
    my @tracking_mat = $parser->read_file($file);

    print "Gathered ", scalar(@tracking_mat), " lineages, with ", scalar(@{ $tracking_mat[0] }),
      " timepoints from file ", $cfg{tracking_output_file}, ".\n"
      if $opt{debug};

    return @tracking_mat;
}

1;
