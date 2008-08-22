#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use lib "../lib";
use lib "../lib/perl";

use strict;
use File::Path;
use File::Basename;
use File::Spec::Functions;
use Getopt::Long;
use Data::Dumper;
use Storable qw(nstore);
use Text::CSV;
use IO::File;

use Config::Adhesions;
use Image::Data::Collection;
use Text::CSV::Simple::Extra;
use Emerald;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
$opt{output} = "data.stor";
GetOptions(\%opt, "cfg|config=s", "debug|d", "output|o=s", "image_num=s", 
                  "emerald|e", "emerald_debug|e_d") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Configuration\n" if $opt{debug};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg = $ad_conf->get_cfg_hash;

###############################################################################
#Main Program
###############################################################################

my %data_sets;
if ($opt{emerald} || $opt{emerald_debug}) {
    my @image_nums = &Image::Data::Collection::gather_sorted_image_numbers(\%cfg);
    my $error_folder = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'tracking_data');
    mkpath($error_folder);
    
    my %emerald_opt = ("folder" => $error_folder);

    my @commands;
    foreach (@image_nums) {
        #$0 - the name of the program currently running, used to protect against
        #future file name changes
        push @commands, "$0 -cfg $opt{cfg} -o $opt{output} -image_num $_";
    }
    @commands = &Emerald::create_general_LSF_commands(\@commands,\%emerald_opt);
    if ($opt{emerald_debug}) {
        print join("\n", @commands);
    } else {
        &Emerald::send_LSF_commands(\@commands);
    }
    exit;
} else {
    print "\n\nGathering Data Files\n" if $opt{debug};

    my @data_files;
    push @data_files, @{$cfg{general_data_files}};
    push @data_files, @{$cfg{tracking_files}};

    %data_sets = &Image::Data::Collection::gather_data_sets(\%cfg, \%opt, \@data_files);

    print "\n\nMaking Comparison Matrices\n" if $opt{debug};
    &make_comp_matices;
}

###############################################################################
#Functions
###############################################################################

#######################################
# Process Data Sets
#######################################
sub make_comp_matices {
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    for (0 .. $#data_keys) {
        if (exists $opt{image_num} && $data_keys[$_] != $opt{image_num}) {
            next;
        }

        #The last image can not be compared to a future image, so we skip
        #calculations on it, but still save the image data if the output
        #option is specified
        if ($_ == $#data_keys) {
            if (defined $opt{output}) {
                my $key_1 = $data_keys[$_];
                nstore \%{ $data_sets{$key_1} }, catfile($cfg{individual_results_folder}, $key_1, $opt{output});
                delete $data_sets{$key_1};
            }
            next;
        }

        #These are the keys we will use for all the subsequent matrix creation
        my ($key_1, $key_2) = @data_keys[ $_, $_ + 1 ];

        print "Working on image Number: $key_1 - " if $opt{debug};

        #Gather the Centroid distance matrix
        my @x1 = @{ $data_sets{$key_1}{Centroid_x} };
        my @y1 = @{ $data_sets{$key_1}{Centroid_y} };
        my @x2 = @{ $data_sets{$key_2}{Centroid_x} };
        my @y2 = @{ $data_sets{$key_2}{Centroid_y} };
        @{ $data_sets{$key_1}{Cent_dist} } = &make_dist_mat(\@x1, \@y1, \@x2, \@y2);
        print "Cent_dist Collected - " if $opt{debug};

        #Gather the Area difference matrix
        my @area1 = @{ $data_sets{$key_1}{Area} };
        my @area2 = @{ $data_sets{$key_2}{Area} };
        @{ $data_sets{$key_1}{Area_diff} } = &make_abs_diff_mat(\@area1, \@area2);
        print "Area_diff Collected - " if $opt{debug};

        #Gather the Pixel Similarity matrix
        my @pix_id1 = @{ $data_sets{$key_1}{PixelIdxList} };
        my @pix_id2 = @{ $data_sets{$key_2}{PixelIdxList} };
        @{ $data_sets{$key_1}{Pix_sim} } = &calc_pix_sim(\@pix_id1, \@pix_id2);

        print "Pix_sim Collected" if $opt{debug};
        print "\r"                if $opt{debug};

        if (defined $opt{output}) {
            nstore \%{ $data_sets{$key_1} }, catfile($cfg{individual_results_folder}, $key_1, $opt{output});
            delete $data_sets{$key_1};
        }
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
                if ($a > 1) {
                    die;
                }
                $a;
            } @matching_list;
            push @temp_sim, $match_count / scalar(@pix_list);
        }
        push @sim_percents, \@temp_sim;
    }
    return @sim_percents;
}
