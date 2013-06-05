#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################
use lib "../lib";

use strict;
use File::Path;
use File::Basename;
use File::Find;
use Getopt::Long;
use Data::Dumper;
use File::Spec::Functions;
use Benchmark;
use POSIX;

use Config::ImageSet qw(ParseConfig);

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
$opt{min_ad_size} = 5;
GetOptions(\%opt, "cfg=s", "debug|d", "min_ad_size=s", "opacity", "white_background") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my %cfg = ParseConfig(\%opt);

###############################################################################
#Main Program
###############################################################################

our @files;
find(\&include_in_vis, catdir($cfg{individual_results_folder}));
my @colors = &build_jet_color_map(scalar(@files));
my @hex_colors = &convert_rbg_to_hex(@colors);
@files = sort @files;

my @bmp_files = &build_bitmaps(@files);
print "Done Converting to BMP\n" if $opt{debug};

my @svg_files = &convert_using_potrace(@bmp_files);
print "Done Converting to SVG\n" if $opt{debug};

my @svg_header = &get_svg_header($svg_files[0]);

if ($opt{white_background}) {
	@svg_header = &convert_to_white_background(@svg_header);
}

my @full_svg_file = &build_full_svg_file(@svg_files);

unlink(@bmp_files);
unlink(@svg_files);

my $output_file = catfile($cfg{exp_results_folder}, $cfg{movie_output_folder}, "ghost.svg");

open SVG_OUT, ">$output_file";
print SVG_OUT @svg_header;
print SVG_OUT @full_svg_file;
print SVG_OUT "</g>\n";
print SVG_OUT '</svg>';
close SVG_OUT;

my $output_png = catfile($cfg{exp_results_folder}, $cfg{movie_output_folder}, "ghost.png");
my $output_png_small = catfile($cfg{exp_results_folder}, $cfg{movie_output_folder}, "ghost_small.png");

system "convert -density 100x100 $output_file $output_png_small";
system "convert -density 300x300 $output_file $output_png";

###############################################################################
#Functions
###############################################################################
sub include_in_vis {
    if ($File::Find::name =~ /adhesions_binary.png/) {
        push @files, $File::Find::name;
    }
}

sub build_bitmaps {
    my @files = @_;
    my @bmp_files;
    for (@files) {
        my $bmp_file = $_;
        $bmp_file =~ s/\.png/\.bmp/;
        push @bmp_files, $bmp_file;
        system "convert $_ $bmp_file";
    }
    return @bmp_files;
}

sub convert_using_potrace {
    my @files = @_;
    my @svg_files;
    for (0 .. $#files) {
        my $bmp_file = $files[$_];
        my $color = $hex_colors[$_];
        my $svg_file = $bmp_file;
     	   
        $svg_file =~ s/\.bmp/\.svg/;
        push @svg_files, $svg_file;
        system "potrace -t $opt{min_ad_size} -s --fillcolor=#$color $bmp_file\n";
    }
    return @svg_files;
}

sub get_svg_header {
    my $file = $_[0];

    my @header;
    my $path_count = 0;

    open SVG_FILE, $file;
    while (<SVG_FILE>) {
        if (/\<path/) {
            $path_count++;
        }
        if ($path_count < 2) {
            push @header, $_;
        }
    }
    close SVG_FILE;

    return @header;
}

sub convert_to_white_background {
	my @header = @_;
	
	my @processed_header;
	foreach (@header) {
		if ($_ =~ /path fill/) {
			$_ =~ s/#000000/#ffffff/g;
			push @processed_header, $_;
		} else {
			push @processed_header, $_;
		}
	}
	return @processed_header;
}

sub build_full_svg_file {
    my @svg_files = @_;
    @svg_files = sort @svg_files;

    my @svg_path_data;

    for (0 .. $#svg_files) {
    #for (0 .. 0) {
        my $file = $svg_files[$_];
        my $opacity = ($#svg_files - $_)/$#svg_files;
        my $path_count = 0;
        my $out_of_path = 0;
        open SVG_FILE, $file;
        for my $line (<SVG_FILE>) {
            if ($line =~ /\<path/) {
                $path_count++;
            }
            if ($line =~ /\<\/g\>/) {
                $out_of_path = 1;
            }
            if ($path_count > 1 && not($out_of_path)) {
                if ($opt{opacity}) {
                    $line =~ s/fill/opacity="$opacity" fill/;
                }
                push @svg_path_data, $line;
            }
        }
        close SVG_FILE;
    }

    return @svg_path_data;
}

sub build_jet_color_map {
    #This function is cribbed from the MATLAB jet function
    my $needed_colors = $_[0];

    my $n = ceil($needed_colors/4);
    #print Dumper($n);
    my @u;
    for (1 .. $n) {
        push @u, $_/$n;
    }
    for (1 .. ($n - 1)) {
        push @u, 1;
    }
    for (reverse(1 .. $n)) {
        push @u, $_/$n;
    }
    #print Dumper(\@u);

    my @g;
    my $add_factor = ceil($n/2);
    if ($needed_colors % 4 == 1) {
        $add_factor--;
    }
    for (1 .. scalar(@u)) {
        push @g, $_ + $add_factor;   
    }

    #print Dumper(\@g);

    my @r;
    my @b;
    for (@g) {
        if (($_ + $n) <= $needed_colors) {
            push @r, $_ + $n - 1;
        }
        if (($_ - $n) >= 1) {
            push @b, $_ - $n - 1;
        }
    }

    #print Dumper(\@r);
    #print Dumper(\@b);

    for (0 .. $#g) {
        $g[$_]--;
    }

    my @c_map_r;
    @c_map_r[@r] = @u[0 .. scalar(@r)];

    #print Dumper(\@c_map_r);
    
    my @c_map_g;
    @c_map_g[@g] = @u[0 .. scalar(@g)];
    
    #print Dumper(\@c_map_g);
    
    my @c_map_b;
    @c_map_b[@b] = @u[($#u - scalar(@b) + 1) .. $#u];
    
    my @c_map;
    for my $i (0 .. $needed_colors - 1) {
        for my $j (0 .. 2) {
            if ($j == 0) {
                if (defined $c_map_r[$i]) {
                    $c_map[$i][$j] = $c_map_r[$i];
                } else {
                    $c_map[$i][$j] = 0;
                }
            } elsif ($j == 1) {
                if (defined $c_map_g[$i]) {
                    $c_map[$i][$j] = $c_map_g[$i];
                } else {
                    $c_map[$i][$j] = 0;
                }
            } elsif ($j == 2) {
                if (defined $c_map_b[$i]) {
                    $c_map[$i][$j] = $c_map_b[$i];
                } else {
                    $c_map[$i][$j] = 0;
                }
            }
        }
        #print join(",", @{$c_map[$i]}), "\n";
        #print $c_map[$i][2], "\n";
    }
    return @c_map;
}

sub convert_rbg_to_hex {
    my @colors = @_;
    
    my @hex_colors;
    
    for (@colors) {
        my @this_color = @{$_};
        #my $in_hex = sprintf('%0x%0x%0x', round($this_color[0] * 255), round($this_color[1]*255), round($this_color[2]*255));;
        my $in_hex;
        for (@this_color) {
            my $hex_num = sprintf('%x',round($_ * 255));
            if (length($hex_num) == 1) {
                $hex_num = "0" . $hex_num;
            }
            $in_hex .= $hex_num;
        }
        push @hex_colors, $in_hex;
    }
    return @hex_colors;
}

sub round {
    my($number) = shift;
    return int($number + .5);
}
