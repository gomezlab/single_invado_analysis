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
if (not(defined $opt{input}) || not(-e $opt{input})) {
    print "\n\nGathering Data Files\n" if $opt{debug};
    &gather_data_sets;

    print "\n\nRemoving Excluded Images\n" if $opt{debug};
    &trim_data_sets;

    print "\n\nMaking Comparison Matrices\n" if $opt{debug};
    &make_comp_matices;
} else {
    print "\n\nGathering Data Files/Comparison Matrices from Data File\n";
    %data_sets = %{ retrieve($opt{input}) };
}

if (defined $opt{output}) {
	store \%data_sets, $opt{output};
}

print "\n\nDetermining Tracking Matrix\n" if $opt{debug};
my @tracking_mat;
&make_tracking_mat;

my $num_ad = $#{$data_sets{"001"}{Area}};

mkpath("track_seqs");
for my $start (0..$num_ad) {

	my @seq = $start;
	my @dist_seq;

	my @data_sets = sort { $a <=> $b } keys %data_sets;
	for (0..$#data_sets) {
		next if ($_ == $#data_sets);
		my $i_num = $data_sets[$_];
		push @seq, &min_index(@{${$data_sets{$i_num}{Cent_dist}}[$seq[$_]]});
		push @dist_seq, ${${$data_sets{$i_num}{Cent_dist}}[$seq[$_]]}[$seq[$_ + 1]];
	}

#for (0..$#seq) {
#	print $seq[$_],"->",$dist_seq[$_],"\n";
#}

	open OUTPUT, ">track_seqs/seq_". $start . ".csv";
	print OUTPUT join(",",@seq);
	close OUTPUT;
	open OUTPUT, ">track_seqs/dist_seq" . $start . ".csv";
	print OUTPUT join(",",@dist_seq);
	close OUTPUT;
}
###############################################################################
# Functions
###############################################################################

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
        @{ $cfg{exclude_image_nums} } = split(/\s/, $temp_line);

        shift @{ $cfg{exclude_image_nums} } if (${ $cfg{exclude_image_nums} }[0] =~ //);
    } else {
        $cfg{exclude_image_nums} = 0;
    }
    if ($opt{debug}) {
        print "Image numbers to be excluded:", join(", ", @{ $cfg{exclude_image_nums} }), "\n";
    }

    #check to see if files have been specified for the information used for the
    #tracking, if it has not, die with a message to user to specify these files
    if (not(defined $cfg{tracking_files})) {
        die
          "ERROR: The files that will be used for tracking must be specified in the config file using the name tracking_files\n";
    }

    if (defined $opt{input} && not(-e $opt{input})) {
        warn "Warning: The data file specified after \'-i\' ($opt{input}) ",
          "does not exist, input data will be gathered from standard data files.\n";
    }

    return %cfg;
}

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

sub make_comp_matices {
    my @data_keys = sort { $a <=> $b } keys %data_sets;

    for (0 .. $#data_keys) {

        #The last image can not be compared to a future image, so we skip calculations on it
        next if ($_ == $#data_keys);

        #These are the keys we will use for all the subsequent matrix creation
        my ($key_1, $key_2) = @data_keys[ $_, $_ + 1 ];

        #Gather the Centroid distance matrix
        my @x1 = @{ $data_sets{$key_1}{Centroid_x} };
        my @y1 = @{ $data_sets{$key_1}{Centroid_y} };
        my @x2 = @{ $data_sets{$key_2}{Centroid_x} };
        my @y2 = @{ $data_sets{$key_2}{Centroid_y} };
        @{ $data_sets{$key_1}{Cent_dist} } = &make_dist_mat(\@x1, \@y1, \@x2, \@y2);

        #Gather the Area difference matrix
        my @area1 = @{ $data_sets{$key_1}{Area} };
        my @area2 = @{ $data_sets{$key_2}{Area} };
        @{ $data_sets{$key_1}{Area_diff} } = &make_abs_diff_mat(\@area1, \@area2);

        #Gather the Average Signal Intensity difference matrix
        my @sig1 = @{ $data_sets{$key_1}{Average_adhesion_signal} };
        my @sig2 = @{ $data_sets{$key_2}{Average_adhesion_signal} };
        @{ $data_sets{$key_1}{Sig_diff} } = &make_abs_diff_mat(\@sig1, \@sig2);
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
            #$diff_mat[$i][$j] = sprintf('%02d',sqrt(($x1[$i] - $x2[$j])**2 + ($y1[$i] - $y2[$j])**2));
            $diff_mat[$i][$j] = sqrt(($x1[$i] - $x2[$j])**2 + ($y1[$i] - $y2[$j])**2);
        }
    }
    return @diff_mat;
}

sub make_abs_diff_mat {
    my ($ref_1, $ref_2) = @_;
    my @area1 = @$ref_1;
    my @area2 = @$ref_2;

    my @diff_mat;

    for my $i (0 .. $#area1) {
        for my $j (0 .. $#area2) {
            #$diff_mat[$i][$j] = sprintf('%02d',abs($area1[$i] - $area2[$j]));
            $diff_mat[$i][$j] = abs($area1[$i] - $area2[$j]);
        }
    }
    return @diff_mat;
}

sub make_tracking_mat {
    &initialize_tracking_mat;

}

sub initialize_tracking_mat {
    my $first_key = (sort { $a <=> $b } keys %data_sets)[0];
    print "Intializing the tracking matrix on image #: $first_key\n" if $opt{debug};

    my $num_adhesions = $#{ $data_sets{$first_key}{Area} };

    for (0 .. $num_adhesions) {
        push @{ $tracking_mat[$_] }, $_;
    }
}

sub min_index {
    my $min_i = 0;
    for (0 .. $#_) {
        if ($_[$_] < $_[$min_i]) {
            $min_i = $_;
        }
    }
    return $min_i;
}
