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

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg=s", "debug|d");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Configuration:\n" if $opt{debug};
my %cfg = &get_config;

###############################################################################
# Main Program
###############################################################################

print "Gathering Data Files:\n" if $opt{debug};
my %data_sets = &gather_data_sets;

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

    #check to see if files have been specified for the information used for the
    #tracking, if it has not, die with a message to user to specify these files
    if (defined $cfg{tracking_files}) {

    } else {
        die
          "The files that will be used for tracking must be specified in the config file using the name tracking_files\n";
    }

    if ($opt{debug}) {
        print "Image numbers to be excluded:\n";
        print join(" ", @{ $cfg{exclude_image_nums} }), "\n";
    }

    return %cfg;
}

sub gather_data_sets {
    my %data;

    my @folders = <$cfg{results_folder}/*/$cfg{raw_data_folder}>;

    foreach my $folder (@folders) {
        foreach my $file (@{ $cfg{tracking_files} }) {
            print $file, " ";
        }
    }
}

sub gather_data_from_matlab_file {
    my ($file) = @_;

    open INPUT, "$file" or die "Problem opening $file";
    my @in = split("   ", <INPUT>);
    close INPUT;

    shift @in if $in[0] =~ //;
}
