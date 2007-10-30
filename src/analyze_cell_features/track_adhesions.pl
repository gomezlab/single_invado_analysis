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
GetOptions(\%opt, "cfg=s","debug|d");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my %cfg = &get_config;


###############################################################################
# Main Program
###############################################################################



###############################################################################
# Functions
###############################################################################

sub get_config {
	if ($opt{debug}) {
		print "Collecting Configuration:\n";
	}

	my %cfg = ParseConfig(-ConfigFile => $opt{cfg},
						  -MergeDuplicateOptions => 1,
	);
	

	if (defined $cfg{exclude_file}) {
		open EX_INPUT, $cfg{exclude_file} or die "Can't open the specified exclude file: $cfg{exclude_file}";
		my $temp_line = <EX_INPUT>;
		chomp($temp_line);
		@{$cfg{exclude_image_nums}} = split(/\s/,$temp_line);
		close EX_INPUT;
	}
	else {
		$cfg{exclude_image_nums} = 0;
	}

	if ($opt{debug}) {
		print "Image numbers to be excluded:\n";
		print join(" ",@{$cfg{exclude_image_nums}}),"\n";
	}

	return %cfg;
}
