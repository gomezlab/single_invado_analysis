#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Path;
use File::Find;
use File::Spec::Functions;
use Getopt::Long;
use IO::File;
use Benchmark;
use Data::Dumper;

use lib "../lib";
use Config::Adhesions;
use Image::Data::Collection;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|config=s", "debug|d") or die;
die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Overall Configuration\n\n" if $opt{debug};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg = $ad_conf->get_cfg_hash;

###############################################################################
# Main Program
###############################################################################
my ($t1, $t2);

#Collecting Visualizations
chdir "../visualize_cell_features";

#Build Movies 
my @image_numbers = &Image::Data::Collection::gather_sorted_image_numbers(\%cfg);
my $image_num_length = length(scalar(@image_numbers));

my @commands = (
	# "ffmpeg -v 0 -y -r $cfg{movie_frame_rate} -i $cfg{individual_results_folder}/%0" . $image_num_length . "d/gel.png -sameq $cfg{exp_results_folder}/gel.mov 2>&1",
	# "ffmpeg -v 0 -y -r $cfg{movie_frame_rate} -i $cfg{individual_results_folder}/%0" . $image_num_length . "d/binary_shift.png -sameq $cfg{exp_results_folder}/binary_shift.mov 2>&1",
	# "ffmpeg -v 0 -y -r $cfg{movie_frame_rate} -i $cfg{individual_results_folder}/%0" . $image_num_length . "d/puncta_highlight.png -sameq $cfg{exp_results_folder}/puncta_highlight.mov 2>&1",
	# "ffmpeg -v 0 -y -r $cfg{movie_frame_rate} -i $cfg{individual_results_folder}/%0" . $image_num_length . "d/degradation_highlight.png -sameq $cfg{exp_results_folder}/degradation_highlight.mov 2>&1",
	"ffmpeg -v 0 -y -r $cfg{movie_frame_rate} -i $cfg{individual_results_folder}/%0" . $image_num_length . "d/invado_and_not.png -sameq $cfg{exp_results_folder}/invado_and_not.mov 2>&1"
);

for (@commands) {
	if ($opt{debug}) {
		print $_, "\n";
	} else {
		print "\n\nBuild Movies\n\n" if $opt{debug};
		$t1 = new Benchmark;
		system "$_ > /dev/null 2> /dev/null";
		$t2 = new Benchmark;
		print "Runtime: ",timestr(timediff($t2,$t1)), "\n";
	}
}

###############################################################################
# Functions
###############################################################################
