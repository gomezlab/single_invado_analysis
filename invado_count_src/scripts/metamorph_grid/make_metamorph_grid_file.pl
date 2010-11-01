#!/usr/bin/perl -w

################################################################################
# Global Variables and Modules
################################################################################
use strict;
use Data::Dumper;
use Getopt::Long;

my %opt;
$opt{dish_count} = 8;
$opt{debug} = 0;
GetOptions(\%opt, "corners=s", "dish_count=s") or die;

my $start_x = -1600;
my $start_y = -1600;
my $grid_size = 5;
my $move_increment = 800;

my $standard_z = 3000;

my $standard_label = "POS_";


my @stage_header;
push @stage_header, "\"Stage Memory List\", Version 5.0";
push @stage_header, "0, 0, 0, 0, 0, 0, 0, \"um\", \"um\"";
push @stage_header, "0";
push @stage_header, "NUMBER_OF_POS";

################################################################################
# Main
################################################################################

my @stage_lines;
if (not defined $opt{corners}) {
	@stage_lines = &build_corners_input_positions;
	&output_stage_lines("stage_corners.STG");
} else {
	system("R CMD BATCH --vanilla extrapolate_focus.R /dev/null");
	open INPUT, "stage_pos_predicted.csv";
	while (<INPUT>) {
		chomp($_);
		push @stage_lines, $_;
	}
	close INPUT;
	&output_stage_lines("final_stage_positions.STG");
}

################################################################################
# Functions
################################################################################
sub output_stage_lines {
	$stage_header[3] = scalar(@stage_lines);

	unshift @stage_lines, @stage_header;
	
	open OUTPUT, ">$_[0]";
	print OUTPUT join("\n",@stage_lines);
	close OUTPUT;
}

sub build_corners_input_positions {
	my $position_count = 0;
	my @lines;

	for my $this_dish (1..$opt{dish_count}) {
		for my $x_increment (0..($grid_size-1)) {
			my $this_x = $x_increment*$move_increment+$start_x;
			
			#to build a grid that goes in a zig zag, reverse the y increment sequence every other time
			my @y_increment_seq = (0..($grid_size-1));
			if ($x_increment % 2 == 0) {
				@y_increment_seq = reverse(@y_increment_seq);
			}

			for my $y_increment (@y_increment_seq) {
				my $this_y = $y_increment*$move_increment+$start_y;
				$position_count++;
				
				my $position_z;
				
				#determine if we are at a corner, if so, we put in the standard
				#z because we'd like to start close to the real focus,
				#otherwise, put a dummy value in
				if (($x_increment == 0 && $y_increment == 0) ||
					($x_increment == ($grid_size - 1) && $y_increment == 0) ||
					($x_increment == 0 && $y_increment == ($grid_size - 1)) ||
					($x_increment == ($grid_size - 1) && $y_increment == ($grid_size - 1))
					) {
					$position_z = $standard_z;
				} else {
					$position_z = 1111;
				}
				push @lines, "\"$standard_label$position_count\",$this_x,$this_y,$position_z,0,$position_z,FALSE,$this_dish,TRUE,TRUE,0";
			}
		}
	}
	
	return @lines;
}
