#!/usr/bin/perl -w

################################################################################
# Global Variables and Modules
################################################################################
use strict;
use Data::Dumper;
use Getopt::Long;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "corners=s") or die;

# "Stage Memory List", Version 5.0
# 0, 0, 0, 0, 0, 0, 0, "um", "um"
# 0
# 7
# "POSITION_LABLE", 32.7, -38.2, 2500, 0, 2500, FALSE, 1, TRUE, TRUE, 0
# "POSITION_LABEL", 32.7, 29.9, 2500, 0, 2500, FALSE, 1, TRUE, TRUE, 0
# "POSITION_LABEL_1", -10, 29.9, 2500, 0, 2500, FALSE, 1, TRUE, TRUE, 0
# "POSITION_LABEL_2", 10, -0.2, 2500, 0, 2500, FALSE, 2, TRUE, TRUE, 0
# "POSITION_LABEL_3", 10, -20.3, 2500, 0, 2500, FALSE, 2, TRUE, TRUE, 0
# "POSITION_LABEL_4", 10, -0.2, 2500, 0, 2500, FALSE, 3, TRUE, TRUE, 0
# "POSITION_LABEL_5", -0, -0.2, 2500, 0, 2500, FALSE, 4, TRUE, TRUE, 0

my $start_x = -400;
my $start_y = -400;
my $grid_size = 5;
my $move_increment = 200;

my $standard_z = 3000;

my $standard_label = "POS_";

my $dish_count = 8;

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

	for my $this_dish (1..$dish_count) {
		for my $x_increment (0..($grid_size-1)) {
			my $this_x = $x_increment*$move_increment+$start_x;
			for my $y_increment (0..($grid_size-1)) {
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

# sub build_complete_stage_list {
# 	my $position_count = 0;
# 	my @lines;
# 
# 	for my $this_dish (1..$dish_count) {
# 		for my $x_increment (0..($grid_size-1)) {
# 			my $this_x = $x_increment*$move_increment+$start_x;
# 			for my $y_increment (0..($grid_size-1)) {
# 				my $this_y = $y_increment*$move_increment+$start_y;
# 				$position_count++;
# 				push @lines, "\"$standard_label$position_count\", $this_x, $this_y, $standard_z, 0, $standard_z, FALSE, $this_dish, TRUE, TRUE, 0";
# 			}
# 		}
# 	}
# 	
# 	return @lines;
# }
# 
# sub read_in_corner_values {
# 	my %corner_focus;
# 	
# 	open INPUT, "$opt{corners}";
# 	#the first four lines contain header data, skip it
# 	for (1..4) {
# 		my $temp = <INPUT>;
# 	}
# 	my @all_lines = <INPUT>;
# 	close INPUT;
# 	
# 	my @xyz_dish;
# 	
# 	# my $current_line = 0;
# 	# my @corner;
# 	# for (1..$dish_count) {
# 	# 	@corner = split(",", $all_lines[$current_line]);
# 	# 	$current_line += ($grid_size - 2);
# 	# 	@corner = split(",", $all_lines[$current_line]);
# 	# 	$current_line += $grid_size*($grid_size - 2);
# 	# 	@corner = split(",", $all_lines[$current_line]);
# 	# 	$current_line += ($grid_size - 2);
# 	# 	@corner_3 = split(",", $all_lines[$current_line]);
# 	# }
# 
# 	# while (<INPUT>) {
# 	# 	my @split_line = split(",",$_);
# 	# 	#have to round the x and y because the scope sometimes decides to land
# 	# 	#on a slightly +-0.1 micron different spot that specified, shouldn't
# 	# 	#make a very big difference in the calculations
# 	# 	my $x = round($split_line[1]);
# 	# 	my $y = round($split_line[2]);
# 	# 	push @xyz_dish, ([$x,$y,$split_line[3],$split_line[7]]);
# 	# }
# 
# 	return @xyz_dish;
# }
