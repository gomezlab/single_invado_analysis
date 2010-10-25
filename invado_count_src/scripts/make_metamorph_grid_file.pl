#!/usr/bin/perl -w

use strict;
use Data::Dumper;
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
my $x_grid_size = 5;
my $y_grid_size = 5;
my $move_increment = 200;

my $standard_z = 2700;

my $standard_label = "POS_";

my $dish_count = 8;

my @stage_lines;
push @stage_lines, "\"Stage Memory List\", Version 5.0";
push @stage_lines, "0, 0, 0, 0, 0, 0, 0, \"um\", \"um\"";
push @stage_lines, "0";
push @stage_lines, "NUMBER_OF_POS";

my $position_count = 0;

for my $this_dish (1..$dish_count) {
	for my $x_increment (0..($x_grid_size-1)) {
		my $this_x = $x_increment*$move_increment+$start_x;
		for my $y_increment (0..($y_grid_size-1)) {
			my $this_y = $y_increment*$move_increment+$start_y;
			$position_count++;
			push @stage_lines, "\"$standard_label$position_count\", $this_x, $this_y, $standard_z,0,$standard_z, FALSE, $this_dish, TRUE, TRUE, 0";
		}
	}
}

$stage_lines[3] = $position_count;

open OUTPUT, ">stage_mem_list.STG";
print OUTPUT join("\n",@stage_lines);
close OUTPUT;
