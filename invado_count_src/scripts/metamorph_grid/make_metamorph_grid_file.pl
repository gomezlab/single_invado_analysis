#!/usr/bin/perl -w

################################################################################
# Global Variables and Modules
################################################################################
use strict;
use Data::Dumper;
use Getopt::Long;
use Statistics::Regression;

my %opt;
$opt{dish_count} = 8;
$opt{debug} = 0;
GetOptions(\%opt,"corners=s","dish_count=s","debug|d") or die;

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
	my %corner_positions = &read_in_corner_positions($opt{corners});
	my %corner_models = &build_field_models(%corner_positions);
	
	@stage_lines = &extrapolate_stage_positions(\$opt{corners}, \%corner_models);
	
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

sub read_in_corner_positions {
	my $file = shift;
	
	my %corner_hash;
	
	open INPUT, "$file";
	#burn the first four lines, they are headers
	for (1..4) {
		my $temp = <INPUT>;
	}
	while (<INPUT>) {
		my @line_split = split(",", $_);
		if ($line_split[3] != 1111) {
			push @{$corner_hash{$line_split[7]}{z}}, $line_split[3]; 
			push @{$corner_hash{$line_split[7]}{x}}, $line_split[1]; 
			push @{$corner_hash{$line_split[7]}{y}}, $line_split[2]; 
		}
	}
	close INPUT;
	
	return %corner_hash;
}

sub build_field_models {
	my %corner_positions = @_;

	my %field_models;

	foreach (keys %corner_positions) {
		my %corner_set = %{$corner_positions{$_}};
		
		my @z = @{$corner_set{z}};
		my @x = @{$corner_set{x}};
		my @y = @{$corner_set{y}};
		
		my $reg = Statistics::Regression->new( "corner regression", ["const", "X_POS", "Y_POS" ] );
		
		for my $i (0..$#z) {
			$reg->include($z[$i], [1.0, $x[$i], $y[$i]]);
		}
	
		$field_models{$_} = $reg;
	}
	
	return %field_models;
}

sub extrapolate_stage_positions {
	my $corner_file = ${$_[0]};
	my %corner_models = %{$_[1]};

	my @final_positions;
	
	open INPUT, "$corner_file";
	#burn the first four lines, they are headers
	for (1..4) {
		my $temp = <INPUT>;
	}
	while (<INPUT>) {
		chomp($_);
		my @line_split = split(",", $_);
		if ($line_split[3] != 1111) {
			push @final_positions, $_; 
		} else {
			my $stage_pos = $line_split[7];

			my @theta = $corner_models{$stage_pos}->theta();

			my $new_z = $theta[0] + $theta[1]*$line_split[1] + $theta[2]*$line_split[2];
			
			$line_split[3] = $new_z;
			$line_split[5] = $new_z;
			
			push @final_positions, join(",", @line_split);
		}
	}
	close INPUT;
	
	return @final_positions;
}
