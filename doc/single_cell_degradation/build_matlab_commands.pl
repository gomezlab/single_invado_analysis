#!/usr/bin/perl -w

use Data::Dumper;

#Get input file
my $gold_standard_file = 'gold_standards/Matt.txt';
open INPUT, $gold_standard_file or die "$!";
my @standard = <INPUT>;
close INPUT;

#build matlab commands
my %prop_location_lookup = (
"2011_12_02-BB-time_series_02-028" => '2011_12_02/BB94/time_series_02/cell_props/lin_time_series/',
"2011_12_02-co-time_series_01-001" => '2011_12_02/control/time_series_01/cell_props/lin_time_series/',
"2011_12_02-co-time_series_01-002" => '2011_12_02/control/time_series_01/cell_props/lin_time_series/',
"2011_12_02-co-time_series_01-003" => '2011_12_02/control/time_series_01/cell_props/lin_time_series/',
"2011_12_05-BB-time_series_03-003" => '2011_12_05/BB94/time_series_03/cell_props/lin_time_series/',
"2011_12_05-BB-time_series_04-006" => '2011_12_05/BB94/time_series_04/cell_props/lin_time_series/',
"2011_12_05-co-time_series_01-005" => '2011_12_05/control/time_series_01/cell_props/lin_time_series/',
"2011_12_05-co-time_series_02-002" => '2011_12_05/control/time_series_02/cell_props/lin_time_series/',
"2011_12_05-co-time_series_02-004" => '2011_12_05/control/time_series_02/cell_props/lin_time_series/',
"2011_12_05-co-time_series_02-012" => '2011_12_05/control/time_series_02/cell_props/lin_time_series/',
"2011_12_05-co-time_series_03-001" => '2011_12_05/control/time_series_03/cell_props/lin_time_series/',
"2011_12_05-co-time_series_03-002" => '2011_12_05/control/time_series_03/cell_props/lin_time_series/',
"2011_12_05-co-time_series_03-004" => '2011_12_05/control/time_series_03/cell_props/lin_time_series/',
"2011_12_05-co-time_series_05-001" => '2011_12_05/control/time_series_05/cell_props/lin_time_series/',
"2011_12_05-co-time_series_05-002" => '2011_12_05/control/time_series_05/cell_props/lin_time_series/',
"2011_12_16-BB-time_series_01-004" => '2011_12_16/BB94/time_series_01/cell_props/lin_time_series/',
"2011_12_16-BB-time_series_02-004" => '2011_12_16/BB94/time_series_02/cell_props/lin_time_series/',
"2011_12_16-BB-time_series_02-005" => '2011_12_16/BB94/time_series_02/cell_props/lin_time_series/',
"2011_12_16-BB-time_series_02-006" => '2011_12_16/BB94/time_series_02/cell_props/lin_time_series/',
"2011_12_16-BB-time_series_03-001" => '2011_12_16/BB94/time_series_03/cell_props/lin_time_series/',
"2011_12_16-BB-time_series_03-005" => '2011_12_16/BB94/time_series_03/cell_props/lin_time_series/',
"2011_12_16-BB-time_series_04-002" => '2011_12_16/BB94/time_series_04/cell_props/lin_time_series/',
"2011_12_16-BB-time_series_04-004" => '2011_12_16/BB94/time_series_04/cell_props/lin_time_series/',
"2011_12_16-BB-time_series_04-005" => '2011_12_16/BB94/time_series_04/cell_props/lin_time_series/',
"2011_12_16-BB-time_series_05-001" => '2011_12_16/BB94/time_series_05/cell_props/lin_time_series/',
"2011_12_16-co-time_series_01-001" => '2011_12_16/control/time_series_01/cell_props/lin_time_series/',
"2011_12_16-co-time_series_02-001" => '2011_12_16/control/time_series_02/cell_props/lin_time_series/',
"2011_12_16-co-time_series_02-005" => '2011_12_16/control/time_series_02/cell_props/lin_time_series/',
"2011_12_16-co-time_series_02-006" => '2011_12_16/control/time_series_02/cell_props/lin_time_series/',
"2011_12_16-co-time_series_03-001" => '2011_12_16/control/time_series_03/cell_props/lin_time_series/',
"2011_12_16-co-time_series_03-002" => '2011_12_16/control/time_series_03/cell_props/lin_time_series/',
"2011_12_16-co-time_series_03-003" => '2011_12_16/control/time_series_03/cell_props/lin_time_series/',
"2011_12_16-co-time_series_03-005" => '2011_12_16/control/time_series_03/cell_props/lin_time_series/',
"2011_12_16-co-time_series_03-006" => '2011_12_16/control/time_series_03/cell_props/lin_time_series/',
"2011_12_17-BB-time_series_01-001" => '2011_12_17/BB94/time_series_01/cell_props/lin_time_series/',
"2011_12_17-BB-time_series_01-002" => '2011_12_17/BB94/time_series_01/cell_props/lin_time_series/',
"2011_12_17-BB-time_series_02-002" => '2011_12_17/BB94/time_series_02/cell_props/lin_time_series/',
"2011_12_17-BB-time_series_02-006" => '2011_12_17/BB94/time_series_02/cell_props/lin_time_series/',
"2011_12_17-BB-time_series_03-013" => '2011_12_17/BB94/time_series_03/cell_props/lin_time_series/',
"2011_12_17-BB-time_series_03-017" => '2011_12_17/BB94/time_series_03/cell_props/lin_time_series/',
"2011_12_17-Co-time_series_01-002" => '2011_12_17/Control/time_series_01/cell_props/lin_time_series/',
"2011_12_17-Co-time_series_01-008" => '2011_12_17/Control/time_series_01/cell_props/lin_time_series/',
"2011_12_17-Co-time_series_02-002" => '2011_12_17/Control/time_series_02/cell_props/lin_time_series/',
"2011_12_17-Co-time_series_02-012" => '2011_12_17/Control/time_series_02/cell_props/lin_time_series/',
);


for (@standard) {
	chomp;
	if (/(.*): (.*)/) {
		my $file_name = $1;
		my $degrade_nums = $2;
		my $command_start = "process_results('../../results/Invado_count_20000/'";
		$file_name =~ s/\..*//;
		
		my $cell_num;
		if ($file_name =~ /(\d+)$/) {
			$cell_num = $1;
		} else {
			die "Couldn't find cell num: $file_name";
		}

		my $prop_location = $prop_location_lookup{$file_name};
		if (not defined $prop_location) {
			die $file_name;
		}

		my $command = "$command_start,'$prop_location',[$degrade_nums],$cell_num,'$file_name')";

		print "$command\n";
	} else {
		die "messed up line: $_";
	}
}
