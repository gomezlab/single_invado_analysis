#!/usr/bin/perl

use threads;
use threads::shared;
use File::Spec::Functions;
use File::Basename;

$| = 1;

my $max_processes = 4;

my @config_files = <../../data/time_series_*/*.cfg>;

my @data_folders = map catfile(dirname($_),"run.txt"), @config_files;
my @time_series_list = map { /(time_series_\d*)/; $1; } @config_files;

my @processes :shared;

@processes = map { "nice -20 ./build_results.pl -cfg $config_files[$_] -d > $data_folders[$_]" } (0..$#config_files);

my @started;
while (@processes) {
	while (@processes && &gather_running_status(@started) < $max_processes) {
		my $command = pop @processes;
		my $config = pop @config_files;
		push @started, $config;

		print "Executing $command\n";
		threads->create('execute_process',$command)->detach;
	}
}

while(&gather_running_status(@started)) {
}

sub print_sys_line {
	my $process_string = shift;
	sleep rand(10)+5;
	print "From Func: $process_string\n";
}

sub execute_process {
	my $process_string = shift;
	system($process_string);
	print "Finished $process_string\n";
}

sub gather_running_status {
	my @started_configs = @_;
	
	my $running = 0;
	foreach (@started_configs) {
		my $status_file = catfile(dirname($_),'status.txt');
		if (-e $status_file) {
			open STATUS, $status_file;
			my $line = <STATUS>;
			close STATUS;

			chomp($line);
			if (not($line =~ /DONE/)) {
				$running++;
			}
		} else {
			$running++;
		}
	}
	return $running;
}
