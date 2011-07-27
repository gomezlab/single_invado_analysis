#!/usr/bin/perl -w

use strict;
use File::Path;
use File::Copy;
use Data::Dumper;
use Storable;

my %exp_nums;
my @exp_folders = <../../results/Invadopodia_testing/control/*/>;
my $count = 1;
open OUTPUT, ">exp_nums.txt" or die;
for (@exp_folders) {
	if ($_ =~ /\/(\d.*?pos.*?)\//) {
		print OUTPUT "$1 => $count\n";
		$exp_nums{$1} = $count;
		$count++;
	} 
}
close OUTPUT;

store \%exp_nums, 'exp_nums.stor'; 

my @invado_sm_folders = <../../results/Invadopodia_testing/control/*/visualizations/single_invadopodia/no_highlights/>;
my @notinvado_sm_folders = <../../results/Invadopodia_testing/control/*/visualizations/single_notinvadopodia/no_highlights/>;

mkpath('sm_sets');

#invadopodia processing
my %file_sets_invado;
for (@invado_sm_folders) {
	if ($_ =~ /\/(\d.*?pos.*?)\//) {
		my $pos_code = $1;
		my @files = <$_/*>;
		for my $file (@files) {
			if ($file =~ /\/(\d+).png/) {
				my $output_name = "sm_sets/$1_$exp_nums{$pos_code}.png";
				$file_sets_invado{$file} = $output_name;
			}
		}
	}
}
my @keys = keys %file_sets_invado;
print "Found $#keys\n";
for (0..99) {
	my $rand_key_num = sprintf("%.0f", rand($#keys));
	my $file = $keys[$rand_key_num];
	# print "$file\n$file_sets_invado{$file}\n\n";
	copy($file,$file_sets_invado{$file}) or die "$!";
	
	delete $file_sets_invado{$file};
	@keys = keys %file_sets_invado;
}

#not invadopodia processing
my %file_sets_notinvado;
for (@notinvado_sm_folders) {
	if ($_ =~ /\/(\d.*?pos.*?)\//) {
		my $pos_code = $1;
		my @files = <$_/*>;
		for my $file (@files) {
			if ($file =~ /\/(\d+).png/) {
				my $output_name = "sm_sets/$1_$exp_nums{$pos_code}.png";
				$file_sets_notinvado{$file} = $output_name;
			}
		}
	}
}

@keys = keys %file_sets_notinvado;
print "Found $#keys\n";
for (0..99) {
	my $rand_key_num = sprintf("%.0f", rand($#keys));
	my $file = $keys[$rand_key_num];
	# print "$file\n$file_sets_notinvado{$file}\n\n";
	copy($file,$file_sets_notinvado{$file}) or die "$!";
	
	delete $file_sets_notinvado{$file};
	@keys = keys %file_sets_notinvado;
}
