#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Spec::Functions;
use File::Basename;
use File::Path;
use File::Copy;
use Getopt::Long;
use Cwd;
use Data::Dumper;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "src=s", "target=s", "copy_mode", "debug|d") or die;
die "Can't find src folder on command line." if (not $opt{src});
die "Can't find target folder on command line." if (not $opt{target});

###############################################################################
# Main Program
###############################################################################

my @src_files = <$opt{src}/*>;
die "Couldn't find any files in $opt{src}." if scalar(@src_files) == 0;
my %file_sets;

#bundle each of the src files into the file_sets hash based on the position code
for (@src_files) {
	if ($_ =~ /Position (\d+)_XY/) {
		my $pos_str = "pos" . sprintf("%02d", $1);
		push @{$file_sets{$pos_str}}, $_;
	} else {
		die "Unrecoginable filename ($_), can't find position code.";
	}
}

#check that there are 3 entries in each file_sets hash array
for (keys %file_sets) {
	die "Wrong number of files in position code $_\n" . Dumper($file_sets{$_}) 
	  if (scalar(@{$file_sets{$_}}) != 3);
}

#find the date for each position code
my %pos_code_to_date;
for (keys %file_sets) {
	my @log_file = grep /.*.log/, @{$file_sets{$_}};
	die "Found more than one log file in set, postion code: $_" if (scalar(@log_file) > 1);
	die "Didn't find any log files, postion code: $_" if (scalar(@log_file) < 1);

	$pos_code_to_date{$_} = &get_data_from_log_file($log_file[0]);
}

#sort through each of the file sets and move everything into place
for (keys %file_sets) {
	my $exp_dir = catdir($opt{target}, $pos_code_to_date{$_} . "_" . $_);
	
	my $config_file = catfile($exp_dir, "analysis.cfg");
	my $exp_name = catdir(basename($opt{target}), basename($exp_dir));

	my $puncta_file = catfile($exp_dir, "Images", "puncta", "puncta.tiff");
	mkpath(dirname($puncta_file));
	my @src_puncta_file = grep /C0.*.tiff/, @{$file_sets{$_}};
	die "Found more than one puncta file in set, postion code: $_" if (scalar(@src_puncta_file) > 1);
	die "Didn't find any puncta files, postion code: $_" if (scalar(@src_puncta_file) < 1);

	my $gel_file = catfile($exp_dir, "Images", "gel", "gel.tiff");
	mkpath(dirname($gel_file));
	my @src_gel_file = grep /C1.*.tiff/, @{$file_sets{$_}};
	die "Found more than one gel file in set, postion code: $_" if (scalar(@src_gel_file) > 1);
	die "Didn't find any gel files, postion code: $_" if (scalar(@src_gel_file) < 1);
	
	my $log_file = catfile($exp_dir, "collection_properties.log");
	mkpath(dirname($log_file));
	my @src_log_file = grep /.*.log/, @{$file_sets{$_}};
	die "Found more than one log file in set, postion code: $_" if (scalar(@src_log_file) > 1);
	die "Didn't find any log files, postion code: $_" if (scalar(@src_log_file) < 1);
	
	if ($opt{debug}) {
		print "mkpath($exp_dir)\n";
		if ($opt{copy_mode}) {
			print "copy($src_gel_file[0], $gel_file)\n";
			print "copy($src_puncta_file[0], $puncta_file)\n";
			print "copy($src_log_file[0], $log_file)\n";
		} else {
			print "move($src_gel_file[0], $gel_file)\n";
			print "move($src_puncta_file[0], $puncta_file)\n";
			print "move($src_log_file[0], $log_file)\n";
		}
		print "&write_standard_config_file($config_file, $exp_name)\n";
		print "\n";
	} else {
		mkpath($exp_dir);
		if ($opt{copy_mode}) {
			copy($src_gel_file[0], $gel_file);
			copy($src_puncta_file[0], $puncta_file);
			copy($src_log_file[0], $log_file);
		} else {
			move($src_gel_file[0], $gel_file);
			move($src_puncta_file[0], $puncta_file);
			move($src_log_file[0], $log_file);
		}
		&write_standard_config_file($config_file, $exp_name);
	}
}
###############################################################################
# Functions
###############################################################################

sub get_data_from_log_file {
	my $file_name = shift;
	open INPUT, $file_name;
	
	my @file_lines = <INPUT>;

	close INPUT;

	my $date;

	if ($file_lines[1] =~ /Capture Date-Time: (\d+)\/(\d+)\//) {
		$date = sprintf("%02d",$1) . "_" . sprintf("%02d", $2);
	} else {
		die "Could not file capture date in file: $file_name";
	}
	
	return $date;
}

sub write_standard_config_file {
	my $file_name = shift;
	my $exp_name = shift;

	open OUTPUT, ">$file_name";
	print OUTPUT 
"<<include ../../../config/Invado_default.cfg>>

###############################################################################
# Required Experiment Parameters
###############################################################################

exp_name = $exp_name 

###############################################################################
# Experiment Specific Parameters
###############################################################################
";

	close OUTPUT;


	
	
}
