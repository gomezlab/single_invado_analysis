#!/usr/bin/env perl

################################################################################
# Global Variables and Modules
################################################################################

use lib "../lib";

use strict;
use File::Path;
use File::Spec::Functions;
use File::Basename;
use File::Copy;
use Math::Matlab::Local;
use Getopt::Long;
use Data::Dumper;

use Config::ImageSet qw(ParseConfig);

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|c=s", "debug|d", "lsf|l")
  or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Gathering Config\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

################################################################################
# Main Program
################################################################################

mkpath($cfg{individual_results_folder});

my @image_sets = ([qw(puncta_image_folder puncta_image_file)],
				  [qw(gel_image_folder gel_image_file)],
				  [qw(ECM_image_folder ECM_image_file)]);

my @matlab_code;
my $all_images_empty = 1;

foreach (@image_sets) {
    my $folder   = $cfg{ $_->[0] };
    my $out_file = $cfg{ $_->[1] };
    
	next if (not(defined($folder)));

	if ($opt{debug}) {
		print "Searching: $cfg{exp_data_folder}/$folder/\n";
	}
    my @image_files = sort <$cfg{exp_data_folder}/$folder/*>;
    my @image_files = map { $_ =~ s/\'//g; $_; } @image_files;
	if (scalar(@image_files) == 1 &&
		&is_file_TIFF($image_files[0])) {
		&convert_from_tiff($image_files[0]);
    	@image_files = sort <$cfg{exp_data_folder}/$folder/*>;
	}

    $all_images_empty = 0 if (@image_files);

    if ($opt{debug}) {
        if (scalar(@image_files) > 1) {
            print "Image files found: $image_files[0] - $image_files[$#image_files]\n";
        } elsif (scalar(@image_files) == 0) {
            print "No image files found matching $cfg{exp_data_folder}/$folder, moving onto next image set.\n";
            next;
        } else {
            print "Image file found: $image_files[0]\n";
        }
        print "For Config Variable: ", $_->[0], "\n\n";
    } else {
        next if (not @image_files);
    }
	
	my $f_count = split("",@image_files);
	for (1..scalar(@image_files)) {
		my $source_file = $image_files[$_ - 1];
		my $output_folder_num = sprintf('%0'.$f_count.'d',$_);
		my $output_folder = "$cfg{individual_results_folder}/$output_folder_num";
		my $output_file = "$output_folder/$out_file";

		my $dir_command = "mkpath(dirname(\"$output_file\"))";
		my $copy_command = "copy(\"$source_file\",\"$output_file\")";
		if ($opt{debug}) {
			print "$dir_command\n$copy_command\n";
		} else {
			eval($dir_command);
			eval($copy_command);
		}
	}
}
die "Unable to find any images to include in the new experiment" if $all_images_empty;

my $flat_field_file = catfile($cfg{exp_data_folder},$cfg{flat_field_file});
if (-e $flat_field_file) {
	my $target_dir = catdir($cfg{results_folder}, $cfg{exp_name},'puncta_props');
	mkpath($target_dir);

	copy($flat_field_file,$target_dir);
}

################################################################################
#Functions
################################################################################

sub convert_from_tiff {
	my $tiff_file = $_[0];
	my $command = "convert $tiff_file " . catdir(dirname($tiff_file),"data_%05d.png") . ";";
	system($command);
	unlink $tiff_file;
}

sub is_file_TIFF {
     my $file = shift @_;
 
     my $type_output = `file $file`;
     if ($type_output =~ /TIFF/) {
         return 1;
     } else {
         return 0;
     }
}
