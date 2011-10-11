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
use Config::General qw(ParseConfig);

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
$opt{small_output} = 0;
$opt{movie} = 0;
GetOptions(\%opt, "src=s", "target=s", "debug|d", 
	"default_config|cfg|default_cfg=s","small_output","movie") or die;
die "Can't find src folder on command line." if (not $opt{src});
die "Can't find target folder on command line." if (not $opt{target});
die "Can't find default_config on command line." if (not $opt{default_config});

my $default_config_file_location = $opt{default_config};

my $relative_config_file;
if ($default_config_file_location =~ m#\.\./\.\./data/(.*)#) {
	$relative_config_file = $1;
} else {
	die "Unable to find ../../data/ in the default config file location: $default_config_file_location";
}

my @field_translation = (
	['FITC', 'puncta'],
	['TxRd', 'gel'],
);

###############################################################################
# Main Program
###############################################################################

my @src_files = <$opt{src}/*>;
die "Couldn't find any files in $opt{src}." if scalar(@src_files) == 0;
my %file_sets;

my $log_file;

#bundle each of the src files into the file_sets hash based on the position code
my %stage_position_sets;
for (@src_files) {
	if ($_ =~ /.*\.nd/) {
		$log_file = $_;
		next;
	}
	if ($_ =~ /_w\d+(.*)_s(\d+)_t(\d+).TIF/) {
		#which kind of image is this? texas red/GFP/etc?
		my $image_type = $1;

		#which stage position was this image taken at?
		my $stage_position = $2;
	
		#when was the image taken
		my $image_time = $3;
			
		$stage_position_sets{$stage_position}{$image_type}[$image_time] = $_;
	} else {
		warn "Unrecoginable filename ($_), can't find position code.";
	}
}

&confirm_stage_hash(%stage_position_sets);

print "Done with bundling the source files into sets\n\n";

#now we sort through each stage position to move the files into place and create the config files
my $current_max_exp_num = &determine_max_exp_num($opt{target});
my $data_folder = &determine_data_folder_name($opt{default_config});
foreach my $position (sort {$a<=>$b} keys %stage_position_sets) {
	my $position_directory = catdir($opt{target}, "time_series_" . 
		sprintf('%02d',$position + $current_max_exp_num));
	
	my $exp_name;
	if ($position_directory =~ /$data_folder(.*)/) {
		$exp_name = $1;
	} else {
		die "Unable to find exp name in $position_directory";
	}
	
	if ($opt{debug}) {
		print "mkpath($position_directory);\n";
		print "copy($log_file, catfile($position_directory, 'log.txt'));\n";
		print "&write_standard_config_file(catfile($position_directory,'analysis.cfg'),$exp_name)\n";
	} else {
		mkpath($position_directory);
		copy($log_file, catfile($position_directory, 'log.txt'));
		&write_standard_config_file(catfile($position_directory,'analysis.cfg'),$exp_name);
	}

	foreach my $field (keys %{$stage_position_sets{$position}}) {
		my $image_type = &find_image_type($field,\@field_translation);
		
		my $no_spaces_field = $field;
		$no_spaces_field =~ s/ /_/g;

		my @temp = split("",@{$stage_position_sets{$position}{$field}});
		my $image_num_length = scalar(@temp);
		
		my $image_directory = catdir($position_directory, 'Images', $image_type);

		#we start the scan through each position field combination at one
		#because it will always be empty as confirmed by the earlier hash check
		foreach my $image_num (1 .. $#{$stage_position_sets{$position}{$field}}) {
			my $source_file = ${$stage_position_sets{$position}{$field}}[$image_num];
			
			my $suffix;
			if ($source_file =~ /.*\.(.*)/) {
				$suffix = $1;
			} else {
				die "Unable to find file suffix in: $source_file";
			}

			my $image_filename = $no_spaces_field . "_" . 
				sprintf('%0' . $image_num_length . 'd', $image_num) . 
				".$suffix";
			
			my $target_file = catfile($image_directory, $image_filename);
			
			if ($image_num == 1) {
				if ($opt{debug}) {
					print "mkpath(dirname($target_file));\n";
				} else {
					mkpath(dirname($target_file));
				}
			}
			
			$target_file =~ s/tif/png/i;
			
			my $resize_str = '';
			if ($opt{small_output}) {
				$resize_str = "-resize 50% -normalize ";
			}
			my $command = "convert \"$source_file\" $resize_str\"$target_file\"";
			if ($opt{debug}) {
				#print "$command\n\n";
			} else {
				system "$command 2>/dev/null";
			}
		}

		if ($opt{movie}) {
			my $movie_source_dir = catdir($image_directory,"$no_spaces_field" ."_%03d.png");
			my $output_file = catdir($position_directory,"$no_spaces_field.mp4");
			my $command = "ffmpeg -i $movie_source_dir $output_file > /dev/null";
			system($command);
		}
		# print catdir($image_directory,"$no_spaces_field" ."_%03d.png"),"\n";
		# print catdir($position_directory,"$no_spaces_field.mp4"),"\n";
	}

}

###############################################################################
# Functions
###############################################################################

sub find_image_type {
	my $field_name = $_[0];
	my @field_translation = @{$_[1]};
	
	my $image_type;
	foreach my $trans_set_ref (@field_translation) {
		if ($field_name =~ /${$trans_set_ref}[0]/) {
			die "Already found field for $field_name -> $image_type" if defined $image_type;
			$image_type = ${$trans_set_ref}[1];
		}
	}
	die "Didn't find any appropriate translation for $field_name" if not defined $image_type;
	
	return $image_type;
}

sub confirm_stage_hash {
	my %stage_hash = @_;
	
	#make sure the number of fields is the same in each position
	my $number_of_fields;
	foreach my $position (sort {$a<=>$b} keys %stage_hash) {
		if (defined $number_of_fields) {
			die "Different number of fields in stage position $position." 
				if $number_of_fields != scalar(keys %{$stage_hash{$position}});
		} else {
			$number_of_fields = scalar(keys %{$stage_hash{$position}});
		}
	}
	
	#make sure the number of images in each field in each position is the same
	foreach my $position (sort {$a<=>$b} keys %stage_hash) {
		my $image_num;
		foreach my $field (keys %{$stage_hash{$position}}) {
			if (defined $image_num) {
				die "Different number of images in each field in stage position $position." 
				if $image_num != scalar(@{$stage_hash{$position}{$field}});
			} else {
				$image_num = scalar(@{$stage_hash{$position}{$field}});
			}

			if (defined ${$stage_hash{$position}{$field}}[0]) {
				die "Expected the 0 index to be undefined, instead: ${$stage_hash{$position}{$field}}[0]";
			}
		}
	}
}

sub write_standard_config_file {
	my $file_name = shift;
	my $exp_name = shift;

	my $default_config_depth = &find_main_config_depth($file_name);
	my $config_prefix = "../" x $default_config_depth;
	

	open OUTPUT, ">$file_name";
	print OUTPUT 
"<<include " . $config_prefix . $relative_config_file . ">>

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

sub determine_data_folder_name {
	my $config_file = $_[0];

	my %cfg = ParseConfig(
        -ConfigFile            => $config_file,
        -MergeDuplicateOptions => 1,
        -IncludeRelative       => 1,
	);
	
	return $cfg{data_folder};
}

sub determine_max_exp_num {
	my $target_dir = shift;

	my @files = <$target_dir/*>;

	my $max_num = 0;
	foreach (@files) {
		if ($_ =~ /(\d+)$/) {
			if ($1 > $max_num) {
				$max_num = $1;
			}
		}
	}
	
	return $max_num;
}

sub find_main_config_depth {
	my $file_name = shift;

	my @split_file_parts = split(/\//,$file_name);
	
	while ($split_file_parts[0] eq "..") {
		@split_file_parts = @split_file_parts[1..$#split_file_parts];
	}
	
	return scalar(@split_file_parts) - 2;
}
