#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Temp qw/ tempfile tempdir /;
use File::Path;
use Config::General qw/ ParseConfig /;
use Getopt::Long;
use Data::Dumper;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg=s","debug|d");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my %cfg = &get_config;


###############################################################################
# Main Program
###############################################################################

my @folders = <$cfg{data_folder}/*/raw_data>;

if ($opt{debug}) {
	print "Example Data Folder: ",join("\n",@folders[0..0]),"\n";
	@folders = $folders[0];
}


for (@folders) {
	$_ =~ /($cfg{data_folder}\/(\d+)\/)/;
	my $picture_num = $2;
	
	my $plots_folder = "$1/plots";
	mkpath($plots_folder);

	my $file_name = &write_data_file("$_/Centroid_dist_from_edge","$_/Area");
	foreach (split(/\s/,$cfg{file_ext})) {
		&build_and_execute_gnuplot_file("cent_area",$file_name,"$plots_folder/Cent_dist_vs_area.$_");
	}

	$file_name = &write_data_file("$_/Centroid_dist_from_edge", "$_/Average_adhesion_signal", "$_/Variance_adhesion_signal");
	foreach (split(/\s/,$cfg{file_ext})) {
		&build_and_execute_gnuplot_file("cent_sig",$file_name,"$plots_folder/Cent_dist_vs_sig.$_");
	}

	$file_name = &write_data_file("$_/Area");
	foreach (split(/\s/,$cfg{file_ext})) {
		&build_and_execute_gnuplot_file("area_hist",$file_name,"$plots_folder/Area_hist.$_");
	}

}


###############################################################################
# Functions
###############################################################################

sub get_config {
	if ($opt{debug}) {
		print "Collecting Configuration\n";
	}

	my %default_config = (file_ext => "svg png");
	my %cfg = ParseConfig(-ConfigFile => $opt{cfg},
						  -DefaultConfig => \%default_config,
						  -MergeDuplicateOptions => 1,
	);
	if (defined($cfg{pixel_size_file}) && defined($cfg{target_unit_size})) {
		$cfg{scale_values} = 1;
		if (-e $cfg{pixel_size_file}) {
			open INPUT, "$cfg{pixel_size_file}" or die "Can't open specified pixel_size_file: $cfg{pixel_size_file}.";
			$cfg{pixel_size} = <INPUT>;
			chomp($cfg{pixel_size});
			close INPUT;
		} else {
			die "Can't find specified pixel_size_file: $cfg{pixel_size_file}.";
		}
		if($opt{debug}) {
			print "Pixel Size -> $cfg{pixel_size}, Target Size -> $cfg{target_unit_size}.\n";
		}
	}

	return %cfg;
}

sub gather_data_from_matlab_file {
	my ($file) = @_;
	
	open INPUT, "$file" or die "Problem opening $file"; 
	my @in = split("   ",<INPUT>);
	close INPUT;

	if(defined($cfg{scale_values})) {
		my $scale_factor = $cfg{pixel_size}/$cfg{target_unit_size};
		if ($file =~ /area/i) {
			$scale_factor = ($cfg{pixel_size}/$cfg{target_unit_size})**2;
		}

		for (1..$#in) {
			if ($file =~ /adhesion_signal/) {
				next;
			}
			$in[$_] = $in[$_] *$scale_factor;
		}
	}

	return @in;
}

sub build_data_file_from_matlab {
	my @data;

	foreach (@_) {
		my @temp = &gather_data_from_matlab_file($_);
		push @data, \@temp;
	}
	
	if ($opt{debug}) {
		my $data_size = scalar(@{$data[0]});
		foreach (0..$#data) {
			if (scalar(@{$data[$_]} != $data_size)) {
				print "Problem with files: $_[0], $_[$_], the number of entries do ",
					  "not match, they are: ", scalar(@{$data[0]}), " and ", scalar(@{$data[$_]}), "\n";
			}
			if (${$data[$_]}[0] =~ /\s+/) {
				print "First entry in file ''$_[$_] is not empty, but is expected to be.\n";
			}
		}	
	}

	my @output;
	foreach my $i (1..scalar(@{$data[0]}) - 1) {
		my $temp = "";
		foreach my $set_num (0..$#data) {
			if ($set_num != $#data) {
				$temp = $temp . sprintf("%f",${$data[$set_num]}[$i]). "	";
			} else {
				$temp = $temp . sprintf("%f",${$data[$set_num]}[$i]);
			}
		}
		push @output, $temp;
	}
	return @output;
}

sub write_data_file {
	my @header_line;
	foreach (@_) {
		$_ =~ /.*\/(.*)/;
		push @header_line, $1;
	}
	
	my @data_lines = &build_data_file_from_matlab(@_);

	if ($opt{debug}) {
		open OUTPUT, ">$_/".join("_",@header_line)."_data_file" or die "Unable to open the output file";
		print OUTPUT "#", join("	",@header_line), "\n", join("\n",@data_lines);
		close OUTPUT;
	}

	my ($temp_h,$file_name) = tempfile();
	print $temp_h "#", join("	",@header_line), "\n", join("\n",@data_lines);
	close $temp_h;
	return $file_name;
}

sub build_and_execute_gnuplot_file {
	my @gnu_commands = &build_gnuplot_file(@_);
	
	my ($temp_h,$gnuplot_file_name) = tempfile();
	print $temp_h join("\n",@gnu_commands);
	close $temp_h;

	system "gnuplot $gnuplot_file_name";
}

sub build_gnuplot_file {
	my ($plot_type,$data_file_name,$out_file) = @_;
	
	my ($title,$xlabel,$ylabel) = &get_title_and_labels($plot_type);
	
	my $term_type = &get_output_type($out_file);
	
	my @commands;
	if ($plot_type =~ /hist/) {
		push @commands, ("set terminal $term_type","set key off","set style data histogram",
						 "set title \"$title\"","set xlabel \"$xlabel\"",
						 "set ylabel \"$ylabel\"","set output '$out_file'",
						 "plot '$data_file_name'");
	} else {
		push @commands, ("set terminal $term_type","set key off",
						 "set title \"$title\"","set xlabel \"$xlabel\"",
						 "set ylabel \"$ylabel\"","set output '$out_file'",
						 "plot '$data_file_name'");
	}
	return @commands;
}

sub get_title_and_labels {
	my $plot_type = $_[0];
	my %t_hash = (cent_area => ["Centroid Distance from Cell Edge vs. Area of Identified Adhesions",
							    "Centroid Distance from Cell Edge",
							    "Area of Identified Adhesions",],
				  cent_sig  => ["Centroid Distance from Cell Edge vs. Average Adhesion Signal",
							    "Centroid Distance from Cell Edge",
							    "Average Normalized Fluorescence Signal",],
			   	  area_hist => ["Area of Identified Adhesions",
				  				"",
								"",],
				 );
	
	if (not(defined($t_hash{$plot_type}))) {
		die "Plot type, $plot_type, unrecognized, possible values: ", join(", ", keys %t_hash), "\n";
	}

   	return @{$t_hash{$plot_type}};
}

sub get_output_type {
	my $out_file = $_[0];
	if (not($out_file =~ /.*\.(.*)/)) {
		die "Can't figure out the file extension for determining the gnuplot terminal type, file - $out_file\n";
	}
	my $file_type = $1;
	my %f_types_to_terms = ("png" => "png", "svg" => "svg enhanced", "jpg" => "jpeg", "jpeg" => "jpeg", "pdf" => "pdf enhanced", "eps" => "postscript eps enhanced");

	if (not(defined($f_types_to_terms{$file_type}))) {
		die "In file: $out_file\n", "Can't match determined file extension - $file_type - to a gnuplot ",
			"term name, possible file extensions: ", 
			join(", ", keys %f_types_to_terms), "\n";
	}

	return $f_types_to_terms{$file_type};
}
