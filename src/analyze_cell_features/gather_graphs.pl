#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Temp qw/ tempfile tempdir /;
use File::Path;
use Config::General;
use Getopt::Long;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg=s","debug|d");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my %cfg = Config::General::ParseConfig($opt{cfg});


###############################################################################
# Main Program
###############################################################################

my @folders = <$cfg{data_folder}/*/raw_data>;

if ($opt{debug}) {
	print "Example file: ",join("\n",@folders[0..0]),"\n";
	#@folders = $folders[0];
}


for (@folders) {
	$_ =~ /($cfg{data_folder}\/(\d+)\/)/;
	my $picture_num = $2;
	
	my $plots_folder = "$1/plots";
	mkpath($plots_folder);

	my $file_name = &write_data_file("$_/Centroid_dist_from_edge","$_/Area");

 	&build_and_execute_gnuplot_file("cent_area",$file_name,"$plots_folder/Cent_dist_vs_area.png");

	$file_name = &write_data_file("$_/Centroid_dist_from_edge", "$_/Average_adhesion_signal", "$_/Variance_adhesion_signal");
	
 	&build_and_execute_gnuplot_file("cent_sig",$file_name,"$plots_folder/Cent_dist_vs_sig.png");
}


###############################################################################
# Functions
###############################################################################

sub gather_data_from_matlab_file {
	my ($file) = @_;
	
	open INPUT, "$file" or die "Problem opening $file"; 
	my @in = split("   ",<INPUT>);
	close INPUT;
	return @in;
}

sub build_data_file_from_matlab {
	my ($file_1, $file_2) = @_;

	my @in_1 = &gather_data_from_matlab_file($file_1);
	my @in_2 = &gather_data_from_matlab_file($file_2);
	
	if ($opt{debug}) {
		if (scalar(@in_1) != scalar(@in_2)) {
		print "Problem with files: $file_1, $file_2, the number of entries do ",
			  "not match, they are: ", scalar(@in_1), " and ", scalar(@in_2);
		}
		if ($in_1[0] =~ /\s+/ || $in_2[0] =~ /\s+/) {
			print "First entry in one of the files $file_1 or $file_2 is not empty\n";
		}
	}

	my @output;
	for (1..$#in_1) {
		push @output, sprintf("%f",$in_1[$_]). "	".sprintf("%f",$in_2[$_]);
	}
	return @output;
}

sub build_data_file_from_matlab_with_var {
	my ($file_1, $file_2, $file_3) = @_;

	my @first_out = &build_data_file_from_matlab($file_1,$file_2);

	my @var_in = &gather_data_from_matlab_file($file_3);
	
	if ($opt{debug}) {
		if (scalar(@first_out) != scalar(@var_in) - 1) {
			print "Problem with data the number of entries do not match, they are: ",
				  scalar(@first_out), " and ", scalar(@var_in);
		}
		if ($var_in[0] =~ /\s+/) {
			print "First entry in one of the files $file_3 is not empty\n";
		}
	}

	my @output;
	for (0..$#first_out) {
		my $temp = "$first_out[$_]	". sprintf("%f",sqrt($var_in[$_ + 1]));
		push @output, $temp;
	}
	return @output;
}

sub write_data_file {
	my ($file_1, $file_2, $file_3);
	if (scalar(@_) == 2) {
		($file_1, $file_2) = @_;
	} elsif (scalar(@_) == 3) {
		($file_1, $file_2, $file_3) = @_;
	} else {
		die "Wrong number of parameters(",scalar(@_),") to 'write_data_file', expected 2 or 3.";
	}

	my ($temp_h,$file_name) = tempfile();
	my (@header_line,@data_lines);
	
	$file_1 =~ /.*\/(.*)/;
	push @header_line, $1;
	$file_2 =~ /.*\/(.*)/;
	push @header_line, $1;
	if (scalar(@_) == 3) {
		$file_3 =~ /.*\/(.*)/;
		push @header_line, $1;
	}
	
	if (scalar(@_) == 2) {
		@data_lines = &build_data_file_from_matlab($file_1,$file_2);
	}
	elsif (scalar(@_) == 3) {
		@data_lines = &build_data_file_from_matlab_with_var($file_1,$file_2,$file_3);
	}

	if ($opt{debug}) {
		open OUTPUT, ">$_/".join("_",@header_line)."_data_file" or die "Unable to open the output file";
		print OUTPUT "#", join("	",@header_line), "\n", join("\n",@data_lines);
		close OUTPUT;
	}

	print $temp_h "#", join("	",@header_line), "\n", join("\n",@data_lines);
	close $temp_h;
	return $file_name;
}

sub build_and_execute_gnuplot_file {
	my ($plot_type,$data_file_name,$out_file) = @_;
	
	my ($title,$xlabel,$ylabel);

	if ($plot_type =~ /cent_area/i) {
		$title  = "Centroid Distance from Cell Edge vs. Area of Identified Adhesions";
		$xlabel = "Centroid Distance from Cell Edge";
		$ylabel = "Area of Identified Adhesions";
	} elsif ($plot_type =~ /cent_sig/) {
		$title  = "Centroid Distance from Cell Edge vs. Average Adhesion Signal";
		$xlabel = "Centroid Distance from Cell Edge";
		$ylabel = "Average Normalized Fluorescence Signal";
	}
	else {
		die "$plot_type not defined for building a plot file";
	}
	
	my @out = ("set terminal png","set key off","set size square",
			   "set title \"$title\"","set xlabel \"$xlabel\"",
			   "set ylabel \"$ylabel\"","set output '$out_file'",
			   "plot '$data_file_name'");

	my ($temp_h,$gnuplot_file_name) = tempfile();
	print $temp_h join("\n",@out);
	close $temp_h;

	system "gnuplot $gnuplot_file_name";
}	
