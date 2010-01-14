#!/usr/bin/perl -w

use Getopt::Long;
use File::Find;

my %opt;
$opt{debug} = 0;
$opt{resample} = 300;
$opt{output_format} = "tiff";
$opt{dir} = '.';
GetOptions(\%opt, "debug|d", "resample=s", "dir=s", "output_format=s");

find(\&find_vect_files, $opt{dir});

sub find_vect_files {
    if (   $_ =~ /.*\.svg/
        && not($File::Find::name =~ /no_conv/)
        && not($File::Find::name =~ /\.svn/)) {
       &convert_vect_images($_);
    }
}

sub convert_vect_images {
    my $image_name = $_[0];
    my $output_name = $image_name;
    $output_name =~ s/(.*)\..*/$1.$opt{output_format}/;
    
    my $command = "time inkscape $image_name -d $opt{resample} --export-png=$output_name --export-background-opacity=1.0";
    
	if ($opt{debug}) {
        print $command, "\n";
    } else {
        system($command);
    }
}
