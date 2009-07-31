#!/usr/bin/perl -w

use Getopt::Long;
use File::Find;

my %opt;
$opt{debug} = 0;
$opt{resample} = 300;
$opt{dir} = '.';
GetOptions(\%opt, "debug|d", "resample=s", "dir=s");

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
    my $png_name = $image_name;
    $png_name =~ s/(.*)\..*/$1.png/;
    
    my $command = "time inkscape $image_name -d $opt{resample} --export-png=$png_name --export-background-opacity=1.0";
    
	if ($opt{debug}) {
        print $command, "\n";
    } else {
        system($command);
    }
}
