#!/usr/bin/perl -w

use File::Find

$| = 1;

our @lines;

find(\&perl_file_add, ('../'));

my %use_hits;
foreach (@lines) {
	if ($_ =~ /^use/) {
		$use_hits{$_}++;
	}
}

die join("",sort keys %use_hits);

sub perl_file_add {
	if ($File::Find::name =~ /\.pl$/ &&
	    not($File::Find::name =~ /.*lib\/perl.*/)) {
		open INPUT, "$_" or die "$!";
		push @lines, <INPUT>;
		close INPUT;
	}
}
