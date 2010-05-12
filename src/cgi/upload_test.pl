#!/usr/bin/perl -wT

use strict;
use File::Path;
use File::Basename;
use File::Spec::Functions;
use POSIX;
use CGI qw/:all/;
use CGI::Carp;

my $q = CGI->new();
print $q->header,                    # create the HTTP header
	  $q->start_html('hello world'), # start the HTML
	  $q->h1('hello world'),         # level 1 header
	  $q->end_html;                  # end the HTML
