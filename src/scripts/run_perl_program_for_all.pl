#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Basename;
use File::Spec::Functions;
use Cwd 'abs_path';
use Getopt::Long;

use lib "../lib";
use Config::Adhesions qw(ParseConfig);

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
$opt{extra} = "";
GetOptions(\%opt, "cfg|config=s", "debug|d", "program|p=s", "extra|e=s", "run_all_debug");

die "Can't find cfg file ($opt{cfg}) specified on the command line" if not(exists $opt{cfg});
die "Can't find program to execute on the command line" if not(exists $opt{program});

print "Collecting Configuration\n" if $opt{debug};
my %cfg = ParseConfig(\%opt);

###############################################################################
#Main Program
###############################################################################
chdir(dirname($opt{program}));

my $program_base = basename($opt{program});

my $debug_string = ($opt{debug}) ? "-d" : "";

my $cfg_suffix = basename($opt{cfg});
$cfg_suffix =~ s/.*\.(.*)/$1/;

my @exp = <$cfg{data_folder}/*/*$cfg_suffix>;

foreach (@exp) {
    next if /config\/default/; 
    if ($opt{run_all_debug}) {
        print("./$program_base -cfg $_ $debug_string $opt{extra}\n");
    } else {
        system("./$program_base -cfg $_ $debug_string $opt{extra}");
    }
}
