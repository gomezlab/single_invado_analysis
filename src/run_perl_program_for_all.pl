#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use File::Basename;
use File::Spec::Functions;
use Cwd 'abs_path';
use Getopt::Long;

use lib "lib";
use Config::Adhesions;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg|config=s", "debug|d", "program|p=s");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};
die "Can't find program to execute on the command line" if not exists $opt{program};

print "Collecting Configuration\n" if $opt{debug};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg = $ad_conf->get_cfg_hash;


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
    system("./$program_base -cfg $_ $debug_string");
}
