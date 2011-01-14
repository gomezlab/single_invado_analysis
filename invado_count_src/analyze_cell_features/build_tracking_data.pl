#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use lib "../lib";
use lib "../lib/perl";

use strict;
use File::Path;
use File::Basename;
use File::Spec::Functions;
use Getopt::Long;
use Data::Dumper;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);
use Text::CSV;
use IO::File;
use Math::Matrix;

use Config::Adhesions;
use Text::CSV::Simple::Extra;
use Image::Data::Collection;
use Emerald;
use FA_job;

#Perl built-in variable that controls buffering print output, 1 turns off
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
$opt{output} = "data.stor";
GetOptions(\%opt, "cfg|config=s", "debug|d", "output|o=s", "image_num=s", 
                  "lsf|l") or die;

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

print "Collecting Configuration\n" if $opt{debug};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg = $ad_conf->get_cfg_hash;

###############################################################################
#Main Program
###############################################################################

my %data_sets;
if ($opt{lsf}) {
    my @commands = ("$0 -cfg $opt{cfg} -o $opt{output}");
    
    $opt{error_folder} = catdir($cfg{exp_results_folder}, $cfg{errors_folder}, 'tracking_data');
    $opt{resource} = "mem32";
    $opt{queue} = "week";
    if (defined $cfg{job_group}) {
        $opt{job_group} = $cfg{job_group};
    }
    
    &FA_job::send_general_program(\@commands,\%opt);

    exit(0);
} else {
    print "\n\nGathering Data Files\n" if $opt{debug};

    my @data_files;
    push @data_files, @{$cfg{general_data_files}};
    push @data_files, @{$cfg{tracking_files}};

    %data_sets = &Image::Data::Collection::gather_data_sets(\%cfg, \%opt, \@data_files);

    print "\n\nPackaging Comparison Matrices\n" if $opt{debug};
    &make_comp_matices;
}

###############################################################################
#Functions
###############################################################################

#######################################
# Process Data Sets
#######################################
sub make_comp_matices {
    my @data_keys = sort { $a <=> $b } keys %data_sets;
    for (0 .. $#data_keys) {
        if (exists $opt{image_num} && $data_keys[$_] != $opt{image_num}) {
            next;
        }

        #The last image can not be compared to a future image, so we skip
        #calculations on it, but still save the image data if the output
        #option is specified
        if ($_ == $#data_keys) {
            if (defined $opt{output}) {
                my $key_1 = $data_keys[$_];
        		delete $data_sets{$key_1}{PixelIdxList};
                nstore \%{ $data_sets{$key_1} }, catfile($cfg{individual_results_folder}, $key_1, $opt{output});
                delete $data_sets{$key_1};
            }
            next;
        }

        #These are the keys we will use for all the subsequent matrix creation
        my ($key_1, $key_2) = @data_keys[ $_, $_ + 1 ];

        if (defined $opt{output}) {
            nstore \%{ $data_sets{$key_1} }, catfile($cfg{individual_results_folder}, $key_1, $opt{output});
            delete $data_sets{$key_1};
        }
    }
}
