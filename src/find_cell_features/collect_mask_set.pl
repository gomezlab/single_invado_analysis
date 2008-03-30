#!/usr/bin/env perl

=head1 Name

collect_mask_set.pl - Finds the cell bodies in a set of mRFP images

=head1 Synopsis

collect_mask_set.pl -cfg adhesion_config_file -debug

=cut

###############################################################################
# Global Variables and Modules
###############################################################################

use strict;
use File::Path;
use Image::ExifTool;
use Math::Matlab::Local;
use Getopt::Long;
use File::Spec::Functions;

use lib "../lib";
use Config::Adhesions;
use Image::Stack;
use Math::Matlab::Extra;

#Perl built-in variable that controls buffering print output, 1 turns off 
#buffering
$| = 1;

my %opt;
$opt{debug} = 0;
GetOptions(\%opt, "cfg=s", "debug|d");

die "Can't find cfg file specified on the command line" if not exists $opt{cfg};

my $ad_conf = new Config::Adhesions(\%opt);
my %cfg = $ad_conf->get_cfg_hash;

my $matlab_wrapper;
if (defined $cfg{matlab_executable}) {
    $matlab_wrapper = Math::Matlab::Local->new({ cmd => "$cfg{matlab_executable} -nodisplay -nojvm -nosplash", });
} else {
    $matlab_wrapper = Math::Matlab::Local->new();
}

###############################################################################
# Main Program
###############################################################################

my @image_files = <$cfg{individual_results_folder}/*/$cfg{cell_mask_file}>;

if ($opt{debug}) {
    if (scalar(@image_files) > 1) {
        print "Cell mask files found: $image_files[0] - $image_files[$#image_files]\n";
    } else {
        print "Cell mask file found: $image_files[0]\n";
    }
}

my @matlab_code = &create_matlab_code;

my $error_file = catdir($cfg{exp_results_folder},$cfg{matlab_errors_folder},$cfg{cell_mask_errors_file});

&Math::Matlab::Extra::execute_commands($matlab_wrapper,\@matlab_code,$error_file);

###############################################################################
#Functions
###############################################################################

sub create_matlab_code {
    my @matlab_code;

    foreach my $file_name (@image_files) {
        my $i_num;
        if ($file_name =~ /$cfg{individual_results_folder}\/(\d+)\//) {
            $i_num = $1;
        } else {
            die "Skipping file: $file_name\n",
                "Unable to find image number.";
            next;
        }

        next if grep $i_num == $_, @{ $cfg{exclude_image_nums} };
        
        my $padded_num = sprintf("%0" . length(scalar(@image_files)) . "d", $i_num);

        my $out_file = catfile(dirname($file_name), "cell_mask.png");
        $matlab_code[0] .= "find_cell_mask('$file_name','$out_file')\n";
    }
    
    return @matlab_code;
}
