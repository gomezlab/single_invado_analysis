#!/usr/bin/env perl

###############################################################################
# Global Variables and Modules
###############################################################################

use strict;
use warnings;
use File::Path;
use File::Basename;
use Math::Matlab;

package Math::Matlab::Extra;

###############################################################################
# Functions
###############################################################################

sub execute_commands {
    my $matlab_object = $_[0];
    my @matlab_code = @{$_[1]};
    my $error_file = $_[2];

    unlink($error_file) if (-e $error_file);

    foreach my $command (@matlab_code) {
        if (not($matlab_object->execute($command))) {
            &File::Path::mkpath(&File::Basename::dirname($error_file));
            open ERR_OUT, ">>$error_file" or die "Error in opening Matlab Error file: $error_file.";
            print ERR_OUT $matlab_object->err_msg;
            print ERR_OUT "\n\nMATLAB COMMANDS\n\n$command";
            close ERR_OUT;

            $matlab_object->remove_files;
        }
    }
}

1;
