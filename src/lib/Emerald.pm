#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
{
package Emerald;
use strict;
use warnings;
use File::Spec;
use Data::Dumper;

###############################################################################
# Module Definition
###############################################################################
my %opt = ("queue" => "idle", "output_file" => "out.txt", 
           "error_file" => "error.txt", "folder" => "./",
           "runtime" => "24",);

sub send_LSF_commands {
    my @commands = @{$_[0]};
    
    foreach (@commands) {
        system("$_");
    }
}

sub create_LSF_Matlab_commands {
    my @commands = @{$_[0]};
    if (scalar(@_) > 1) {
        my %temp = %{$_[1]};
        $opt{$_} = $temp{$_} foreach (keys %temp);
    }
    $opt{output_file} = File::Spec->catfile($opt{folder},$opt{output_file});
    $opt{error_file}  = File::Spec->catfile($opt{folder},$opt{error_file});
    unlink($opt{output_file}, $opt{error_file});

    my $bsub_command   = "bsub -R RH5 -q $opt{queue} -o $opt{output_file} -e $opt{error_file} -We $opt{runtime}";
    my $matlab_command = "/afs/isis.unc.edu/pkg/matlab-2008a/matlab -nodisplay -nojvm -nosplash -r";
    $matlab_command = "matlab -nodisplay -nojvm -nosplash -logfile $opt{output_file} -r";

    @commands = map { split(/\n/, $_) } @commands;
    @commands = map { "$bsub_command $matlab_command \"$_\""} @commands;
    return @commands;   
}

sub create_general_LSF_commands {
    my @commands = @{$_[0]};
    if (scalar(@_) > 1) {
        my %temp = %{$_[1]};
        $opt{$_} = $temp{$_} foreach (keys %temp);
    }
    $opt{output_file} = File::Spec->catfile($opt{folder},$opt{output_file});
    $opt{error_file}  = File::Spec->catfile($opt{folder},$opt{error_file});
    unlink($opt{output_file}, $opt{error_file});

    my $bsub_command   = "bsub -q $opt{queue} -o $opt{output_file} -e $opt{error_file} -We $opt{runtime}";
    
    @commands = map { split(/\n/, $_) } @commands;
    @commands = map { "$bsub_command \"$_\""} @commands;
    return @commands;   
}

1;
}
