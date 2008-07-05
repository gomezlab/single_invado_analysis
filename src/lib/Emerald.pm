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

our @EXPORT = qw(send_emerald_commands);
use Exporter;
our @ISA = qw(Exporter);


###############################################################################
# Module Definition
###############################################################################

sub send_emerald_commands {
    my @commands = @{$_[0]};
    
    foreach (@commands) {
        system("$_");
    }
}

sub create_emerald_Matlab_commands {
    my @commands = @{$_[0]};
    my %opt = ("queue" => "week", "output_file" => "out.txt", 
               "error_file" => "error.txt", "folder" => "./", "stdout" => 0);
    if (scalar(@_) > 1) {
        my %temp = %{$_[1]};
        $opt{$_} = $temp{$_} foreach (keys %temp);
    }
    $opt{output_file} = File::Spec->catfile($opt{folder},$opt{output_file});
    $opt{error_file}  = File::Spec->catfile($opt{folder},$opt{error_file});

    my $bsub_command   = "bsub -R blade -q $opt{queue} -o $opt{output_file} -e $opt{error_file}";
    my $matlab_command = "/afs/isis/pkg/matlab/matlab -nodisplay -nojvm -nosplash -r";

    @commands = map { split(/\n/, $_) } @commands;
    @commands = map { "$bsub_command $matlab_command \"$_\""} @commands;
    return @commands;   
}

sub create_general_emerald_command {
    my @commands = @{$_[0]};
    my %opt = ("queue" => "week", "output_file" => "out.txt", 
               "error_file" => "error.txt", "folder" => "./", "stdout" => 0);
    if (scalar(@_) > 1) {
        my %temp = %{$_[1]};
        $opt{$_} = $temp{$_} foreach (keys %temp);
    }
    $opt{output_file} = File::Spec->catfile($opt{folder},$opt{output_file});
    $opt{error_file}  = File::Spec->catfile($opt{folder},$opt{error_file});

    my $bsub_command   = "bsub -R blade -q $opt{queue} -o $opt{output_file} -e $opt{error_file}";
    
    @commands = map { split(/\n/, $_) } @commands;
    @commands = map { "$bsub_command \"$_\""} @commands;
    return @commands;   
}

1;
}
