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

sub send_emerald_commands {
    my @commands = @{$_[0]};
    my %opt = ("queue" => "week", "output_file" => "out.txt", "error_file" => "error.txt", "folder" => "./", "split_on_newlines" => 1);
    if (scalar(@_) > 1) {
        my %temp = %{$_[1]};
        $opt{$_} = $temp{$_} foreach (keys %temp);
    }
    $opt{output_file} = File::Spec->catfile($opt{folder},$opt{output_file});
    $opt{error_file} = File::Spec->catfile($opt{folder},$opt{error_file});

    my $initial_command_num = length(@commands);
    
    @commands = map { split(/\n/, $_) } @commands;
    @commands = join("; ", @commands) if not($opt{split_on_newlines});
    @commands = map { "bsub -R blade -q $opt{queue} -o $opt{output_file} -e $opt{error_file} /afs/isis/pkg/matlab/matlab -nodisplay -nojvm -nosplash -r \"$_\""} @commands;
    
    foreach (@commands) {
        system("$_");
    }
}

1;
}
