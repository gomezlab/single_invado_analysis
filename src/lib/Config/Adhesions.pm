#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
use strict;
use warnings;
use File::Spec;

use base qw(Config::General);

my %derived_vars = (
    individual_results_folder => [qw(results_folder exp_name single_image_folder)],
    exp_results_folder        => [qw(results_folder exp_name)],
    exp_data_folder           => [qw(data_folder exp_name)],
);

###############################################################################
# Module Definition
###############################################################################
package Config::Adhesions;

sub new {
    my $class        = $_[0];
    my %opt          = %{ $_[1] };
    my @needed_files = ();
    @needed_files = @{ $_[2] } if (scalar(@_) > 2);

    my $conf = new Config::General(
        -ConfigFile            => $opt{cfg},
        -MergeDuplicateOptions => 1,
        -IncludeRelative       => 1,
    );
    my %cfg = $conf->getall;
    %{ $cfg{opt} }          = %opt;
    @{ $cfg{needed_files} } = @needed_files;

    my $cfg_ref = \%cfg;

    bless $cfg_ref, $class;

    $cfg_ref->build_derived_parameters;
    $cfg_ref->collect_cfg_info_from_files;

    return $cfg_ref;
}

sub collect_cfg_info_from_files {
    my $cfg = shift;

    #check to see if a file for the frames that should be excluded from the
    #analysis is included, if it is, collect the data from it, otherwise, set
    #exclude_image_nums to 0

    if (defined $cfg->{exclude_file}) {
        my $exclude_file = File::Spec->catfile($cfg->{data_folder}, $cfg->{exp_name}, $cfg->{exclude_file});
        open EX_INPUT, $exclude_file or die "Can't open the specified exclude file: $exclude_file";
        my $temp_line = <EX_INPUT>;
        close EX_INPUT;

        if (not($temp_line)) {
            @{ $cfg->{exclude_image_nums} } = (0);
        } else {
            chomp($temp_line);
            @{ $cfg->{exclude_image_nums} } = split(",", $temp_line);
        }
    } else {
        @{ $cfg->{exclude_image_nums} } = (0);
    }
}

sub build_derived_parameters {
    my $cfg = $_[0];
    foreach my $this_key (keys %derived_vars) {
        my $all_present = 1;
        foreach my $var_name (@{ $derived_vars{$this_key} }) {
            $all_present = 0 if not defined $cfg->{$var_name};
        }

        if ($all_present) {
            $cfg->{$this_key} = File::Spec->catdir(map $cfg->{$_}, @{ $derived_vars{$this_key} });
        }
    }
}

sub get_cfg_hash {
    my $self = shift;
    return map { $_ => ${$self}{$_} } keys %{$self};
}

1;
