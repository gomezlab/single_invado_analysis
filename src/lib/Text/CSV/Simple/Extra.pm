#!/usr/bin/perl -w

###############################################################################
# Global Variables and Modules
###############################################################################
{ 
package Text::CSV::Simple::Extra;
use strict;
use warnings;
use Text::CSV;


our @EXPORT = qw( output_mat_csv );
use Exporter;
our @ISA = qw(Exporter);
###############################################################################
# Functions
###############################################################################

sub output_mat_csv {
    my @mat = @{$_[0]};
    my $file = $_[1];

    my $out_hand = new IO::File ">" . $file or die "Unable to create csv file: $file";
    my $csv = Text::CSV->new();

    for my $i (0 .. $#mat) {
        $csv->print($out_hand, \@{ $mat[$i] });
        print $out_hand "\n";
    }

    $out_hand->close;
}

}
1;
