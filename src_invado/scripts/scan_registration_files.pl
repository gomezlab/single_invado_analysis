#!/usr/bin/perl -w

use Text::CSV::Simple;

my $base_dir = "../../results/Invadopodia/";

my @files = <$base_dir/*/*/ind*/*/affine_matrix.csv>;

my %x_trans;
my %y_trans;

for (@files) {
	my $parser = Text::CSV::Simple->new;
	my @data = $parser->read_file($_);

	die "$_" if (scalar(@data) != 3);
	die "$_" if (scalar(@{$data[2]}) != 2);
	
	$x_trans{$data[2][0]}++;
	$y_trans{$data[2][1]}++;

	if (abs($data[2][0]) == 7 || abs($data[2][1]) == 7) {
		# print "$_\n\n";
	}
}
my $under_x = 0;
my $all = 0;
for (sort {$a <=> $b} keys %x_trans) {
	print "$_ $x_trans{$_}\n";
	$all += $x_trans{$_};

	if (abs($_) < 5) {
		$under_x += $x_trans{$_};
	}
}

print "\n\n";

my $under_y = 0;
my $all_2 = 0;
for (sort  {$a <=> $b} keys %y_trans) {
	print "$_ $y_trans{$_}\n";
	$all_2 += $y_trans{$_};

	if (abs($_) < 5) {
		$under_y += $y_trans{$_};
	}
}

print "$under_x $all: ", $under_x/$all, "\n";
print "$under_y $all_2: ", $under_y/$all_2;
