#!/usr/bin/perl

use strict;
use warnings;

my %regression = ();
my $file = "deepb1_gatk_polynomial_regressions.out";
open my $info, $file or die "Could not open $file: $!";
while(defined(my $line = <$info>)) {
	chomp $line;
	my @fields = split(/\t/, $line);
	for my $i (2..$#fields) {
		$regression{$fields[0]}{$fields[1]}{$i-2} = $fields[$i];
	}
}
my %recal = ();
foreach my $nucs (sort keys %regression) {
	foreach my $qual (sort {$a <=> $b} keys %{ $regression{$nucs} }) {
		foreach my $pos (2..149) {
			$recal{$nucs}{$qual}{$pos} = 0;
			foreach my $order (sort {$a <=> $b} keys %{ $regression{$nucs}{$qual} }) {
				my $coeff = $regression{$nucs}{$qual}{$order};
				$recal{$nucs}{$qual}{$pos} += $coeff*($pos**$order);
			}
		}
	}
}
open (OUT, ">deepb1_gatk_recalibrated_interpolated_error_rates.txt");
foreach my $nucs (sort keys %recal) {
	foreach my $qual (sort {$a <=> $b} keys %{ $recal{$nucs} }) {
		foreach my $pos (2..149) {
			print OUT "$nucs\t$qual\t$pos\t$recal{$nucs}{$qual}{$pos}\n";
		}
	}
}
