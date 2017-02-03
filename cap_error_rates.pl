#!/usr/bin/perl

use strict;
use warnings;

my %gatk = ();
my $file = 'deepb1_gatk_recalibrated_interpolated_error_rates.txt';
open my $info, $file or die "Could not open $file: $!";
while(defined(my $line = <$info>)) {
	chomp $line;
	my @fields = split(/\t/, $line);
	$gatk{$fields[0]}{$fields[1]}{$fields[2]} = $fields[3];
}
open(OUT, ">deepb1_gatk_recalibrated_interpolated_capped_error_rates.txt");
foreach my $nucs (sort keys %gatk) {
	foreach my $qual (sort keys %{ $gatk{$nucs} }) {
		foreach my $pos (sort {$a <=> $b} keys %{ $gatk{$nucs}{$qual} }) {
			if ($gatk{$nucs}{$qual}{$pos} < .0000011) {
				print OUT "$nucs\t$qual\t$pos\t.0000011\n";
			} elsif ($gatk{$nucs}{$qual}{$pos} > 0.1) {
				print OUT "$nucs\t$qual\t$pos\t0.1\n";
			} else {
				print OUT "$nucs\t$qual\t$pos\t$gatk{$nucs}{$qual}{$pos}\n";
			}
		}
	}
}	
