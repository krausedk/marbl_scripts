#!/usr/bin/perl

use strict;
use warnings;

my %combined = ();
my @files = glob("deepb1_pe_to_M1627_novoalign.sam.*");
foreach my $file (@files) {
	open my $info, $file or die "Could not open $file: $!";
	while(defined(my $line = <$info>)) {
		chomp $line;
		my @fields = split(/\t/, $line);
		$combined{$fields[0]}{$fields[1]}{$fields[2]}{$fields[3]}{$fields[4]}{'pos'} += $fields[5];
		$combined{$fields[0]}{$fields[1]}{$fields[2]}{$fields[3]}{$fields[4]}{'neg'} += $fields[6];
	}
}
open(OUT, ">combined_deepb1_gatk_error_rates.out");
foreach my $qual (sort keys %combined) {
	foreach my $pos (sort {$a <=> $b} keys %{ $combined{$qual} }) {
		foreach my $dinuc (sort keys %{ $combined{$qual}{$pos} }) {
			foreach my $nuc1 (sort keys %{ $combined{$qual}{$pos}{$dinuc} }) {
				foreach my $nuc2 (sort keys %{ $combined{$qual}{$pos}{$dinuc}{$nuc1} }) {
					print OUT "$qual\t$pos\t$dinuc\t$nuc1\t$nuc2\t$combined{$qual}{$pos}{$dinuc}{$nuc1}{$nuc2}{'pos'}\t$combined{$qual}{$pos}{$dinuc}{$nuc1}{$nuc2}{'neg'}\n";
				}
			}
		}
	}
}


