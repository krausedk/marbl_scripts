#!/usr/bin/perl

use strict;
use warnings;

my %gatk = ();
my $file = 'combined_deepb1_gatk_error_rates_noins.out';
open my $info, $file or die "Could not open $file: $!";
while(defined(my $line = <$info>)) {
	chomp $line;
	my @fields = split(/\t/, $line);
	my $qual = ord($fields[0])-33;
	$gatk{$qual}{$fields[1]}{$fields[2]}{$fields[3]}{$fields[4]} = {'pos' => $fields[5], 'neg' => $fields[6]};
}
my @dinucs = ('0', '1', 'AA', 'AT', 'AG', 'AC', 'TA', 'TT', 'TG', 'TC', 'GA', 'GT', 'GG', 'GC', 'CA', 'CT', 'CG', 'CC', 'NA', 'NT', 'NG', 'NC', 'NN', 'AN', 'TN', 'CN', 'GN');
my @nucs = ('A', 'T', 'G', 'C');
my %rates = ();
foreach my $qual (30..41) {
	foreach my $pos (2..149) {
		foreach my $dinuc (@dinucs) {
			foreach my $refnuc (@nucs) {
				foreach my $readnuc (@nucs) {
					next if $refnuc eq $readnuc;
					my $neg_refnuc = $refnuc;
					$neg_refnuc =~ tr/ATCG/TAGC/;
					my $neg_readnuc = $readnuc;
					$neg_readnuc =~ tr/ATCG/TAGC/;
					my $sum_calls = 0;
					my $sum_error = 0;
					foreach my $sub_ref_nuc (@nucs) {
						if (exists $gatk{$qual}{$pos}{$dinuc}{$refnuc}{$sub_ref_nuc}) {
							$sum_calls += $gatk{$qual}{$pos}{$dinuc}{$refnuc}{$sub_ref_nuc}{'pos'};
						}
						if (exists $gatk{$qual}{$pos}{$dinuc}{$neg_refnuc}{$sub_ref_nuc}) {
							$sum_calls += $gatk{$qual}{$pos}{$dinuc}{$neg_refnuc}{$sub_ref_nuc}{'neg'};
						}
					}
					if (exists $gatk{$qual}{$pos}{$dinuc}{$refnuc}{$readnuc}) {
						$sum_error += $gatk{$qual}{$pos}{$dinuc}{$refnuc}{$readnuc}{'pos'};
					}
					if (exists $gatk{$qual}{$pos}{$dinuc}{$neg_refnuc}{$neg_readnuc}) {
						$sum_error += $gatk{$qual}{$pos}{$dinuc}{$neg_refnuc}{$neg_readnuc}{'neg'};
					}
					my $errortype = join('', $refnuc, $readnuc);
					if ($sum_calls > 0) {
						my $qual_trans = chr($qual+33);
						$rates{$qual_trans}{$pos}{$dinuc}{$errortype} = $sum_error/$sum_calls;
					}
				}
			}
		}
	}
}
open (OUT, ">deepb1_gatk_rates_for_ancova.txt");
print OUT "qual\tpos\tnucs\trate\n";
foreach my $qual (sort keys %rates) {
	foreach my $pos (sort {$a <=> $b} keys %{ $rates{$qual} }) {
		foreach my $dinuc (sort keys %{ $rates{$qual}{$pos} }) {
			foreach my $type (sort keys %{ $rates{$qual}{$pos}{$dinuc} }) {
				my $newname = join('', $dinuc, $type);
				my $number = sprintf("%.8f", $rates{$qual}{$pos}{$dinuc}{$type});
				print OUT "$qual\t$pos\t$newname\t$number\n";
			}
		}
	}
}	
