#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

# recalibrate the error rates in a sample using the snp_database_file #

my %internal_pos = ();
my $snpfile = 'same_snp_database_final.out';
open my $snpinfo, $snpfile or die "Could not open $snpfile: $!";
while(defined(my $line = <$snpinfo>)) {
	chomp $line;
	my @fields = split(/\t/, $line);
	$internal_pos{$fields[0]} = $fields[1];
}
my $file = "$ARGV[0]";
open my $info, $file or die "Could not open $file: $!";
my %snps = ();
my %calls = ();
my %deletions = ();
my %insertions = ();
my $loopcounter = 1;
while(defined(my $line = <$info>)) {
	next if $. < $ARGV[2];
	last if $. > $ARGV[3];
	print "$loopcounter";
	chomp($line);
	my @samline = split(/\s/, $line);
	my $leftpos = $samline[3];
	my $flag = $samline[1];
	if ($flag !~ m/99|147|83|163/) {
		print "impropermap\n";
		$loopcounter ++;
		next;
	}
	my $mapq = $samline[4];
	if($mapq <41) {
		print "mapqlow\n";
		$loopcounter ++;
		next;
	}
	my $cigar = $samline[5];
	if ($cigar =~ m/^(\d+)I/) {
		print "damn you novoalign!\n";
		$loopcounter ++;
		next;
	}
	my $nucs = $samline[9];
	my $qual = $samline[10];
	my @nucarray = split(//, $nucs);
	my @dinucarray;
	my @pos_in_read;
	if ($flag & 16) {
		my $revnucstring = reverse($nucs);
		$revnucstring =~ tr/ATCG/TAGC/;
		push(@dinucarray, '0');
		push(@dinucarray, '1');
		for (my $i = 0; $i < $#nucarray-1; $i++) {
			push(@dinucarray, substr($revnucstring, $i, 2));
		}
		@dinucarray = reverse(@dinucarray);
		for my $i (0..$#dinucarray) {
			push(@pos_in_read, $i);
		}
		@pos_in_read = reverse(@pos_in_read);
	} else {
		push(@dinucarray, '0');
		push(@dinucarray, '1');
		for (my $i = 0; $i < $#nucarray-1; $i++) {
			push(@dinucarray, substr($nucs, $i, 2));
		}
		for my $i (0..$#dinucarray) {
			push(@pos_in_read, $i);
		}
	}
	my @qualarray = split(//, $qual);
	my @posarray = ();
	my $nuccounter = 0;
	my $counter = $leftpos;
	my @derived_nucarray = ();
	my @derived_qualarray = ();
	my @derived_dinucarray = ();
	my @derived_posarray = ();
	my $mdtag;
	my $inscount = 0;
	my $delcount = 0;
	my $snpcount = 0;
	my $Ncount = 0;
	for my $i(11..$#samline) {
		if($samline[$i] =~ m/MD/) {
			my @md = split(/:/, $samline[$i]);
			$mdtag = $md[2];
		}
	}
	$cigar =~ s/([MID])/$1,/g;
	my @cigedit = split(/,/, $cigar);
	for my $i(0..$#cigedit) {
		if($cigedit[$i] =~ m/(\d+)M/) {
			for my $j (1..$1) {
				push(@posarray, $counter);
				push(@derived_nucarray, $nucarray[$nuccounter]);
				push(@derived_dinucarray, $dinucarray[$nuccounter]);
				push(@derived_qualarray, $qualarray[$nuccounter]);
				push(@derived_posarray, $pos_in_read[$nuccounter]);
				$nuccounter ++;
				$counter ++;
			}
		} elsif ($cigedit[$i] =~ m/(\d+)I/) {
			my @insertion;
			$inscount ++;
			for my $k (1..$1) {
				push(@insertion, $nucarray[$nuccounter]);
				$nuccounter ++;
			}
			my $joinedarray = join('', @insertion);
			my $ins = join('', $derived_nucarray[-1], '+', $1, $joinedarray);
			splice(@derived_nucarray, -1);
			push(@derived_nucarray, $ins);
		} elsif ($cigedit[$i] =~ m/(\d+)D/) {
			$delcount ++;
			for my $l (1..$1) {
				push(@posarray, $counter);
				push(@derived_nucarray, 0);
				push(@derived_qualarray, '#');
				push(@derived_dinucarray, 0);
				push(@derived_posarray, $pos_in_read[$nuccounter]);
				$counter ++;
			}
		}
	}
	#snps stored by basequality score #
	for my $i(0..$#posarray) {
		if (exists $internal_pos{$posarray[$i]}) {
			if($derived_nucarray[$i] eq $internal_pos{$posarray[$i]}) {
				if ($flag & 16) {
					$calls{$derived_qualarray[$i]}{$derived_posarray[$i]}{$derived_dinucarray[$i]}{$derived_nucarray[$i]}{$derived_nucarray[$i]}{'neg'} ++;
				} else {
					$calls{$derived_qualarray[$i]}{$derived_posarray[$i]}{$derived_dinucarray[$i]}{$derived_nucarray[$i]}{$derived_nucarray[$i]}{'pos'} ++;
				}
			} else {
				if ($flag & 16) {
					$calls{$derived_qualarray[$i]}{$derived_posarray[$i]}{$derived_dinucarray[$i]}{$internal_pos{$posarray[$i]}}{$derived_nucarray[$i]}{'neg'} ++;
				} else {
					$calls{$derived_qualarray[$i]}{$derived_posarray[$i]}{$derived_dinucarray[$i]}{$internal_pos{$posarray[$i]}}{$derived_nucarray[$i]}{'pos'} ++;
				}
			}
		}
	}
	print "success\n";
	$loopcounter ++;
}
open(OUT, ">$ARGV[1]");
foreach my $qual (sort keys %calls) {
	foreach my $readpos (sort {$a <=> $b} keys %{ $calls{$qual} }) {
		foreach my $dinuc (sort keys %{ $calls{$qual}{$readpos} }) {
			foreach my $refbase (sort keys %{ $calls{$qual}{$readpos}{$dinuc} }) {
				foreach my $readbase (sort keys %{ $calls{$qual}{$readpos}{$dinuc}{$refbase} }) {
					print OUT "$qual\t$readpos\t$dinuc\t$refbase\t$readbase\t$calls{$qual}{$readpos}{$dinuc}{$refbase}{$readbase}{'pos'}\t$calls{$qual}{$readpos}{$dinuc}{$refbase}{$readbase}{'neg'}\n";
				}
			}
		}
	}
}
