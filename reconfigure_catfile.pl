#!/usr/bin/perl

use strict;
use warnings;

open(OUT, ">$ARGV[1]");
my %revised=();
my $file = "$ARGV[0]";
open my $info, $file or die "Could not open $file: $!";
while(defined(my $line = <$info>)) {
	chomp $line;
	my @fields = split(/\t/, $line);
	my @refnucs = split(/\//, $fields[2]);
	my @donornucs = split(/\//, $fields[3]);
	my @readnucs = split(/\//, $fields[4]);
	my @dinucs = split(/\//, $fields[6]);
	my @newnucs;
	my @altnewnucs;
	for my $i (0..$#dinucs) {		
		if ($dinucs[$i] =~ m/,/) {
			my @split = split(/,/, $dinucs[$i]);
			my $nucs1 = join('', $split[0], $refnucs[$i], $donornucs[$i]);
			my $altnucs1 = join('', $split[0], $donornucs[$i], $refnucs[$i]);
			my $nucs2 = join('', $split[1], $refnucs[$i], $donornucs[$i]);
			my $altnucs2 = join('', $split[1], $donornucs[$i], $refnucs[$i]); 
			my $joined = join(',', $nucs1, $nucs2);
			my $altjoined = join(',', $altnucs1, $altnucs2);
			push(@newnucs, $joined);
			push(@altnewnucs, $altjoined);
		} else {
			my $joined = join('', $dinucs[$i], $refnucs[$i], $donornucs[$i]);
			my $altjoined = join('', $dinucs[$i], $donornucs[$i], $refnucs[$i]);
			push(@newnucs, $joined);
			push(@altnewnucs, $altjoined);
		}
	}
	my $newnucstring = join('/', @newnucs);
	my $altnewnucstring = join('/', @altnewnucs);
	print OUT "$fields[0]\t$fields[1]\t$newnucstring\t$altnewnucstring\t$fields[5]\t$fields[7]\t$fields[8]\t$fields[9]\t$fields[10]\n";
}
