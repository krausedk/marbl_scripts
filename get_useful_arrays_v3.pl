#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $file = $ARGV[0];
open my $info, $file or die "Could not open $file: $!";
my %snps = ();
my %calls = ();
my %deletions = ();
my %insertions = ();
my $loopcounter = 1;
my %snplist = ();
my $snpfile = 'snp_database_final_evencov_nodels.out';
open my $stuff, $snpfile or die "Could not open $snpfile: $!";
while(defined(my $line = <$stuff>)) {
	chomp($line);
	my @fields = split(/\t/, $line);
#	$donorlist{$fields[0]} = {'nuc' => $fields[1], 'reads' => $fields[2]};
	$snplist{$fields[0]} = {'donor' => $fields[2], 'ref' => $fields[1]};
}
open(OUT, ">$ARGV[1]");
while(defined(my $line1 = <$info>)) {
	next if $. < $ARGV[2];
	last if $. > $ARGV[3];
	defined(my $line2 = <$info>) or last;
	print "$loopcounter";
	chomp($line1);
	chomp($line2);
	my @samline1 = split(/\s/, $line1);
	my @samline2 = split(/\s/, $line2);
	my $leftpos1 = $samline1[3];
	my $leftpos2 = $samline2[3];
	my $title1 = $samline1[0];
	my $flag1 = $samline1[1];
	my $flag2 = $samline2[1];
	if(($flag1 =~ m/97|145|81|161|65|129|113|177/) or ($flag2 =~ m/97|145|81|161|65|129|113|177/)) {
		print "impropermap\n";
		$loopcounter ++;
		next;
	}
	my $mapq1 = $samline1[4];
	my $mapq2 = $samline2[4];
	if(($mapq1 < 41) or ($mapq2 < 41)) {
		print "mapqlow\n";
		$loopcounter ++;
		next;
	}
	my $nucs1 = $samline1[9];
	my $nucs2 = $samline2[9];
	my $qual1 = $samline1[10];
	my $qual2 = $samline2[10];
	my @nucarray1 = split(//, $nucs1);
	my @nucarray2 = split(//, $nucs2);
	my @qualarray1 = split(//, $qual1);
	my @qualarray2 = split(//, $qual2);
	### make the dinucarrays ###
	my @dinucarray1;
	my @dinucarray2;
	my @pos_in_read1;
	my @pos_in_read2;
	if ($flag1 & 16) {
		my $rev_nuc1string = reverse($nucs1);
		$rev_nuc1string =~ tr/ATCG/TAGC/;
		push(@dinucarray1, '0');
		push(@dinucarray1, '1');
		for (my $i = 0; $i < $#nucarray1-1; $i++) {
			push(@dinucarray1, substr($rev_nuc1string, $i, 2));
		}
		@dinucarray1 = reverse(@dinucarray1);
		for my $i (0..$#dinucarray1) {
			push(@pos_in_read1, $i);
		}
		@pos_in_read1 = reverse(@pos_in_read1);
		push(@dinucarray2, '0');
		push(@dinucarray2, '1');
		for (my $i = 0; $i < $#nucarray2-1; $i++) {
			push(@dinucarray2, substr($nucs2, $i, 2));
		}
		for my $i (0..$#dinucarray2) {
			push(@pos_in_read2, $i);
		}
	} else {
		my $rev_nuc2string = reverse($nucs2);
		$rev_nuc2string =~ tr/ATCG/TAGC/;
		push(@dinucarray2, '0');
		push(@dinucarray2, '1');
		for (my $i = 0; $i < $#nucarray2-1; $i++) {
			push(@dinucarray2, substr($rev_nuc2string, $i, 2));
		}
		@dinucarray2 = reverse(@dinucarray2);
		for my $i(0..$#dinucarray2) {
			push(@pos_in_read2, $i);
		}
		@pos_in_read2 = reverse(@pos_in_read2);
		push(@dinucarray1, '0');
		push(@dinucarray1, '1');
		for (my $i = 0; $i < $#nucarray1-1; $i++) {
			push(@dinucarray1, substr($nucs1, $i, 2));
		}
		for my $i(0..$#dinucarray1) {
			push(@pos_in_read1, $i);
		}
	}
	my @posarray1 = ();
	my @posarray2 = ();
	my $nuccounter1 = 0;
	my $nuccounter2 = 0;
	my $counter1 = $leftpos1;
	my $counter2 = $leftpos2;
	my @derived_nucarray1 = ();
	my @derived_nucarray2 = ();
	my @derived_qualarray1 = ();
	my @derived_qualarray2 = ();
	my @derived_dinucarray1 = ();
	my @derived_dinucarray2 = ();
	my @derived_posarray1 = ();
	my @derived_posarray2 = ();
	my @derived_nuc_context = ();
	my $cigar1 = $samline1[5];
	my $cigar2 = $samline2[5];
	#starting a cigar with an insertion is a new novoalign problem, so skip these lines for now#
	if ($cigar1 =~ m/^(\d+)I/) {
		print "damn you novoalign\n";
		$loopcounter ++;
		next;
	}
	if ($cigar2 =~ m/^(\d+)I/) {
		print "damn you novoalign\n";
		$loopcounter ++;
		next;
	}
	my $mdtag1;
	my $mdtag2;
	my $inscount1 = 0;
	my $delcount1 = 0;
	my $snpcount1 = 0;
	my $Ncount1 = 0;
	my $inscount2 = 0;
	my $delcount2 = 0;
	my $snpcount2 = 0;
	my $Ncount2 = 0;
	for my $i(11..$#samline1) {
		if($samline1[$i] =~ m/MD/) {
			my @md1 = split(/:/, $samline1[$i]);
			$mdtag1 = $md1[2];
		}
	}
	for my $i(11..$#samline2) {
		if($samline2[$i] =~ m/MD/) {
			my @md2 = split(/:/, $samline2[$i]);
			$mdtag2 = $md2[2];
		}
	}
	$cigar1 =~ s/([MID])/$1,/g;
	my @cigedit1 = split(/,/, $cigar1);
	for my $i(0..$#cigedit1) {
		if($cigedit1[$i] =~ m/(\d+)M/) {
			for my $j (1..$1) {
				push(@posarray1, $counter1);
				push(@derived_nucarray1, $nucarray1[$nuccounter1]);
				push(@derived_qualarray1, $qualarray1[$nuccounter1]);
				push(@derived_dinucarray1, $dinucarray1[$nuccounter1]);
				push(@derived_posarray1, $pos_in_read1[$nuccounter1]);
				$nuccounter1 ++;
				$counter1 ++;
			}
		} elsif ($cigedit1[$i] =~ m/(\d+)I/) {
			my @insertion1;
			$inscount1 ++;
			for my $k (1..$1) {
				push(@insertion1, $nucarray1[$nuccounter1]);
				$nuccounter1 ++;
			}
			my $joinedarray1 = join('', @insertion1);
			my $ins1 = join('', $derived_nucarray1[-1], '+', $1, $joinedarray1);
			splice(@derived_nucarray1, -1);
			push(@derived_nucarray1, $ins1);
		} elsif ($cigedit1[$i] =~ m/(\d+)D/) {
			$delcount1 ++;
			for my $l (1..$1) {
				push(@posarray1, $counter1);
				push(@derived_nucarray1, 0);
				push(@derived_qualarray1, 'H');
				push(@derived_dinucarray1, $dinucarray1[$nuccounter1]);
				push(@derived_posarray1, $pos_in_read1[$nuccounter1]);
				$counter1 ++;
			}
		}
	}
	$cigar2 =~ s/([MID])/$1,/g;
	my @cigedit2 = split(/,/, $cigar2);
	for my $i(0..$#cigedit2) {
		if($cigedit2[$i] =~ m/(\d+)M/) {
			for my $j (1..$1) {
				push(@posarray2, $counter2);
				push(@derived_nucarray2, $nucarray2[$nuccounter2]);
				push(@derived_qualarray2, $qualarray2[$nuccounter2]);
				push(@derived_dinucarray2, $dinucarray2[$nuccounter2]);
				push(@derived_posarray2, $pos_in_read2[$nuccounter2]);
				$nuccounter2 ++;
				$counter2 ++;
			}
		} elsif ($cigedit2[$i] =~ m/(\d+)I/) {
			my @insertion2;
			$inscount2 ++;
			for my $k (1..$1) {
				push(@insertion2, $nucarray2[$nuccounter2]);
				$nuccounter2 ++;
			}
			my $joinedarray2 = join('', @insertion2);
			my $ins2 = join('', $derived_nucarray2[-1], '+', $1, $joinedarray2);
			splice(@derived_nucarray2, -1);
			push(@derived_nucarray2, $ins2);
		} elsif ($cigedit2[$i] =~ m/(\d+)D/) {
			$delcount2 ++;
			for my $l (1..$1) {
				push(@posarray2, $counter2);
				push(@derived_nucarray2, 0);
				push(@derived_qualarray2, 'H');
				push(@derived_dinucarray2, $dinucarray2[$nuccounter2]);
				push(@derived_posarray2, $pos_in_read2[$nuccounter2]);
				$counter2 ++;
			}
		}
	}
	my @derived_mdarray1;
	my @derived_mdarray2;
	if (defined $mdtag1) {
		#snps by position
		$mdtag1 =~ s/\^//g;
		$mdtag1 =~ s/^0//g;
		$mdtag1 =~ s/([A-Z])0([A-Z])/$1$2/g;
		$mdtag1 =~ s/(\d+)/$1,/g;
		$mdtag1 =~ s/([A-Z])/$1,/g;	
		my @mdarray1 = split(/,/, $mdtag1);
		for my $i(0..$#mdarray1) {
			if ($mdarray1[$i] =~ m/(\d+)/) {
				for my $j(1..$1) {
					push(@derived_mdarray1, 1);
				}
			} elsif ($mdarray1[$i] =~ m/([A-Z])/) {
				push(@derived_mdarray1, $1);
			}
		}
	}
	if (defined $mdtag2) {
		#snps by position
		$mdtag2 =~ s/\^//g;
		$mdtag2 =~ s/^0//g;
		$mdtag2 =~ s/([A-Z])0([A-Z])/$1$2/g;
		$mdtag2 =~ s/(\d+)/$1,/g;
		$mdtag2 =~ s/([A-Z])/$1,/g;	
		my @mdarray2 = split(/,/, $mdtag2);
		for my $i(0..$#mdarray2) {
			if ($mdarray2[$i] =~ m/(\d+)/) {
				for my $j(1..$1) {
					push(@derived_mdarray2, 1);
				}
			} elsif ($mdarray2[$i] =~ m/([A-Z])/) {
				push(@derived_mdarray2, $1);
			}
		}
	}
	my %CIGMD1 = ();
	my %CIGMD2 = ();
	my %fusion = ();
	for my $i(0..$#posarray1) {
		if (((ord($derived_qualarray1[$i]) - 33) > 29) and ($derived_posarray1[$i] > 1) and ($derived_posarray1[$i] < 248)) {
			$CIGMD1{$posarray1[$i]} = {'nuc' => $derived_nucarray1[$i], 'qual' => $derived_qualarray1[$i], 'md' => $derived_mdarray1[$i], 'flag' => $flag1, 'dinuc' => $derived_dinucarray1[$i], 'readpos' => $derived_posarray1[$i]};
		}
	}
	for my $i(0..$#posarray2) {
		if (((ord($derived_qualarray2[$i]) - 33) > 29) and ($derived_posarray2[$i] > 1) and ($derived_posarray2[$i] < 248)) {
			$CIGMD2{$posarray2[$i]} = {'nuc' => $derived_nucarray2[$i], 'qual' => $derived_qualarray2[$i], 'md' => $derived_mdarray2[$i], 'flag' => $flag2, 'dinuc' => $derived_dinucarray2[$i], 'readpos' => $derived_posarray2[$i]};
		}
	}
	foreach my $key (sort {$a <=> $b} keys %CIGMD1) {
		if (exists $CIGMD2{$key}) {
			if ($CIGMD1{$key}{'nuc'} eq $CIGMD2{$key}{'nuc'}) {
				$fusion{$key} = {'nuc' => $CIGMD1{$key}{'nuc'}, 'qual' => join(",", $CIGMD1{$key}{'qual'}, $CIGMD2{$key}{'qual'}), 'md' => $CIGMD1{$key}{'md'}, 'flag' => join(",", $CIGMD1{$key}{'flag'}, $CIGMD2{$key}{'flag'}), 'count' => '2', 'dinuc' => join(",", $CIGMD1{$key}{'dinuc'}, $CIGMD2{$key}{'dinuc'}), 'readpos' => join(",", $CIGMD1{$key}{'readpos'}, $CIGMD2{$key}{'readpos'})};
				delete $CIGMD2{$key};
			} else {
				delete $CIGMD2{$key};
			}
		} else {
			$fusion{$key} = {'nuc' => $CIGMD1{$key}{'nuc'}, 'qual' => $CIGMD1{$key}{'qual'}, 'md' => $CIGMD1{$key}{'md'}, 'flag' => $CIGMD1{$key}{'flag'}, 'count' => '1', 'dinuc' => $CIGMD1{$key}{'dinuc'}, 'readpos' => $CIGMD1{$key}{'readpos'}};
		}
	}
	foreach my $key (sort {$a <=> $b} keys %CIGMD2) {
		$fusion{$key} = {'nuc' => $CIGMD2{$key}{'nuc'}, 'qual' => $CIGMD2{$key}{'qual'}, 'md' => $CIGMD2{$key}{'md'}, 'flag' => $CIGMD2{$key}{'flag'}, 'count' => '1', 'dinuc' => $CIGMD2{$key}{'dinuc'}, 'readpos' => $CIGMD2{$key}{'readpos'}};
	}
	my $snpcounter = 0;
	my $indelcounter = 0;
	my @positions = ();
	my @origins = ();
	my @refcalls = ();
	my @donorcalls = ();
	my @read_nucs = ();
	my @read_quals = ();
	my @read_flags = ();
	my @read_dinucs = ();
	my @read_readpos = ();
	my @snp_reads = ();
	my $donorcount = 0;
	my $refcount = 0;
	foreach my $key (sort {$a <=> $b} keys %fusion) {
		if (exists $snplist{$key}) {
			push (@positions, $key);
			if ($fusion{$key}{'nuc'} eq $snplist{$key}{'donor'}) {
				$donorcount += $fusion{$key}{'count'};
				push (@origins, 'donor');
				push (@refcalls, $snplist{$key}{'ref'});
				push (@donorcalls, $snplist{$key}{'donor'});
			} elsif ($fusion{$key}{'nuc'} eq $snplist{$key}{'ref'}) {
				$refcount += $fusion{$key}{'count'};
				push (@origins, 'ref');
				push (@refcalls, $snplist{$key}{'ref'});
				push (@donorcalls, $snplist{$key}{'donor'});
			} else {
				push (@origins, 'error');
				push (@refcalls, $snplist{$key}{'ref'});
				push (@donorcalls, $snplist{$key}{'donor'});
				if ($fusion{$key}{'nuc'} =~ m/\+|0/) {
					$indelcounter ++;
				} else {
					$snpcounter ++;
				}
			}
			push (@read_nucs, $fusion{$key}{'nuc'});
			push (@read_quals, $fusion{$key}{'qual'});
			push (@read_flags, $fusion{$key}{'flag'});
			push (@read_dinucs, $fusion{$key}{'dinuc'});
			push (@read_readpos, $fusion{$key}{'readpos'});
#			push (@snp_reads, $donorlist{$key}{'reads'});
		} else {
			if ($fusion{$key}{'md'} eq '1') {
				if ($fusion{$key}{'nuc'} =~ m/\+/) {
					$indelcounter ++;
				}
			} elsif ($fusion{$key}{'nuc'} eq '0') {
				$indelcounter ++;
			} else {
				$snpcounter ++;
			}
		}
	}
	my $callcount = $refcount + $donorcount;
	if ($callcount < 4) {
		print "fewSNPs\n";
		$loopcounter ++;
		next;
	} else {
		print OUT join("/", @positions);
		print OUT "\t";
		print OUT join("/", @origins);
		print OUT "\t";
		print OUT join("/", @refcalls);
		print OUT "\t";
		print OUT join("/", @donorcalls);
		print OUT "\t";
		print OUT join("/", @read_nucs);
		print OUT "\t";
		print OUT join("/", @read_quals);
		print OUT "\t";
		print OUT join("/", @read_dinucs);
		print OUT "\t";
		print OUT join("/", @read_flags);
		print OUT "\t";
		print OUT join("/", @read_readpos);
		print OUT "\t";
#		print OUT join("/", @snp_reads);
#		print OUT "\t";
		print OUT "$snpcounter\t$indelcounter\t$title1\n";
		print "success!\n";
		$loopcounter ++;
	}
}
