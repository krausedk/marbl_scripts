#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

open(STDERR, '>', 'marbl_v3_errorlog.txt') or die "Can't open log";
my $loopcounter = 1;
my $recocounter = 0;
my $covcounter = 0;
my %marbl_coverage = ();
my %MARBL = ();
my %likelihood_map = ();
my $likelihood_file = 'deepb1_gatk_recalibrated_interpolated_capped_error_rates.txt';
open my $stuff, $likelihood_file or die "Could not open $likelihood_file: $!";
while(defined(my $line2 = <$stuff>)) {
	chomp($line2);
	my @fields = split(/\t/, $line2);
	$likelihood_map{$fields[0]}{$fields[1]}{$fields[2]} = $fields[3];
}
my $file = $ARGV[0];
open my $info, $file or die "Could not open $file: $!";
while(defined(my $line = <$info>)) {
	chomp($line);
	my @fields = split(/\t/, $line);
	my @positions = split(/\//, $fields[0]);
	my @origins = split(/\//, $fields[1]);
	my @nucs = split(/\//, $fields[2]);
	my @altnucs = split(/\//, $fields[3]);
	my @quals = split(/\//, $fields[4]);
	my @flags = split(/\//, $fields[5]);
	my @read_positions = split(/\//, $fields[6]);
	my $read_snps = $fields[7];
	my $read_indels = $fields[8];
	my $title = $fields[9];
	#same as earlier versions - store the likelihoods for assessing junctions#
	my @likelihoods;
	my @null_likelihoods_R2D;
	my @null_likelihoods_D2R;
	for my $i (0..$#quals) {
		if ($quals[$i] =~ m/,/) {
			my @sub_qual = split(/,/, $quals[$i]);
			my @sub_flag = split(/,/, $flags[$i]);
			my @sub_nuc = split(/,/, $nucs[$i]);
			my @sub_altnuc = split(/,/, $altnucs[$i]);
			my @sub_read_position = split(/,/, $read_positions[$i]);
			my $null_likelihood_R2D = 1;
			my $null_likelihood_D2R = 1;
			my $likelihood = 1;
			for my $j(0..$#sub_qual) {
				if ($sub_flag[$j] & 16) {
					# fix the nucs field so that it represents a reverse complement read, i.e. AATG becomes AAAC #
					my @sub_sub_nuc = split(//, $sub_nuc[$j]);
					my @sub_sub_altnuc = split(//, $sub_altnuc[$j]);
					for my $q (2..3) {
						$sub_sub_nuc[$q] =~ tr/ATCG/TAGC/;
						$sub_sub_altnuc[$q] =~ tr/ATCG/TAGC/;
					}
					$sub_nuc[$j] = join('', @sub_sub_nuc);
					$sub_altnuc[$j] = join('', @sub_sub_altnuc);
					if ($nucs[$i] =~ m/\+|0/) {
						# insertions/deletions get a 0.001 value #
						$likelihood *= .001;
						$null_likelihood_R2D *= .001;
						$null_likelihood_D2R *= .001;
					} else {
						if($likelihood_map{$sub_nuc[$j]}{$sub_qual[$j]}{$sub_read_position[$j]} lt $ARGV[1]) {
							$null_likelihood_R2D *= $likelihood_map{$sub_nuc[$j]}{$sub_qual[$j]}{$sub_read_position[$j]};
							$null_likelihood_D2R *= $likelihood_map{$sub_altnuc[$j]}{$sub_qual[$j]}{$sub_read_position[$j]};
						}
						if ($origins[$i] eq 'donor') {
							if ($likelihood_map{$sub_nuc[$j]}{$sub_qual[$j]}{$sub_read_position[$j]} lt $ARGV[1]) {
								$likelihood *= $likelihood_map{$sub_nuc[$j]}{$sub_qual[$j]}{$sub_read_position[$j]};
							}
						} else {
							if ($likelihood_map{$sub_altnuc[$j]}{$sub_qual[$j]}{$sub_read_position[$j]} lt $ARGV[1]) {
								$likelihood *= $likelihood_map{$sub_altnuc[$j]}{$sub_qual[$j]}{$sub_read_position[$j]};
							}
						}
					}
				} else {
					if ($nucs[$i] =~ m/\+|0/) {
						$likelihood *= .001;
						$null_likelihood_R2D *= .001;
						$null_likelihood_D2R *= .001;
					} else {
						if($likelihood_map{$sub_nuc[$j]}{$sub_qual[$j]}{$sub_read_position[$j]} lt $ARGV[1]) {
							$null_likelihood_R2D *= $likelihood_map{$sub_nuc[$j]}{$sub_qual[$j]}{$sub_read_position[$j]};
							$null_likelihood_D2R *= $likelihood_map{$sub_altnuc[$j]}{$sub_qual[$j]}{$sub_read_position[$j]};
						}
						if ($origins[$i] eq 'donor') {
							if($likelihood_map{$sub_nuc[$j]}{$sub_qual[$j]}{$sub_read_position[$j]} lt $ARGV[1]) {
								$likelihood *= $likelihood_map{$sub_nuc[$j]}{$sub_qual[$j]}{$sub_read_position[$j]};
							}
						} else {
							if($likelihood_map{$sub_altnuc[$j]}{$sub_qual[$j]}{$sub_read_position[$j]} lt $ARGV[1]) {
								$likelihood *= $likelihood_map{$sub_altnuc[$j]}{$sub_qual[$j]}{$sub_read_position[$j]};
							}
						}
					}
				}
			}
			push (@null_likelihoods_R2D, $null_likelihood_R2D);
			push(@null_likelihoods_D2R, $null_likelihood_D2R);
			push (@likelihoods, $likelihood);
		} else {
			my $likelihood = 1;
			my $null_likelihood_R2D = 1;
			my $null_likelihood_D2R = 1;
			if ($flags[$i] & 16) {
				# fix the nucs field so that it represents a reverse complement read, i.e. AATG becomes AAAC #
					my @sub_sub_nuc = split(//, $nucs[$i]);
					my @sub_sub_altnuc = split(//, $altnucs[$i]);
					for my $q (2..3) {
						$sub_sub_nuc[$q] =~ tr/ATCG/TAGC/;
						$sub_sub_altnuc[$q] =~ tr/ATCG/TAGC/;
					}
					$nucs[$i] = join('', @sub_sub_nuc);
					$altnucs[$i] = join('', @sub_sub_altnuc);
				if ($nucs[$i] =~ m/\+|0/) {
					$likelihood *= .001;
					$null_likelihood_R2D *= .001;
					$null_likelihood_D2R *= .001;
				} else {
					if($likelihood_map{$nucs[$i]}{$quals[$i]}{$read_positions[$i]} lt $ARGV[1]) {
						$null_likelihood_R2D *= $likelihood_map{$nucs[$i]}{$quals[$i]}{$read_positions[$i]};
						$null_likelihood_D2R *= $likelihood_map{$altnucs[$i]}{$quals[$i]}{$read_positions[$i]};
					}
					if ($origins[$i] eq 'donor') {
						if ($likelihood_map{$nucs[$i]}{$quals[$i]}{$read_positions[$i]} lt $ARGV[1]) {
							$likelihood *= $likelihood_map{$nucs[$i]}{$quals[$i]}{$read_positions[$i]};
						}
					} else {
						if($likelihood_map{$altnucs[$i]}{$quals[$i]}{$read_positions[$i]} lt $ARGV[1]) {
							$likelihood *= $likelihood_map{$altnucs[$i]}{$quals[$i]}{$read_positions[$i]};
						}
					}
				}
			} else {
				if ($nucs[$i] =~ m/\+|0/) {
					$likelihood *= .001;
					$null_likelihood_R2D *= .001;
					$null_likelihood_D2R *= .001;
				} else {
					if($likelihood_map{$nucs[$i]}{$quals[$i]}{$read_positions[$i]} lt $ARGV[1]) {
						$null_likelihood_R2D *= $likelihood_map{$nucs[$i]}{$quals[$i]}{$read_positions[$i]};
						$null_likelihood_D2R *= $likelihood_map{$altnucs[$i]}{$quals[$i]}{$read_positions[$i]};
					}
					if ($origins[$i] eq 'donor') {
						if($likelihood_map{$nucs[$i]}{$quals[$i]}{$read_positions[$i]} lt $ARGV[1]) {
							$likelihood *= $likelihood_map{$nucs[$i]}{$quals[$i]}{$read_positions[$i]};
						}
					} else {
						if($likelihood_map{$altnucs[$i]}{$quals[$i]}{$read_positions[$i]} lt $ARGV[1]) {
							$likelihood *= $likelihood_map{$altnucs[$i]}{$quals[$i]}{$read_positions[$i]};
						}
					}
				}
			}
			push (@null_likelihoods_R2D, $null_likelihood_R2D);
			push (@null_likelihoods_D2R, $null_likelihood_D2R);
			push (@likelihoods, $likelihood);
		}
	}
	### delete all the errors from the list, since they don't affect anything ###
	my @origins_error;
	my @likelihoods_error;
	my @null_likelihoods_R2D_error;
	my @null_likelihoods_D2R_error;
	my @positions_error;
	foreach my $origin(@origins) {
		push(@origins_error, $origin);
	}
	foreach my $like(@likelihoods) {
		push(@likelihoods_error, $like);
	}
	foreach my $position(@positions) {
		push(@positions_error, $position);
	}
	foreach my $null_R2D(@null_likelihoods_R2D) {
		push(@null_likelihoods_R2D_error, $null_R2D);
	}
	foreach my $null_D2R(@null_likelihoods_D2R) {
		push(@null_likelihoods_D2R_error, $null_D2R);
	}
	my %readinfo = ();
	# store R2D as 1 and D2R as -1 so that we can phase_switch easily during max_theoretical_breakpoint calculation #
	for my $i(0..$#positions_error) {
		if (($origins_error[$i] eq 'donor') or ($origins_error[$i] eq 'ref')) {
			$readinfo{$positions_error[$i]} = {'ori' => $origins_error[$i], '1' => $null_likelihoods_R2D[$i], '-1' => $null_likelihoods_D2R[$i], 'pval' => $likelihoods_error[$i]};
		}
	}
	### calculate the numbers of tests and theoretical breakpoints ###
	my %tests = ();
	my %theo_breaks = ();
	# calculate the leftside=reference junctions #
	foreach my $pos (sort {$a <=> $b} keys %readinfo) {
		my $p_val_before = $readinfo{$pos}{'-1'};
		foreach my $before (keys %readinfo) {
			if ($before >= $pos) {
				next;
			} else {
				$p_val_before *= $readinfo{$before}{'-1'};
			}
		}
		if ($p_val_before < $ARGV[2]) {
			my $p_val_after = 1;
			foreach my $after (keys %readinfo) {
				if ($after <= $pos) {
					next;
				} else {
					$p_val_after *= $readinfo{$after}{'1'};
				}
			}
			# get immediate next key #
			my $adjacentpos;
			foreach my $after (sort {$a <=> $b} keys %readinfo) {
				if ($after <= $pos) {
					next;
				} else {
					$adjacentpos = $after;
					last;
				}
			}
			if ($p_val_after < $ARGV[2]) {
				$tests{$pos}{$adjacentpos}{'leftref'} = 1;
				$tests{$pos}{$adjacentpos}{'either'} = 1;
			}
		}
	}
	# calculate the leftside=donor junctions #
	foreach my $pos (sort {$a <=> $b} keys %readinfo) {
		my $p_val_before = $readinfo{$pos}{'1'};
		foreach my $before (keys %readinfo) {
			if ($before >= $pos) {
				next;
			} else {
				$p_val_before *= $readinfo{$before}{'1'};
			}
		}
		if ($p_val_before < $ARGV[2]) {
			my $p_val_after = 1;
			foreach my $after (keys %readinfo) {
				if ($after <= $pos) {
					next;
				} else {
					$p_val_after *= $readinfo{$after}{'-1'};
				}
			}
			# get immediate next key #
			my $adjacentpos;
			foreach my $after (sort {$a <=> $b} keys %readinfo) {
				if ($after <= $pos) {
					next;
				} else {
					$adjacentpos = $after;
					last;
				}
			}
			if ($p_val_after < $ARGV[2]) {
				$tests{$pos}{$adjacentpos}{'leftdonor'} = 1;
				if (exists $tests{$pos}{$adjacentpos}{'either'}) {
					$tests{$pos}{$adjacentpos}{'both'} = 1;
				} else {
					$tests{$pos}{$adjacentpos}{'either'} = 1;
				}
			}
		}
	}
	# look for max theoretical breakpoints starting with left_reference #
	my @keys = sort {$a <=> $b} keys %readinfo;
	my $i = 0;
	my $j = 0;
	my $k = 0;
	my $p_val_left;
	my $p_val_right;
	# start phase at -1, D2R #
	my $phase = -1;
	while($i < $#keys) {
		if ($k == 0) {
			$p_val_left = 1;
			$p_val_left *= $readinfo{$keys[$i]}{$phase};
#			print OUT2 "$p_val_left\n";
#			print OUT2 "LEFT\n";
		}
		$j++;
		$p_val_right = 1;
		while($p_val_left >= $ARGV[2]) {
			$i++;
			if (exists $keys[$i]) {
				$p_val_left *= $readinfo{$keys[$i]}{$phase};
				$j++;
			} else {
				last;
			}
#			print OUT2 "$p_val_left\n";
#			print OUT2 "LEFT\n";
		}
#		print OUT2 "LEFTDONE!\n";
		$phase *= -1;
		if (exists $keys[$j]) {
			$p_val_right *= $readinfo{$keys[$j]}{$phase};
#			print OUT2 "$p_val_right\n";
#			print OUT2 "RIGHT\n";
		} else {
#			print OUT2 "NORIGHT\n";
			last;
		}
		while($p_val_right >= $ARGV[2]) {
			$j++;
			if(exists $keys[$j]) {
				$p_val_right *= $readinfo{$keys[$j]}{$phase};
#				print OUT2 "$p_val_right\n";
#				print OUT2 "RIGHT\n";
			} else {
#				print OUT2 "NORIGHT\n";
				last;
			}
		}
#		print OUT2 "RIGHTDONE!\n";
		if (($p_val_left < $ARGV[2]) and ($p_val_right < $ARGV[2])) {
			$theo_breaks{$keys[$i]}{$keys[$i+1]}{'leftref'} = 1;
			$i = $j;
			$k ++;
			$p_val_left = $p_val_right;
		} else {
			last;
		}
	}
	# look for max theoretical breakpoints starting with left_donor #
	$i = 0;
	$j = 0;
	$k = 0;	
	# start phase at 1, R2D #
	$phase = 1;
	while($i < $#keys) {
		if ($k == 0) {
			$p_val_left = 1;
			$p_val_left *= $readinfo{$keys[$i]}{$phase};
		}
		$j++;
		$p_val_right = 1;
		while($p_val_left >= $ARGV[2]) {
			$i++;
			if(exists $keys[$i]) {
				$p_val_left *= $readinfo{$keys[$i]}{$phase};
				$j++;
			} else {
				last;
			}
		}
		$phase *= -1;
		if (exists $keys[$j]) {
			$p_val_right *= $readinfo{$keys[$j]}{$phase};
		} else {
			last;
		}
		while($p_val_right >= $ARGV[2]) {
			$j++;
			if(exists $keys[$j]) {
				$p_val_right *= $readinfo{$keys[$j]}{$phase};
			} else {
				last;
			}
		}
		if (($p_val_left < $ARGV[2]) and ($p_val_right < $ARGV[2])) {
			$theo_breaks{$keys[$i]}{$keys[$i+1]}{'leftdonor'} = 1;
			$i = $j;
			$k ++;
			$p_val_left = $p_val_right;
		} else {
			last;
		}
	}
	my $theo_both_count = 0;
	my $theo_leftref_count = 0;
	my $theo_leftdonor_count = 0;
	my $theo_either_count = 0;
	my $theo_break_leftref = 0;
	my $theo_break_leftdonor =  0;
	foreach my $pos1 (keys %tests) {
		foreach my $pos2 (keys %{ $tests{$pos1} }) {
			$theo_both_count++ if exists $tests{$pos1}{$pos2}{'both'};
			$theo_leftref_count++ if exists $tests{$pos1}{$pos2}{'leftref'};
			$theo_leftdonor_count++ if exists $tests{$pos1}{$pos2}{'leftdonor'};
			$theo_either_count++;
		}
	}
	foreach my $pos1 (keys %theo_breaks) {
		foreach my $pos2 (keys %{ $theo_breaks{$pos1} }) {
			$theo_break_leftref++ if exists $theo_breaks{$pos1}{$pos2}{'leftref'};
			$theo_break_leftdonor++ if exists $theo_breaks{$pos1}{$pos2}{'leftdonor'};
		}
	}
	### finally, run marbl and look for breakpoints ###
	$i = 0;
	$j = 0;
	$k = 0;
	my %marbl_breakpoints = ();
	my %origin_tracker = ();
	my $break = 0;
	my $p_chi_left = 1;
	my $p_chi_right = 1;
	my $null_left = 1;
	my $null_right = 1;
	while($i < $#keys) {
		if ($readinfo{$keys[$i]}{'ori'} ne $readinfo{$keys[$i+1]}{'ori'}) {
			my $left_origin = $readinfo{$keys[$i]}{'ori'};
			my $right_origin = $readinfo{$keys[$i+1]}{'ori'};
			$p_chi_left *= $readinfo{$keys[$i]}{'pval'};
			$p_chi_right *= $readinfo{$keys[$i+1]}{'pval'};
			if($p_chi_left < $ARGV[2]) {
				if($p_chi_right < $ARGV[2]) {
					$marbl_breakpoints{$keys[$i]}{$keys[$i+1]} ++;
					$marbl_coverage{$keys[$i]}{$keys[$i+1]} ++;
					$break ++;
					if ($readinfo{$keys[$i]}{'ori'} eq 'ref') {
						$origin_tracker{$keys[$i]}{$keys[$i+1]} ++;
					}
					$p_chi_left = $p_chi_right;
					$p_chi_right = 1;
					$p_chi_left /= $readinfo{$keys[$i+1]}{'pval'};
					$i ++;
					next;
				} else {
					$k = $i + 1;
					while ($k < $#keys) {
						$k++;
						if($right_origin ne $readinfo{$keys[$k]}{'ori'}) {
							last;
						} else {
							$p_chi_right *= $readinfo{$keys[$k]}{'pval'};
						}
					}
					if($p_chi_right < $ARGV[2]) {
						$marbl_breakpoints{$keys[$i]}{$keys[$i+1]} ++;
						$marbl_coverage{$keys[$i]}{$keys[$i+1]} ++;
						$break ++;
						if($readinfo{$keys[$i]}{'ori'} eq 'ref') {
							$origin_tracker{$keys[$i]}{$keys[$i+1]} ++;
						}
						$p_chi_left = $p_chi_right;
						$p_chi_left /= $readinfo{$keys[$k]}{'pval'};
						$p_chi_right = 1;
						$i = $k;
						next;
					} else {
						$p_chi_left = 1;
						$p_chi_right = 1;
						$i = $k;
						next;
					}
				}
			} else {
				$p_chi_left = 1;
				$p_chi_right = 1;
				$i++;
				next;
			}
		} else {
			$p_chi_left *= $readinfo{$keys[$i]}{'pval'};
			$i++;
			next;
		}
	}
	# given the breakpoints identified, recalculate places where breaks can happen #
	# this step is in progress unfortunately #
	# future work on this will possibly be necessary, although it is a complicated problem that is relatively negligible #
	


#	if ($break == 0) {
#		my $theoretical_max_breaks = $theo_break_leftref;
#		$theoretical_max_breaks = $theo_break_leftdonor if $theo_break_leftdonor > $theo_break_leftref;
#		foreach my $pos1 (keys %tests) {
#			foreach my $pos2 (keys %{ $tests{$pos1} }) {
#				$MARBL{$pos1}{$pos2}{'theory'}{'both'} += $theoretical_max_breaks/$theo_either_count if exists $tests{$pos1}{$pos2}{'both'};
#				$MARBL{$pos1}{$pos2}{'theory'}{'leftref'} += $theoretical_max_breaks/$theo_either_count if exists $tests{$pos1}{$pos2}{'leftref'};
#				$MARBL{$pos1}{$pos2}{'theory'}{'leftdonor'} += $theoretical_max_breaks/$theo_either_count if exists $tests{$pos1}{$pos2}{'leftdonor'};
#			}
#		}
#	} else {
	my $theoretical_max_breaks = $theo_break_leftref;
	$theoretical_max_breaks = $theo_break_leftdonor if $theo_break_leftdonor > $theo_break_leftref;
	foreach my $pos1 (keys %tests) {
		foreach my $pos2 (keys %{ $tests{$pos1} }) {
			$MARBL{$pos1}{$pos2}{'theory'}{'both'} += $theoretical_max_breaks/$theo_either_count if exists $tests{$pos1}{$pos2}{'both'};
			$MARBL{$pos1}{$pos2}{'theory'}{'leftref'} += $theoretical_max_breaks/$theo_either_count if exists $tests{$pos1}{$pos2}{'leftref'};
			$MARBL{$pos1}{$pos2}{'theory'}{'leftdonor'} += $theoretical_max_breaks/$theo_either_count if exists $tests{$pos1}{$pos2}{'leftdonor'};
		}
	}
	foreach my $pos1 (keys %marbl_breakpoints) {
		foreach my $pos2 (keys %{ $marbl_breakpoints{$pos1} }) {
			$MARBL{$pos1}{$pos2}{'practice'} ++;
		}
	}
	foreach my $pos1 (keys %origin_tracker) {
		foreach my $pos2 (keys %{ $origin_tracker{$pos1} }) {
			$MARBL{$pos1}{$pos2}{'breakleftref'} ++;
		}
	}
	print "$loopcounter\n";
	$loopcounter ++;
}
open(OUT, ">$ARGV[0]_$ARGV[2]_marbloutput.v3.out");
print OUT "LEFT\tRIGHT\tBOTH\tLEFTREF\tLEFTDONOR\tBREAKS\tBREAKSREFLEFT\n";
foreach my $l (sort {$a <=> $b} keys %MARBL) {
	foreach my $r (sort {$a <=> $b} keys %{ $MARBL{$l} }) {
		print OUT "$l\t$r\t";
		if(exists $MARBL{$l}{$r}{'theory'}{'both'}) {
			print OUT "$MARBL{$l}{$r}{'theory'}{'both'}\t";
		} else {
			print OUT "0\t";
		}
		if(exists $MARBL{$l}{$r}{'theory'}{'leftref'}) {
			print OUT "$MARBL{$l}{$r}{'theory'}{'leftref'}\t";
		} else {
			print OUT "0\t";
		}
		if(exists $MARBL{$l}{$r}{'theory'}{'leftdonor'}) {
			print OUT "$MARBL{$l}{$r}{'theory'}{'leftdonor'}\t";
		} else {
			print OUT "0\t";
		}
		if(exists $MARBL{$l}{$r}{'practice'}) {
			print OUT "$MARBL{$l}{$r}{'practice'}\t";
		} else {
			print OUT "0\t";
		}
		if(exists $MARBL{$l}{$r}{'breakleftref'}) {
			print OUT "$MARBL{$l}{$r}{'breakleftref'}\n";
		} else {
			print OUT "0\n";
		}
	}
}
