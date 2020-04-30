#!/usr/bin/perl
use strict;
use warnings;

open(SAM,">$ARGV[1]") or die;
open(AZIMUTH,">$ARGV[2]") or die;
open(FASTA,">$ARGV[3]") or die;
sub rc { my $r = reverse($_[0]); $r =~ tr/ACGT/TGCA/; return $r; }
sub guidepull {
  my ($chr, $seq, $offset) = @_;
  while($seq =~ /(?=([ACGT]{4}([ACGT]{20})([ACGT]GG)[ACGT]{3}))/g) { #NGG guides
    my $p = $-[2] + $offset +1; #1-indexed 20nt pos
    my $k = "$chr:+:" . scalar($-[2] + $offset +1); #1-indexed 23nt pos
    my $c = $p +17; #1-indexed cut site
    print SAM "$k\t0\t$chr\t$p\t0\t20M\t*\t0\t0\t$2\t*\tCS:i:$c\tAM:Z:$3\n";
    print AZIMUTH "$k\t$1\n";
    print FASTA ">$k\n$2$3\n";
  }
  while($seq =~ /(?=([ACGT]{3}(CC[ACGT])([ACGT]{20})[ACGT]{4}))/g) { #CCN guides
    my $p = $-[3] + $offset +1; #1-indexed 20nt pos
    my $k = "$chr:-:" . scalar($-[2] + $offset +1); #1-indexed 23nt pos
    my $c = $p +3; #1-indexed cut site
    my $az = rc($1);
    my $pam = rc($2);
    my $guide = rc($3);
    print SAM "$k\t16\t$chr\t$p\t0\t20M\t*\t0\t0\t$3\t*\tCS:i:$c\tAM:Z:$pam\n";
    print AZIMUTH "$k\t$az\n";
    print FASTA ">$k\n$guide$pam\n";
  }
}

my ($chr, $seq, $offset) = ("", "", 0);
open(IFP,$ARGV[0]) or die;
while(<IFP>) {
  chomp $_;
  if($_ !~ m/^>(\S+)/) { $seq .= uc($_); next; } #sequence row
  if($chr ne "") { #first chr
    ($chr, $offset) = split(/:/,$chr) if($chr =~ m/:/);
    guidepull($chr, $seq, $offset);
  }
  ($chr, $seq, $offset) = ($1, "", 0);
}
($chr, $offset) = split(/:/,$chr) if($chr =~ m/:/);
guidepull($chr, $seq, $offset);
close IFP;
