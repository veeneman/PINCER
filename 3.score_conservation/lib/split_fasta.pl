#!/usr/bin/perl
#./split_fasta.pl /home/veeneb/hpc/ref/hg38/anno/refseq.aa.gz NP.human.txt

my ($AA, $NP) = @ARGV;

my %h;
open(IFP, $NP) or die;
while(<IFP>) { chomp; $h{$_} = 1; }
close IFP;

open(OFP, ">/dev/null") or die;
open(IFP, "gunzip -c $AA |") or die;
while(<IFP>)
{
  if($_ =~ m/^>(\S+)/)
  {
    close(OFP);
    open(OFP, (defined($h{$1}) ? ">$1.fa" : ">/dev/null"));
  }
  print OFP $_;
}
