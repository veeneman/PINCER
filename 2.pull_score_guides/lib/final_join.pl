#!/usr/bin/perl
#./final_join.pl sam_header sam azm lim hsu

#sam header
open(IFP, $ARGV[0]) or die;
while(<IFP>) { print $_; }
close IFP;

#sam: chr8:+:127734754  0 chr8  127734753 0 20M * 0 0 ACTCCCCCCCCCCCAAAAAA  * CS:i:127734770  AM:Z:AGG
#azm: chr8:+:127734754  3.495460661107916978e-01
#lim: chr8:+:127734754  1
#hsu: chr8:+:127734754  0.108163003873157 1 3 32  140 449 0.991140083942877 0 0 0 0 8
open(SAM, $ARGV[1]) or die;
open(AZM, $ARGV[2]) or die;
open(LIM, $ARGV[3]) or die;
open(HSU, $ARGV[4]) or die;
while(<SAM>)
{
  chomp $_;
  my ($k_sam, @F) = split(/\t/,$_);
  my $l_azm = <AZM>; die if(!defined($l_azm)); chomp $l_azm; my ($k_azm, $azm) = split(/\t/,$l_azm);
  my $l_lim = <LIM>; die if(!defined($l_lim)); chomp $l_lim; my ($k_lim, $lim) = split(/\t/,$l_lim);
  my $l_hsu = <HSU>; die if(!defined($l_hsu)); chomp $l_hsu; my ($k_hsu, @hsu) = split(/\t/,$l_hsu);
  die "ERROR: KEY MISMATCH" if($k_sam ne $k_azm || $k_sam ne $k_lim || $k_sam ne $k_hsu);
  print join("\t", "*", @F, "AZ:f:$azm", "OL:i:$lim", "HG:f:$hsu[0]", "HP:f:$hsu[6]",
             join(",","MG:B:S",@hsu[1..5]),join(",","MP:B:S",@hsu[7..11])) . "\n";
}
die if(!eof(AZM));
die if(!eof(LIM));
die if(!eof(HSU));
close SAM; close AZM; close LIM; close HSU;
