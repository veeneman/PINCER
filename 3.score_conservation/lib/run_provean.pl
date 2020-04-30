#!/usr/bin/perl
use IO::Handle;

my $NP = $ARGV[0];

open($log,">$NP.log") or die;
$log->autoflush;
print $log "$NP Start\t" . `date`;

if(! -e "$NP.fa") { print STDERR "$NP.fa does not exist\n"; exit 0; }
my $seq = "";
open(IFP, "$NP.fa") or die;
my $h = <IFP>; die $h if($h !~ m/^>/);
while(<IFP>) { die if(m/^>/); chomp $_; $seq .= $_; }
close IFP;

my @F = split(//, $seq);
open(VAR,">$NP.var") or die;
for($i = 0; $i <= $#F; $i++) { print VAR "$F[$i]" . scalar($i +1) . "del\n"; }
close VAR;

print $log "$NP Provean\t" . `date`;
system("provean.sh --num_threads 8 -q $NP.fa -v $NP.var > $NP.provean 2>> $NP.log2");

open(BED,">$NP.bed") or die;
open(IFP,"$NP.provean") or die;
while(<IFP>) { last if(m/^# VARIATION/); }
while(<IFP>)
{
  chomp;
  my @F = split(/\t/,$_);
  die if($F[0] !~ m/(\d+)del$/);
  print BED "$NP\t" . scalar($1 -1) . "\t$1\t$F[1]\n";
}
close BED;

unlink "$NP.fa", "$NP.var", "$NP.provean";
print $log "$NP Done\t" . `date`;
close $log;
