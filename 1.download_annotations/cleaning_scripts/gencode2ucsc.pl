#!/usr/bin/perl
my %h;
open(IFP, $ARGV[0]) or die;
while(<IFP>) {
  chomp $_;
  my @F = split(/\t/,$_);
  $h{$F[0]} = $F[1];
}
close IFP;

while(<STDIN>) {
  chomp $_;
  my @F = split(/\t/,$_);
  if($_ =~ m/^#/) { print "$_\n"; next; }
  die if(!defined($h{$F[0]}));
  $F[0] = $h{$F[0]};
  print join("\t",@F) . "\n";
}
