#!/usr/bin/perl
#./clean_sam.pl output.lim offtarget_limit
# takes SAM output, with MD tag from batmis --misplus option (on read strand),
# tells you for each query if it hit an output reporting limit,
# filters non-NGG hits, and
# produces "MM" tag, which is an int array of mismatch locations
open(LIM,">$ARGV[0]") or die;
my ($r, $c) = ("", 0);
while(<STDIN>)
{
  if($_ !~ m/^([^@]\S+)/) { print $_; next; }
  if($1 ne $r) {
    print LIM ($c >= $ARGV[1] ? "$r\t1\n" : "$r\t0\n") if($r ne "");
    ($r, $c) = ($1, 0);
  }
  $c++;
  
  die if($_ !~ m/MD:Z:(\S+)\n$/);
  my $md = $1;
  next if($md =~ m/\D1?$/);
  $md =~ s/\D2$//;
  my @F = split(/\D/,$md);
  my @R;
  my $i = 0;
  foreach(@F[0..($#F -1)]) { $i += $_ +1; push(@R,$i); }
  next if($#R > 3);
  chomp $_;
  print $_;
  print "\t" . join(",","MM:B:S",@R) if($#R > -1);
  print "\n";
}
print LIM ($c >= $ARGV[1] ? "$r\t1\n" : "$r\t0\n") if($r ne "");
close LIM;
