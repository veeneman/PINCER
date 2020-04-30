#!/usr/bin/perl

while(<STDIN>)
{
  chomp $_;
  while(($_ =~ m/^\s+Site\s+order/ || $_ =~ m/^\s+misc_feature\s+order/) && $_ !~ m/\)\s*$/)
  {
   die if($_ =~ m/\s$/);
   my $tmp = <STDIN>;
   chomp $tmp;
   die if($tmp =~ m/\s$/);
   $tmp =~ s/^\s+//;
   $_ = "$_$tmp";
  }
  print "$_\n";
}
