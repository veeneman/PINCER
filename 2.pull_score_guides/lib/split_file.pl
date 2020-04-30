#!/usr/bin/perl
#./split_file.pl file suffix chunk_count unit_lines
#gives [1-#chunks]suffix files, for use with LSF queue, keeping unit_lines together
my ($file, $suffix, $chunk_ct, $unit_lines) = @ARGV;

open(IFP, $file) or die;
my $line_ct = `wc -l < $file`;
die if($line_ct % $unit_lines != 0);
my $units = $line_ct / $unit_lines;

my $chunk_min = int($units / $chunk_ct);
my $chunk_rem = $units % $chunk_ct;
for(my $chunk = 1; $chunk <= $chunk_ct; $chunk++)
{
  open(OFP,">$chunk$suffix") or die;
  my $chunk_size = $chunk_min + ($chunk <= $chunk_rem ? 1 : 0);
  for(my $i = 0; $i < $chunk_size; $i++) {
    for(my $j = 0; $j < $unit_lines; $j++) {
      $_ = <IFP>; print OFP $_;
    }
  }
  close OFP;
}
