#!/usr/bin/perl
#ARGV is ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/GFF/ref_GRCh38.p7_top_level.gff3.gz
#GFF is ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/chr_accessions_GRCh38.p7
open(MAP,"$ARGV[0]") or die;
open(GFF,"zcat $ARGV[1] |") or die;

my %h;
while(<MAP>)
{
  next if($_ =~ m/^\#/);
  chomp $_;
  my @F = split(/\t/,$_); 
  $F[0] =~ s/^MT$/M/;
  $h{$F[1]} = "chr$F[0]";
}

my $header = 1;
while(<GFF>)
{
  if($_ =~ m/^\#/ && $header) { print $_; next; }
  elsif($_ =~ m/^\#/) { next; }
  $header = 0;
  chomp $_;
  my @F = split(/\t/,$_);
  #next if($F[2] eq "region");
  #next if($F[2] eq "match");
  #next if($F[2] eq "cDNA_match");
  #next if($F[2] eq "D_loop");
  #next if($F[2] eq "repeat_region");
  #next if($F[2] eq "C_gene_segment");
  #next if($F[2] eq "V_gene_segment");
  #next if($F[2] eq "D_gene_segment");
  #next if($F[2] eq "J_gene_segment");
  $F[0] = $h{$F[0]} if(defined($h{$F[0]}));
  print join("\t",@F) . "\n";
}

