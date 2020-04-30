#!/usr/bin/env Rscript
#./score_hsu13.R genes.CDS.Rdata ot.cutoff input.bam
suppressMessages(library(GenomicAlignments));
library(stringr);
bam.params = ScanBamParam(tag="MM", what=c("qname", "flag", "rname", "strand", "pos"));

weights = c(0.000,0.000,0.014,0.000,0.000,
            0.395,0.317,0.000,0.389,0.079,
            0.445,0.508,0.613,0.851,0.732,
            0.828,0.615,0.804,0.685,0.583);
die = function(msg) {
  write(paste("Error:", msg), stderr());
  q(status = 1);
}
parse.args = function() {
  argv = commandArgs(trailingOnly = TRUE);
  
  if(length(argv) != 3) {
    write(paste(argv,collapse="\t"),stderr());
    die("Usage: score_hsu13.R genes.CDS.Rdata ot.cutoff input.bam");
  }
  if(!file.exists(argv[1])) {
    die(paste("CDS.Rdata not found:", argv[1]));
  }
  if(is.na(strtoi(argv[2]))) {
    die(paste("not an integer:", argv[2]));
  }
  if(!file.exists(argv[3])) {
    die(paste("input.bam not found:", argv[3]));
  }
  return(argv);
}
hsu = function(x) {
  lx = length(x);
  if(lx == 0) { return(c(100, 0)); }
  score = 100 * prod(1 - weights[x]);
  if(lx > 1) {
    mean_pairwise = sum(x[-1] - x[-lx]) / (lx -1);
    mpw_factor = (((19 - mean_pairwise) / 19) * 4) + 1;
    scl_factor = lx ^ 2;
    score = score / ( mpw_factor * scl_factor );
    score = max(score, 0);
  }
  return(c(score, lx));
}
hsu.agg = function(r) {
  score = 100 / (100 + sum(r[1,], na.rm = T));
  ot.ct = unname(table(factor(r[2,],levels = 0:4)));
  return(c(score, ot.ct));
}

argv = parse.args();
#argv = c("/hpc/grid/oncology/veeneb/ref/hg38/anno/gencode.CDS.Rdata", "5000", "1.bam")
load(argv[1]);
ot.cutoff = strtoi(argv[2]);
input.bam = argv[3];
IFP = open(BamFile(input.bam, obeyQname = TRUE, yieldSize=1)); #SPEED IMPROVEMENT: READ BATCHES
while(TRUE) {
  x = scanBam(IFP, param = bam.params)[[1]];
  x.nrow = length(x$qname);
  if(x.nrow == 0) { break; } #end of file
  if(is.null(x$tag$MM)) { x$tag$MM = vector("list", x.nrow); }
  
  #remove "unaligned" alignments (batmis bug)
  u.r = which(bitwAnd(x$flag, 4) != 0);
  if(length(u.r) != 0) {
    x$qname  = x$qname[-u.r];
    x$flag   = x$flag[-u.r];
    x$rname  = x$rname[-u.r];
    x$strand = x$strand[-u.r];
    x$pos    = x$pos[-u.r];
    x$tag$MM = x$tag$MM[-u.r];
    x.nrow = length(x$qname);
    if(x.nrow == 0) { next; }
  }
  
  #remove on-target hit, not guaranteed to be in output (if reporting cut off)
  o.str = x$qname[1]; #guide origin, eg: chr8:+:127734096
  if(x.nrow == 1) { #single alignment
    write(paste(c(o.str, rep(c(1, rep(0, 5)), 2)), collapse="\t"), stdout());
    next;
  }
  
  o = unlist(strsplit(o.str, ":"));
  o.r = which(x$rname == o[1] & x$strand == o[2] & x$pos == strtoi(o[3]));
  if(length(o.r) != 0) {
    x$qname  = x$qname[-o.r];
    x$flag   = x$flag[-o.r];
    x$rname  = x$rname[-o.r];
    x$strand = x$strand[-o.r];
    x$pos    = x$pos[-o.r];
    x$tag$MM = x$tag$MM[-o.r];
  };
  
  r = vapply(x$tag$MM, function(x) { hsu(x); }, numeric(2));
  r.wg = hsu.agg(r); #genome-wide hsu scores
  cut.pos = x$pos + c(17,6)[(x$strand == "+") +1];
  cut.gr = GRanges(seqnames = x$rname,ranges = IRanges(cut.pos, cut.pos), strand = "*");
  cds.hits = !is.na(findOverlaps(cut.gr, cds, select = "arbitrary", ignore.strand = TRUE));
  r.cds = hsu.agg(r[, cds.hits, drop=F]); #genome-wide hsu scores
  write(paste(c(o.str, r.wg, r.cds), collapse="\t"), stdout());
}
