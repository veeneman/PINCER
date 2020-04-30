#!/usr/bin/env Rscript
#./process.3p.R refseq.clean.Rdata appris.clean.txt output.Rdata
library("stringr", "stringi", "rtracklayer", "GenomicAlignments", "GenomicFeatures");

argv = commandArgs(trailingOnly = TRUE);

load(argv[1]);
appris = read.table(argv[2], header=T);

a = refseq[refseq$type == "CDS" & refseq$Name %in% appris$NP, "Name"];
a = sortSeqlevels(a);
a = a[order(seqnames(a), start(a), end(a))];

b = split(a, paste(seqnames(a), a$Name));
c = unname(unlist(GRangesList(lapply(b, function(i) {
      o = if(as.character(strand(i[1])) == "+") { 1:length(i); } else { length(i):1; }
      i$offset = cumsum(c(0, width(i)[o]))[o];
      i$length = sum(width(i));
      return(i);
    }))));

d = slidingWindows(c, width = 0, step = 1); #slidingWindows(c, width = 1, step = 1)
mcols(d) = mcols(c);
l = lengths(d);
seq2 = Vectorize(seq.default, vectorize.args = c("from", "to"));
i = seq2(rep(0, length(l)), l -1); #seq2(rep(1, length(l)), l);
neg = as.logical(strand(c) == "-");
i[neg] = seq2(l[neg] -1, rep(0, length(l[neg]))); #seq2(l[neg], rep(1, length(l[neg])));
i = do.call("c", i);
e = as.data.frame(d, use.outer.mcols = TRUE);
e = GRanges(e$seqnames,
            IRanges(e$start, e$end),
            e$strand,
            DataFrame("Name" = e$Name,
                      "ThreePrime.nt" = e$offset + i,
                      "ThreePrime.pct" = (e$offset + i) / e$length));
e = sortSeqlevels(e);
e = e[order(seqnames(e), start(e), end(e))];

ThreePrime = e;
save(ThreePrime, file = argv[3])

#f = data.frame(e)[,c("seqnames", "start", "end", "ThreePrime.nt")]
#write.table(f, file = "debug.bedgraph", quote = F, sep = "\t", row.names = F, col.names = F);
