## Data Availability
https://github.com/veeneman/PINCER/releases/tag/1.1

## How to load guide database in R
```R
#Run this code once to install libraries
install.packages("BiocManager");
library(BiocManager);
install("GenomicRanges");
install("Biostrings");

#Load libraries
library("GenomicRanges");
library("Biostrings");

#Load guide database
load("database.hg38.Rdata");

#Look at guide database
head(guides,2);
```
