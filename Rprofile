#set R terminal width to same as bash
WIDTH = system("stty size 2> /dev/null | awk '{print $2}'", intern=T);
if(length(WIDTH) > 0) {
  if(!is.na(suppressWarnings(as.integer(WIDTH)))) {
    options(width = as.integer(WIDTH));
  }
}
rm(WIDTH);

#set CRAN repository permanently
local({
  r <- getOption("repos")
  r["CRAN"] <- "http://ftp.ussg.iu.edu/CRAN/"
  options(repos = r)
})

#These are very helpful, but may damage code portability and are therefore disabled
##library wrapping function which suppresses messages and will load multiple libraries
library = function(...) {
  libs = str_replace_all(sapply(eval(substitute(alist(...))), deparse),"\"","")
  invisible(suppressMessages(lapply(libs,  base::library, character.only = T)));
}
#use cairo by default; X11 isn't installed here
if(capabilities("cairo")) {
  png = function(...) { grDevices::png(..., type = "cairo"); }
}
#libraries I use constantly
#suppressMessages(source("http://bioconductor.org/biocLite.R")); #slow; OLD
#BiocManager::install("Biostrings"); #current bioclite
options("stringsAsFactors" = FALSE);
base::library(stringr);
library(openxlsx);
#library(gdata);
