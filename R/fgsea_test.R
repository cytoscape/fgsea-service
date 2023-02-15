# FGSEA tutorial
# http://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
# Manual
# https://bioconductor.org/packages/release/bioc/manuals/fgsea/man/fgsea.pdf
# Paper
# biorxiv.org/content/10.1101/060012v3.full

# Docs for read.gmt()
# https://www.rdocumentation.org/packages/qusage/versions/2.4.0/topics/read.gmt

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("fgsea")


source('R/plumber_fgsea.R')


writeFgseaRes <- function(fgseaRes, fileName="enrichments_results.txt") {
  write.table(fgseaRes, 
              file=paste(here(), "/data/", fileName, sep=""),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
}


testRanks1 <- function() {
  rnk.file <- here("data", "brca_hd_tep_ranks.rnk")
  
  ranks <- read.table(rnk.file, header=TRUE, colClasses = c("character", "numeric"))
  colnames(ranks) <- c("g","r")
  ranks <- setNames(ranks$r, ranks$g)
  
  runFgseaPreranked(ranks)
}


testRnaSeq1 <- function() {
  exp.file <- here("data", "FakeExpression.txt")
  
  # Fails here if the table can't be read...
  colClasses <- c("character", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric")
  RNASeq <- read.table(exp.file, header=TRUE, colClasses = colClasses)
  RNASeq <- RNASeq[ , !(tolower(names(RNASeq)) %in% c("description"))] # remove description column
  
  classes <- c('A', 'A', 'A', 'B', 'B', 'B')
  
  runFgseaRnaseq(RNASeq, classes)
}




