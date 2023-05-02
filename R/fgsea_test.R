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


readExpFile <- function(fileName, colClasses, colsToDrop=c()) {
  RNASeq <- read.table(file = here("data", fileName),
                       colClasses = colClasses, 
                       header = TRUE, 
                       sep = '\t')
  RNASeq[, !(names(RNASeq) %in% colsToDrop)] 
}

readExpFileAndRunFGSEA <- function(fileName, colClasses, colsToDrop=c(), classes) {
  RNASeq <- readExpFile(fileName, colClasses, colsToDrop)
  runFgseaRnaseq(RNASeq, classes)
}

writeFgseaRes <- function(fgseaRes, fileName="enrichments_results.txt") {
  write.table(fgseaRes, 
              file=paste(here(), "/data/", fileName, sep=""),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
}


testPreranked1 <- function() {
  rnk.file <- here("data", "brca_hd_tep_ranks.rnk")
  
  ranks <- read.table(rnk.file, header=TRUE, colClasses = c("character", "numeric"))
  colnames(ranks) <- c("g","r")
  ranks <- setNames(ranks$r, ranks$g)
  
  runFgseaPreranked(ranks)
}

testCalculateRanks1 <- function() {
  RNASeq <- readExpFile(
    fileName = "GSE129943_rsem_counts_2016.txt", 
    colClasses = c("character", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"),
    colsToDrop = c("HGNC")
  )
  classes <- c('A', 'A', 'A', 'B', 'B', 'B')
  calculateRanks(RNASeq, classes)  
}

testRnaSeq1 <- function() {
  readExpFileAndRunFGSEA(
    fileName = "FakeExpression.txt", 
    colClasses = c("character", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"),
    colsToDrop = c("description"),
    classes = c('A', 'A', 'A', 'B', 'B', 'B')
  )
}


testRnaSeq2 <- function() {
  readExpFileAndRunFGSEA(
    fileName = "GSE129943_rsem_counts_2016.txt", 
    colClasses = c("character", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"),
    colsToDrop = c("HGNC"),
    classes = c('A', 'A', 'A', 'B', 'B', 'B')
  )
}



