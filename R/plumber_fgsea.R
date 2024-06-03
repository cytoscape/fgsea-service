#* @apiTitle EnrichmentMap Web Service for computing enriched pathways.
#* @apiDescription Computes enriched pathways from RNA-Seq data or a pre-ranked gene list using baderlab gene set databases.

library(plumber)
library(fgsea)
library(edgeR)
library(here)

# Load pathway database into memory.
#gmt.file = here("data", "Human_GOBP_AllPathways_no_GO_iea_June_01_2022_symbol.gmt")
gmt.file = here("data", "Human_GOBP_AllPathways_noPFOCR_no_GO_iea_May_01_2024_symbol.gmt")
pathways <- gmtPathways(gmt.file)

# Ignore gene sets that are smaller or larger than the limits.
fgsea.minSize = 15
fgsea.maxSize = 500

# Columns from FGSEA results to return
result.cols = c("pathway", "size", "pval", "padj", "NES")


asMatrix <- function(x) {
  y <- as.matrix.data.frame(x[,-1])
  rownames(y) <- x[[1]]
  y
}


calculateRanks <- function(RNASeq, classes) {
  # Convert to matrix
  RNASeq <- asMatrix(RNASeq)
  
  # Create data structure to hold counts and subtype information for each sample.
  d <- DGEList(counts=RNASeq, group=classes)
  
  # Filter out genes with low read counts, as these are just noise
  keep <- filterByExpr(d)
  d <- d[keep, , keep.lib.sizes=FALSE]
  
  # Normalize the data
  d <- calcNormFactors(d)
  design <- model.matrix(~0 + classes)
  d <- estimateDisp(d, design)
  
  # Calculate differential expression statistics with a simple design.
  # classesB is the control
  my.contrasts <- makeContrasts(AvsB=classesA-classesB, levels=design)
  fit <- glmQLFit(d, design)
  qlf <- glmQLFTest(fit, coef=2, contrast=my.contrasts)
  
  table_res = topTags(qlf, n = nrow(d))
  
  ranks = sign(table_res$table$logFC) * -log10(table_res$table$PValue)
  
  # Put ranks into format for FGSEA
  genenames <- rownames(table_res$table)
  ranks <- data.frame(genenames, ranks)
  colnames(ranks) <- c("gene","rank")
  ranks <- setNames(ranks$rank, ranks$gene)
  
  return(ranks)
}


runFgseaRnaseq <- function(RNASeq, classes) {
  ranks <- calculateRanks(RNASeq, classes)
  finiteRanks <- ranks[is.finite(ranks)]
  
  # Run FGSEA
  fgseaRes <- fgsea(
    pathways, 
    finiteRanks, 
    minSize = fgsea.minSize, 
    maxSize = fgsea.maxSize
  )
  
  # remove unwanted columns
  fgseaRes <- fgseaRes[, ..result.cols]
  
  result <- list(
    ranks = as.list(finiteRanks), 
    pathways = fgseaRes, 
    messages = list()
  )
  
  # check if there were any non-finite ranks
  diff <- length(ranks) - length(finiteRanks)
  if(diff > 0) {
    result$messages[[1]] <- list(
      level = 'warning',
      type = 'non_finite_ranks',
      text = paste('there are', diff, 'non-finite ranks out of', length(ranks), 'total ranks'),
      data = list(totalRanks = length(ranks), finiteRanks = length(finiteRanks))
    )
  }
  
  result
}


runFgseaPreranked <- function(ranks) {
  fgseaRes <- fgsea(
    pathways, 
    ranks, 
    minSize = fgsea.minSize, 
    maxSize = fgsea.maxSize
  )
  
  # remove unwanted columns
  fgseaRes <- fgseaRes[, ..result.cols]
  
  list(pathways=fgseaRes)
}


# Used to ping this service to see if it responds
#* @get /
#* @serializer unboxedJSON
function(req) {
  list(ok=TRUE)
}


# curl --data-binary @brca_hd_tep_ranks_100.rnk -X POST "http://127.0.0.1:9404/v1/preranked" -H "Content-Type: text/tab-separated-values"
#* @post /v1/preranked
#* @parser tsv
#* @parser csv
#* @serializer unboxedJSON
function(req) {
  ranks <- req$body
  colnames(ranks) <- c("gene","rank")
  ranks <- setNames(ranks$rank, ranks$gene)
  
  runFgseaPreranked(ranks)
}


# curl --data-binary @FakeExpression.txt -X POST "http://127.0.0.1:3723/v1/rnaseq?classes=A,A,A,B,B,B" -H "Content-Type: text/tab-separated-values"
#* @post /v1/rnaseq
#* @parser tsv
#* @parser csv
#* @serializer unboxedJSON
function(req, classes) {
  RNASeq <- req$body
  
  # Parse the 'classes' query parameter into a vector
  classes <- strsplit(classes, split=",",fixed=TRUE)[[1]]
  
  # Boolean vector indicating which columns to keep
  columns.keep <- classes %in% c('A','B')
  
  # Drop columns that are to be ignored
  classes <- classes[columns.keep]
  RNASeq <- RNASeq[, c(TRUE, columns.keep)]
  
  runFgseaRnaseq(RNASeq, classes)
}
