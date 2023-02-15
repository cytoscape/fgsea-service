#* @apiTitle EnrichmentMap Web Service for computing enriched pathways.
#* @apiDescription Computes enriched pathways from RNA-Seq data or a pre-ranked gene list using baderlab gene set databases.

library(plumber)
library(fgsea)
library(edgeR)
library(here)

# Load pathway database into memory.
gmt.file = here("data", "Human_GOBP_AllPathways_no_GO_iea_June_01_2022_symbol.gmt")
pathways <- gmtPathways(gmt.file)

# Ignore gene sets that are smaller or larger than the limits.
fgsea.minSize = 15
fgsea.maxSize = 500

# Columns from FGSEA results to return
result.cols = c("pathway", "size", "pval", "padj", "ES", "NES")


as_matrix <- function(x) {
  y <- as.matrix.data.frame(x[,-1])
  rownames(y) <- x[[1]]
  y
}


runFgseaPreranked <- function(ranks) {
  fgseaRes <- fgsea(
    pathways, 
    ranks, 
    minSize = fgsea.minSize, 
    maxSize = fgsea.maxSize
  )
  
  res <- fgseaRes[, ..result.cols]
  list(pathways=res)
}


runFgseaRnaseq <- function(RNASeq, classes) {
  # Convert to matrix
  RNASeq <- as_matrix(RNASeq)
  
  # Create data structure to hold counts and subtype information for each sample.
  d <- DGEList(counts=RNASeq, group=classes)
  
  # Filter out genes with low read counts, as these are just noise
  keep <- filterByExpr(d)
  d <- d[keep, , keep.lib.sizes=FALSE]
  
  # Normalize the data
  d <- calcNormFactors(d)
  # Calculate dispersion
  d <- estimateCommonDisp(d)
  d <- estimateTagwiseDisp(d)
  
  # Calculate differential expression statistics with a simple design
  de <- exactTest(d, pair=c("A", "B"))
  tt_exact_test <- topTags(de, n=nrow(d))
  tt <- tt_exact_test
  
  # Calculate ranks
  ranks = sign(tt$table$logFC) * -log10(tt$table$PValue)
  
  # Put ranks into format for FGSEA
  genenames <- rownames(tt$table)
  ranks <- data.frame(genenames, ranks)
  colnames(ranks) <- c("gene","rank")
  ranks <- setNames(ranks$rank, ranks$gene)
  
  # Run FGSEA
  fgseaRes <- fgsea(
    pathways, 
    ranks, 
    minSize = fgsea.minSize, 
    maxSize = fgsea.maxSize
  )
  
  res <- fgseaRes[, ..result.cols]
  list(ranks=as.list(ranks), pathways=res)
}


# Used to ping this service to see if it responds
#* @get /
#* @serializer unboxedJSON
function(req) {
  list(ok=TRUE)
}


# curl --data-binary @brca_hd_tep_ranks_100.rnk -X POST "http://127.0.0.1:9404/preranked" -H "Content-Type: text/tab-separated-values"
#* @post /v1/preranked
#* @parser tsv
#* @serializer unboxedJSON
function(req) {
  ranks <- req$body
  colnames(ranks) <- c("gene","rank")
  ranks <- setNames(ranks$rank, ranks$gene)
  
  runFgseaPreranked(ranks)
}


# curl --data-binary @FakeExpression.txt -X POST "http://127.0.0.1:3723/rnaseq?classes=A,A,A,B,B,B" -H "Content-Type: text/tab-separated-values"
#* @post /v1/rnaseq
#* @parser tsv
#* @serializer unboxedJSON
function(req, classes) {
  RNASeq <- req$body
  
  # Parse the 'classes' query parameter into a vector
  classes <- strsplit(classes, split=",",fixed=TRUE)[[1]]
  
  # Boolean vector indicating which columns to keep
  columns.keep <- classes %in% c('A','B')
  
  # Drop description column
  RNASeq <- RNASeq[ , !(tolower(names(RNASeq)) %in% c('description'))]
  
  # Drop columns that are to be ignored
  classes <- classes[columns.keep]
  RNASeq <- RNASeq[, c(TRUE, columns.keep)]
  #print(classes)
  #str(RNASeq)
  
  runFgseaRnaseq(RNASeq, classes)
}
