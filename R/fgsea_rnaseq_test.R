library(fgsea)
library(edgeR)
library(here)

# Ignore gene sets that are smaller or larger than the limits.
fgsea.minSize = 15
fgsea.maxSize = 500

as_matrix <- function(x) {
  y <- as.matrix.data.frame(x[,-1])
  rownames(y) <- x[[1]]
  y
}

# Load pathway database into memory.
gmt.file = here("data", "Human_GOBP_AllPathways_no_GO_iea_June_01_2022_symbol.gmt")
pathways <- gmtPathways(gmt.file)

exp.file = here("data", "FakeExpression_bad.txt")

# Fails here if the table can't be read...
RNASeq <- read.table(exp.file, header=TRUE, colClasses = c("character", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
RNASeq <- RNASeq[ , !(tolower(names(RNASeq)) %in% c("description"))]
RNASeq <- as_matrix(RNASeq)

classes <- c('Exp1A', 'Exp1B', 'Exp1C', 'Exp2A', 'Exp2B', 'Exp2C')

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

res <- fgseaRes[, c("pathway", "size", "pval", "padj", "ES", "NES")]
list(ranks=as.list(ranks), pathways=res)
