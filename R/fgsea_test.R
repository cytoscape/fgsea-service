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

# Initialization

library(fgsea)
library(here)

set.seed(42)

# Load Pathway Commons Tutorial data

gmt.file = here("data", "Human_GOBP_AllPathways_no_GO_iea_February_01_2017_symbol.gmt")
rnk.file = here("data", "brca_hd_tep_ranks.rnk")

ranks <- read.table(rnk.file, header=TRUE, colClasses = c("character", "numeric"))
colnames(ranks) <- c("g","r")
ranks <- setNames(ranks$r, ranks$g)

pathways <- gmtPathways(gmt.file)


# Run FGSEA
#system.time({ fgsea(pathways, ranks, minSize = 15, maxSize = 500) })
fgseaRes <- fgsea(
  pathways, 
  ranks, 
  minSize = 15, 
  maxSize = 500
)

cat("number of enriched pathways:", nrow(fgseaRes))

# Output files for EnrichmentMap
# Output enrichments in "generic" format with just pvalues
write.table(fgseaRes[,c("pathway", "pathway", "pval")], 
            file=paste(here(), "/data/fgsea_PC_enrichments_generic_results.txt", sep=""),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)


