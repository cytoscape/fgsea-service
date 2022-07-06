library(plumber)
library(fgsea)
library(here)

set.seed(42)

#Load pathway database into memory.
gmt.file = here("data", "Human_GOBP_AllPathways_no_GO_iea_February_01_2017_symbol.gmt")
pathways <- gmtPathways(gmt.file)

# curl --data-binary @brca_hd_tep_ranks_100.rnk -X POST "http://127.0.0.1:9404/fgsea" -H "Content-Type: text/tab-separated-values"


#* @post /fgsea
#* @parser tsv
#* @serializer json
function(req) {
  ranks <- req$body
  colnames(ranks) <- c("gene","rank")
  ranks <- setNames(ranks$rank, ranks$gene)
  
  fgseaRes <- fgsea(
    pathways, 
    ranks, 
    minSize = 15, 
    maxSize = 500
  )
  
  # Add full gene sets to results
  fgseaRes$genes <- pathways[match(fgseaRes$pathway, names(pathways))]
  
  fgseaRes[1:100, c("pathway", "size", "pval", "ES", "NES", "leadingEdge", "genes")]
}
