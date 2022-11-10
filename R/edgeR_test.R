# https://bioconductor.org/packages/release/bioc/html/edgeR.html
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6607905/#SM1title

#install.packages("BiocManager")
#BiocManager::install(c('Biobase','limma','edgeR','GSA','locfit'))
#install.packages(c('pheatmap', 'RColorBrewer', 'gProfileR', 'httr', 'RJSONIO'))

library("fgsea")
library("edgeR")
edgeRUsersGuide()

working_dir <- file.path(getwd(), "data")

# Raw RNA-seq counts
RNASeq <- read.table(
  file.path(working_dir, "NIHMS1032148-supplement-Supplementary_Table_10.txt"),
  header = TRUE, sep = "\t", quote="\"", stringsAsFactors = FALSE)

# Class definitions, groups the columns in the RNA-seq counts
classDefinitions_RNASeq <- read.table(
  file.path(working_dir, "NIHMS1032148-supplement-Supplementary_Table_11.txt"),
  header = TRUE, sep = "\t", quote="\"", stringsAsFactors = FALSE)



# Filter out genes with low read counts, as these are just noise

# counts per million
cpms <- cpm(RNASeq)
# determine which rows to keep 
# a gene mush have at least 50 measurements with more than 1 CPM in one of the classes to be included in the analysis
keep <- rowSums(cpms > 1) >= 50
# These are the read counts that pass the filtering
counts <- RNASeq[keep,]



# Normalize Data

# create data structure to hold counts and subtype information for each sample.
d <- DGEList(counts=counts, group=classDefinitions_RNASeq$SUBTYPE)
#Normalize the data
d <- calcNormFactors(d)
#calculate dispersion
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)



# Filter unannotated genes
exclude <- grep("\\?|^LOC", rownames(d), value=T)
d <- d[which(!rownames(d) %in% exclude),]



#calculate differential expression statistics with a simple design
de <- exactTest(d, pair=c("Immunoreactive", "Mesenchymal"))
tt_exact_test <- topTags(de,n=nrow(d))
tt <- tt_exact_test


# Create GSEA input list (.rnk file)
#calculate ranks
ranks_RNAseq = sign(tt$table$logFC) * -log10(tt$table$PValue)

#Separate gene names from IDs
genenames <- unlist(lapply( rownames(tt$table), function(data)
  {unlist(strsplit(data,"\\|"))[1]}))
#geneids <- unlist(lapply( rownames(tt$table),
#  function(data) {unlist(strsplit(data,"\\|"))[2]}))

#create ranks
ranks_RNAseq <- data.frame(genenames, ranks_RNAseq)
colnames(ranks_RNAseq) <- c("gene","rank")
ranks_RNAseq <- setNames(ranks_RNAseq$rank, ranks_RNAseq$gene)


# Run FGSEA
# load GMT file
gmt.file = file.path(working_dir, "Human_GOBP_AllPathways_no_GO_iea_June_01_2022_symbol.gmt")
pathways <- gmtPathways(gmt.file)

# Run FGSEA
fgseaRes <- fgsea(
  pathways, 
  ranks_RNAseq, 
  minSize = 15, 
  maxSize = 500
)

enr.file = file.path(working_dir, "Supplemental_Protocol_Immunoreactive_vs_Mesenchymal.txt")
write.table(fgseaRes[,c("pathway", "pathway", "pval", "padj")], 
            file=enr.file,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
