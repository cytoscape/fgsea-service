##
count= read.delim("/Users/mkucera/git/fgsea-service/data/GSE129943_rsem_counts_2016.txt")
count_mx =  as.matrix(count[ , -c(1,2)])
rownames(count_mx) = count$Ensembl
myGroups = c( "A", "A", "A", "B" ,"B", "B" )

# Convert to matrix
classes = myGroups
RNASeq <- count_mx

# Create data structure to hold counts and subtype information for each sample.
d <- DGEList(counts=RNASeq, group=classes)

# Filter out genes with low read counts, as these are just noise
keep <- filterByExpr(d)
d <- d[keep, , keep.lib.sizes=FALSE]

# Normalize the data
d <- calcNormFactors(d)
design <- model.matrix(~0 + myGroups )
d <- estimateDisp(d,design)

# Calculate differential expression statistics with a simple design
my.contrasts <- makeContrasts(AvsB=myGroupsA-myGroupsB, levels = design )
fit <- glmQLFit(d,design)
qlf <- glmQLFTest(fit,coef=2, contrast = my.contrasts)

table_res = topTags(qlf, n = nrow(d))

ranks =  sign(table_res$table$logFC) * -log10(table_res$table$PValue)

# Calculate ranks
#ranks = sign(tt$table$logFC) * -log10(tt$table$PValue)

# Put ranks into format for FGSEA
genenames <- rownames(table_res)
ranks <- data.frame(genenames, ranks)
colnames(ranks) <- c("gene","rank")
ranks <- setNames(ranks$rank, ranks$gene)

ranks