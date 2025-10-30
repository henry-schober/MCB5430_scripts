# Part 2 Differential Gene Expression Analysis

#Set working Directory
setwd("//wsl.localhost/Ubuntu/home/henrys/MCB_5430/final/DE_tables")

# 2.1

library(DESeq2)
library(ggrepel)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)


samples <- c( "E2_rep1", "E2_rep2", "untr_rep1", "untr_rep2")
conditions <- c("Treated", "Treated", "Untreated", "Untreated")
counts <- as.matrix(read.csv("E2_gene_matrix.csv", row.names = "gene_id"))

head(counts)

# DESeq

dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = as.data.frame(conditions), design = ~ conditions)

dds <- DESeq(dds)

dds

#Results

res <- results(dds, alpha=0.05)

head(res)

res <- res[order(res$padj),]

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by = 'row.names', sort=FALSE)
names(resdata)[1] <- 'gene'

head(resdata,3)

#remove na's
res.narm <- na.omit(resdata)

dim(res.narm)

#select for significant genes
sig.res <- res.narm[res.narm$padj<0.05,]
dim(sig.res)

#differentiate between up and down regulated genes
sig.res.up <- sig.res[sig.res$log2FoldChange>0,]
sig.res.down <- sig.res[sig.res$log2FoldChange<0,]

dim(sig.res.down)
dim(sig.res.up)

# 2.2

#log transform the data

logdds <- rlogTransformation(dds, blind=F)
head(logdds)


PCA_plot <- plotPCA(logdds, intgroup=c('conditions')) + geom_label(aes(label = samples), nudge_y = 5) +
  labs(title = "PCA plot of Treated vs Untreated replicates")
PCA_plot

#2.3

plotMA(dds, alpha=0.05, ylim=c(-28,28))


plot(log2(res.narm$baseMean), res.narm$log2FoldChange, pch=20, cex=.5, main = "MA plot on changed genes", xlab = "Counts", ylab = "Log Fold Change")
points(log2(sig.res$baseMean), sig.res$log2FoldChange, pch=20, col='red')
points(log2(sig.res.up$baseMean), sig.res.up$log2FoldChange, pch=20, col='blue')
points(log2(sig.res.down$baseMean), sig.res.down$log2FoldChange, pch=20, col='green')
dev.copy(jpeg, filename="Final_MA.jpeg")
dev.off()


# 2.4

#reformat this for splitting
E2genes_all <- as.character(sapply(resdata$gene, function(x) unlist(strsplit(x, split="[|]"))[2]))

E2genes_up <- as.character(sapply(sig.res.up$gene, function(x) unlist(strsplit(x, split="[|]"))[2]))

E2genes_down <- as.character(sapply(sig.res.down$gene, function(x) unlist(strsplit(x, split="[|]"))[2]))

E2genes_allc <- as.character(E2genes_all)
E2genes_upc <- as.character(E2genes_up)
E2genes_downc <- as.character(E2genes_down)

# 2.5 converting refseq to entrez

all_entrezID <- unique(mapIds(x= org.Hs.eg.db, column = "ENTREZID", keytype = "REFSEQ", multiVals = "first", keys = E2genes_allc))

up_entrezID <- unique(mapIds(x= org.Hs.eg.db, column = "ENTREZID", keytype = "REFSEQ", multiVals = "first", keys = E2genes_upc))
down_entrezID <- unique(mapIds(x= org.Hs.eg.db, column = "ENTREZID", keytype = "REFSEQ", multiVals = "first", keys = E2genes_downc))

# enrichment from cluster profiler with the 3 options; BP, MF, CC

#BP

up_go_enrich_BP <- enrichGO(gene = up_entrezID, universe = all_entrezID,
                         OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", readable = T,
                         ont = "BP", pAdjustMethod = "fdr", qvalueCutoff = 0.05)

down_go_enrich_BP <- enrichGO(gene = down_entrezID, universe = all_entrezID,
                         OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", readable = T,
                         ont = "BP", pAdjustMethod = "fdr", qvalueCutoff = 0.05)

#MF

up_go_enrich_MF <- enrichGO(gene = up_entrezID, universe = all_entrezID,
                            OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", readable = T,
                            ont = "MF", pAdjustMethod = "fdr", qvalueCutoff = 0.05)

down_go_enrich_MF <- enrichGO(gene = down_entrezID, universe = all_entrezID,
                              OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", readable = T,
                              ont = "MF", pAdjustMethod = "fdr", qvalueCutoff = 0.05)

#CC

up_go_enrich_CC <- enrichGO(gene = up_entrezID, universe = all_entrezID,
                            OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", readable = T,
                            ont = "CC", pAdjustMethod = "fdr", qvalueCutoff = 0.05)

down_go_enrich_CC <- enrichGO(gene = down_entrezID, universe = all_entrezID,
                              OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", readable = T,
                              ont = "CC", pAdjustMethod = "fdr", qvalueCutoff = 0.05)

head(as.data.frame(up_go_enrich_BP))
head(as.data.frame(up_go_enrich_MF))
head(as.data.frame(up_go_enrich_CC))

# Bar plots

barplot(up_go_enrich_BP, showCategory = 20, title = "BP UP regulated enrichment",
        font.size = 4, label_format = 30)

barplot(down_go_enrich_BP, showCategory = 20, title = "BP Down regulated enrichment",
        font.size = 4, label_format = 30)

barplot(up_go_enrich_MF, showCategory = 20, title = "MF UP regulated enrichment",
        font.size = 4, label_format = 30)

barplot(down_go_enrich_MF, showCategory = 20, title = "MF Down regulated enrichment",
        font.size = 4, label_format = 30)

barplot(up_go_enrich_CC, showCategory = 20, title = "CC UP regulated enrichment",
        font.size = 4, label_format = 30)

barplot(down_go_enrich_CC, showCategory = 20, title = "CC Down regulated enrichment",
        font.size = 4, label_format = 30)


# answer to 2.6

#Given what I know about estradiol and estrogen itself, I believe that the Up regulated genes under the BP category best match the known functions of these genes. This is because Estradiol is a primary growth hormone, and the largest Gene ratio that we see has yo do with gland development, as well as mammary gland development was also one of the top clusters. There was a lot of mammary gland development in this enrichment.



