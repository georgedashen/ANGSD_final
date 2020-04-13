library(DESeq2)
library(magrittr)

rownames(gene_count) <- gene_count$Geneid
gene_count <- gene_count[,c(2:16)]
gene_count <- gene_count[rowSums(gene_count)!=0,]

sample_info <- DataFrame(Cells = Info$Cells, Time = Info$Time, row.names = colnames(gene_count))

DESeq.ds <- DESeqDataSetFromMatrix(countData = gene_count,
                                   colData = sample_info,
                                   design = ~ Cells + Time)

DESeq.rlog <- rlog(DESeq.ds, blind = T)

corr_coeff <- cor(assay(DESeq.rlog), method = "pearson")

as.dist(1-corr_coeff) %>% hclust %>% plot(.,labels = colnames(assay(DESeq.rlog)), main = "rlong transformed read counts")

as.dist(1-corr_coeff, upper = T) %>% as.matrix %>% pheatmap::pheatmap(.,main = "Pearson correlation")

par(mfrow=c(1,2))
plotPCA(DESeq.rlog, intgroup = "Cells") 
plotPCA(DESeq.rlog, intgroup = "Time")

library(biomaRt)
require(org.Hs.eg.db)

symbols <- select(org.Hs.eg.db, rownames(gene_count), c("SYMBOL","GENENAME","ENTREZID"),"ENSEMBL")

#library(pcaExplorer)
#pcaExplorer(dds = DESeq.ds, dst = DESeq.rlog)