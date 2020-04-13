Samples <- paste0(rep("SRR87699",15),c(35:49))
Cells <- c(rep("DBTRG",5),rep("U87",5),rep("LN18",5))
Time <- as.character(rep(c(0,6,12,24,48),3))
Info <- data.frame(Samples=Samples, Cells=Cells, Time=Time)
print(Info)

library(ggplot2)
library(tidyr)
library(magrittr)
library(dplyr)

summary <- read.csv("./counts.txt.summary", header = T, sep='\t')
samples <- paste0(Cells, "_", rep(c(1:5),3))
colnames(summary) <- c("Status", samples)
summary <- gather(data = summary,
                  key = "samples",
                  value = "Nreads",
                  "DBTRG_1":"LN18_5")
summary <- summary[which(summary$Nreads!=0),]

gene_count <- read.csv("./counts.txt", header = T, sep = '\t')
dim(gene_count)

gene_count <- gene_count[,-c(2:6)]
colnames(gene_count) <- c(colnames(gene_count)[1], samples)

gene_count_df <- gather(data = gene_count,
                        key = "samples",
                        value = "Nreads",
                        "DBTRG_1" : "LN18_5")
gene_count_df$cells <- gsub("_[0-9]", "", gene_count_df$samples)
gene_count_df <- gene_count_df[which(gene_count_df$Nreads!=0),]

ggplot(summary, aes(fill=Status,y=Nreads,x=samples)) + geom_bar(position="dodge", stat="identity") + coord_flip() + theme(legend.position = "bottom") + labs(title = "featureCounts Summary Statistics")

mean_read <- gene_count_df %>% group_by(cells) %>% summarise(Nreads = mean(Nreads))

ggplot(mean_read, aes(x=cells, y=Nreads)) + geom_bar(stat = "identity") + labs(title = "Mean Number of Read for Each Cell Lines")

#ggplot(gene_count_df, aes(x=cells, y=Nreads)) + geom_point() + labs(title = "Read Distribution for each")

ggplot(gene_count_df, aes(x=log(Nreads+1), fill=cells))+ geom_histogram(alpha = 0.6, position = "identity") + labs(title = "Read Distribution for each")