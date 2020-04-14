TCGs <- as.character(union(rownames(DEG_S), rownames(DEG_R)))
norm.TCG_S <- norm.counts[TCGs, 1:5]
norm.TCG_R <- norm.counts[TCGs, 11:15]

for(i in 1:nrow(norm.TCG_S)) {
  norm.TCG_S[i,] <- norm.TCG_S[i,] / max(norm.TCG_S[i,]+0.1)
}
colnames(norm.TCG_S) <- c("0 h", "6 h","12 h", "24 h", "48 h")

for(i in 1:nrow(norm.TCG_R)) {
  norm.TCG_R[i,] <- norm.TCG_R[i,] / max(norm.TCG_R[i,]+0.1)
}
colnames(norm.TCG_R) <- c("0 h", "6 h","12 h", "24 h", "48 h")

write.csv(norm.TCG_S, file = 'Net_TCG_S.csv')
write.csv(norm.TCG_R, file = 'Net_TCG_R.csv')

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol"), values = TCGs, mart = mart)
#which(duplicated(G_list$ensembl_gene_id)==T)
#1486 CCL3L1
#Through STRING query, CCL3L3 is more commonly used
G_list <- G_list[-(which(duplicated(G_list$ensembl_gene_id)==T)-1),]
G_list_filtered <- G_list[-which(G_list$hgnc_symbol==""),]
#store G_list_filtered and use to attain PPI

TCGs <- read.table("Net_genes.csv", header = F)
PPI <- read.csv("string_interactions.tsv", sep = '\t')
PPI <- PPI[,1:2]
Net_genes <- union(PPI[,1], PPI[,2])
remove <- which(Net_genes %in% TCGs$V2==FALSE)
Net_size <- length(Net_genes)

Interact_mat <- matrix(0, Net_size, Net_size)
rownames(Interact_mat) <- Net_genes
colnames(Interact_mat) <- Net_genes
for(i in 1:nrow(PPI)) {
  Interact_mat[PPI[i,1],PPI[i,2]] = 1
}
Interact_mat <- Interact_mat[-remove, -remove]

Net_genes <- Net_genes[-remove]
Net_size <- length(Net_genes)

#Data interpolated

TCG1 <- t(read.csv('Net_TCG_S_interpolated.csv', header = F))
TCG2 <- t(read.csv('Net_TCG_R_interpolated.csv', header = F))
rownames(G_list) <- G_list$ensembl_gene_id
colnames(TCG1) <- G_list[TCGs,]$hgnc_symbol
colnames(TCG2) <- colnames(TCG1)

TCG1 <- TCG1[,Net_genes]
TCG2 <- TCG2[,Net_genes]
TCG1 <- TCG1[,-which(colSums(TCG1)==0)]
TCG2 <- TCG2[,-which(colSums(TCG2)==0)]

cor_S <- matrix(0, ncol(TCG1), ncol(TCG1))
ppval_S <- cor_S
cor_R <- matrix(0, ncol(TCG2), ncol(TCG2))
ppval_R <- cor_S

for(i in 1:ncol(TCG1)) {
  for(j in 1:ncol(TCG1)){
    partial_test <- cor.test(TCG1[,i], TCG1[,j], method = 'pearson')
    cor_S[i,j] <- partial_test$estimate
    ppval_S[i,j] <- partial_test$p.value
  }
}

for(i in 1:ncol(TCG2)) {
  for(j in 1:ncol(TCG2)){
    partial_test <- cor.test(TCG2[,i], TCG2[,j], method = 'pearson')
    cor_R[i,j] <- partial_test$estimate
    ppval_R[i,j] <- partial_test$p.value
  }
}

edge_S <- matrix(0, ncol(TCG1), ncol(TCG1))
colnames(edge_S) <- colnames(TCG1)
rownames(edge_S) <- colnames(TCG1)
for (i in 1:ncol(TCG1)){
  for (j in 1:ncol(TCG1)){
    if (cor_S[i,j]>0.75 & ppval_S[i,j]*choose(ncol(TCG1),2)<0.05 & Interact_mat[colnames(TCG1)[i],colnames(TCG1)[j]]==1)
      edge_S[i,j]=1
    else if (cor_S[i,j]<-0.75 & ppval_S[i,j]*choose(ncol(TCG1),2)<0.05 & Interact_mat[colnames(TCG1)[i],colnames(TCG1)[j]]==1)
      edge_S[i,j]=-1
  }
}

edge_R <- matrix(0, ncol(TCG2), ncol(TCG2))
colnames(edge_R) <- colnames(TCG2)
rownames(edge_R) <- colnames(TCG2)
for (i in 1:ncol(TCG2)){
  for (j in 1:ncol(TCG2)){
    if (cor_R[i,j]>0.75 & ppval_R[i,j]*choose(ncol(TCG2),2)<0.05 & Interact_mat[colnames(TCG2)[i],colnames(TCG2)[j]]==1)
      edge_R[i,j]=1
    else if (cor_R[i,j]<-0.75 & ppval_R[i,j]*choose(ncol(TCG2),2)<0.01 & Interact_mat[colnames(TCG2)[i],colnames(TCG2)[j]]==1)
      edge_R[i,j]=-1
  }
}

sum(edge_S!=0)
sum(edge_R!=0)

edge_S <- edge_S[which(rowSds(edge_S)!=0), which(rowSds(edge_S)!=0)] #710
edge_R <- edge_R[which(rowSds(edge_R)!=0), which(rowSds(edge_R)!=0)] #611

library(igraph)
#ignore direction
edge_S[which(edge_S == -1)] = 1
edge_R[which(edge_R == -1)] = 1

net_S <- graph_from_adjacency_matrix(edge_S, diag = F, mode="undirected")
deg_S <- degree(net_S, mode="all")

net_R <- graph_from_adjacency_matrix(edge_R, diag = F, mode="undirected")
deg_R <- degree(net_R, mode="all")

par(mfrow=c(1,2))
h <- hist(deg_S, main="Histogram of s_net node degree")
hist(deg_R, breaks = h$breaks, main="Histogram of r_net node degree")

mean(deg_S)
mean(deg_R)