DESeq.ds <- DESeq(DESeq.ds)
norm.counts <- counts(DESeq.ds, normalized=T)

delta <- 5
epsilon <- 10

TCG <- rep(F, nrow(norm.counts))
for(i in 1:nrow(norm.counts)) {
  if(max(norm.counts[i,1:5])>=epsilon){
    isTCG <- FALSE
    fc <- sort(norm.counts[i,1:5])[5] / sort(norm.counts[i,1:5])[1]
    if(fc > delta) {
      isTCG <- TRUE
    }
    TCG[i] <- isTCG
  }
}
DEG_S <- norm.counts[TCG,1:5]

TCG <- rep(F, nrow(norm.counts))
for(i in 1:nrow(norm.counts)) {
  if(max(norm.counts[i,11:15])>=epsilon){
    isTCG <- FALSE
    fc <- sort(norm.counts[i,11:15])[5] / sort(norm.counts[i,11:15])[1]
    if(fc > delta) {
      isTCG <- TRUE
    }
    TCG[i] <- isTCG
  }
}
DEG_R <- norm.counts[TCG,11:15]

norm.DEG_S <- DEG_S
for(i in 1:nrow(DEG_S)) {
  norm.DEG_S[i,] <- DEG_S[i,] / max(DEG_S[i,])
}
colnames(norm.DEG_S) <- c("0 h", "6 h","12 h", "24 h", "48 h")

norm.DEG_R <- DEG_R
for(i in 1:nrow(DEG_R)) {
  norm.DEG_R[i,] <- DEG_R[i,] / max(DEG_R[i,])
}
colnames(norm.DEG_R) <- c("0 h", "6 h","12 h", "24 h", "48 h")

par(mfrow=c(1,2))
pheatmap::pheatmap(norm.DEG_S, main = "Expression Profile of DBTRG TCGs", show_rownames = F, cluster_cols = F)
pheatmap::pheatmap(norm.DEG_R, main = "Expression Profile of LN18 TCGs", show_rownames = F, cluster_cols = F)