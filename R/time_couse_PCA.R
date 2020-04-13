library(MASS)
require(gridExtra)
DBTRG <- t(gene_count[,1:5])
pca <- prcomp(DBTRG %*% t(DBTRG))
pca.df <- pca$x[,1:3]
pca.df <- as.data.frame(scale(pca.df))
pca.df$labels <- factor(Time[1:5], level=c("0","6","12","24","48"))

g1 <- ggplot(pca.df, aes(x=PC1, y=PC2)) + geom_point(size = 5, aes(color=labels)) + geom_path() + labs(title = "DBTRG PC1 v.s. PC2")
g2 <- ggplot(pca.df, aes(x=PC1, y=PC3)) + geom_point(size = 5, aes(color=labels)) + geom_path() + labs(title = "DBTRG PC1 v.s. PC3")
grid.arrange(g1, g2, ncol = 2)

LN18 <- t(gene_count[,11:15])
pca <- prcomp(LN18 %*% t(LN18))
pca.df <- pca$x[,1:3]
pca.df <- as.data.frame(scale(pca.df))
pca.df$labels <- factor(Time[1:5], level=c("0","6","12","24","48"))

g1 <- ggplot(pca.df, aes(x=PC1, y=PC2)) + geom_point(size = 5, aes(color=labels)) + geom_path() + labs(title = "LN18 PC1 v.s. PC2")
g2 <- ggplot(pca.df, aes(x=PC1, y=PC3)) + geom_point(size = 5, aes(color=labels)) + geom_path() + labs(title = "LN18 PC1 v.s. PC3")
grid.arrange(g1, g2, ncol = 2)

U87 <- t(gene_count[,6:10])
pca <- prcomp(U87 %*% t(U87))
pca.df <- pca$x[,1:3]
pca.df <- as.data.frame(scale(pca.df))
pca.df$labels <- factor(Time[1:5], level=c("0","6","12","24","48"))

g1 <- ggplot(pca.df, aes(x=PC1, y=PC2)) + geom_point(size = 5, aes(color=labels)) + geom_path() + labs(title = "U87 PC1 v.s. PC2")
g2 <- ggplot(pca.df, aes(x=PC1, y=PC3)) + geom_point(size = 5, aes(color=labels)) + geom_path() + labs(title = "U87 PC1 v.s. PC3")
grid.arrange(g1, g2, ncol = 2)