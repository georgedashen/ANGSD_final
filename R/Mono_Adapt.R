S_A=apply((abs(DEG_S-DEG_S[,1])),1,max)/(abs(DEG_S[,5]-DEG_S[,1])+0.1)  # Adaptation
S_M=abs(DEG_S[,5]-DEG_S[,1])/apply(DEG_S,1,max)  # Monotonicity

R_A=apply((DEG_R-DEG_R[,1]),1,max)/(abs(DEG_R[,5]-DEG_R[,1])+0.1)  # Adaptation
R_M=abs(DEG_R[,5]-DEG_R[,1])/apply(DEG_R,1,max)  # Monotonicity

S_A_M = cbind(S_A[-which.max(S_A)], S_M[-which.max(S_A)])
R_A_M = cbind(R_A[-which.max(R_A)], R_M[-which.max(R_A)])

dev.new()
plot((S_A_M), pch = 23, col = "aquamarine3", cex = 1, xlab= "Adaptive response", ylab= "Monotoic response")
par(new=T) 
plot((R_A_M), pch = 21, col = "darkgoldenrod",  cex = 1.5, xlab= "", ylab= "")

par(mfrow=c(1,2))
plot(S_A[-which.max(S_A)], ylim=c(0,max(S_A[-which.max(S_A)], R_A[-which.max(R_A)])+100))
plot(R_A[-which.max(R_A)], ylim=c(0,max(S_A[-which.max(S_A)], R_A[-which.max(R_A)])+100))

par(mfrow=c(1,2))
hist(S_A[-which.max(S_A)])
hist(R_A[-which.max(R_A)])

par(mfrow=c(1,2))
h <- hist(log(R_A[-which.max(R_A)], base=20))
hist(log(S_A[-which.max(S_A)], base=20), breaks = h$breaks)

S_A.scaled = log(S_A[-which.max(S_A)], base=20)
R_A.scaled = log(R_A[-which.max(R_A)], base=20)

dataset_A <- data.frame(value = c(S_A.scaled, R_A.scaled), group = factor(rep(c("Sensitive TCGs","Resistant TCGs"), times = c(length(S_A.scaled), length(R_A.scaled)))))
boxplot( t(value) ~ t(group),  notch = F, dataset_A, outline = FALSE, border = c("darkgoldenrod","aquamarine3"),cex = 1, ylab= "Adaptive response", cex.axis=1,pars = list(boxwex = 0.5, staplewex = 0.5, outwex = 0.5))  #,col.axis = "#009E73"

wilcox.test(S_A,R_A,alternative="less")  # Wilcoxon signed rank test with continuity correction;

dataset_M <- data.frame(value = c(S_M[-which.max(S_A)],R_M[-which.max(R_A)]), group = factor(rep(c("Sensitive TCGs","Resistant TCGs"), times = c(length(R_M)-1, length(S_M)-1))))
boxplot( t(value) ~ t(group),  notch = F, dataset_M, outline = FALSE, border = c( "darkgoldenrod","aquamarine3"), cex = 1,  ylab= "Monotonic response", cex.axis=1,pars = list(boxwex = 0.5, staplewex = 0.5, outwex = 0.5))  #,col.axis = "#009E73"

wilcox.test(S_M,R_M,alternative="greater")