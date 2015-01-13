library(ggplot2)
library(gridExtra)
library(grid)

setwd("/Users/nigel/git/cafe-quality/NotTracked/master_full")

getld <-function(fname) {
  d= read.csv(fname)
  #remove duds
  d = d[!is.na(d$Reference) & d$Reference != "SmrtBellSequence",]
  d$ErrorRate = d$NumErrors/d$Length
  d$SNPErrorRate = d$NumSNPErrors / d$Length
  d$IndelErrorRate = d$NumIndelErrors / d$Length
  d$MeanGC = (d$SnrG + d$SnrC) / 2
  #d=d[d$Reference=="lambda_NEB3011",]
  #d=d[d$Reference=="HP.V1.02",]
  return(d)
}
d = getld("master_5ba5286_combined_reads.csv")

head(d)

q = grep("Snr",colnames(d))
snr = d[,q]
head(snr)

b = prcomp(snr, center=TRUE,scale=FALSE)
summary(b)
plot(b)
head(b)
str(b)
head(b$x)
hist(b$x[,1])
plot(b$x)

t=as.matrix(snr)%*%b$rotation
# Let's look at loadings
t=as.matrix(snr)%*%b$rotation