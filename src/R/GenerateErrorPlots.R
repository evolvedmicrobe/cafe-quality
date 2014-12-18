library(ggplot2)
library(gridExtra)
library(grid)

d= read.csv("/Users/nigel/git/cafe-quality/master_full/local_read_combined_reads.csv")
#remove duds
d = d[!is.na(d$Reference) & d$Reference != "SmrtBellSequence",]
head(d)
d[69830,]
k = levels(d$Reference)
sum(d$Reference=="lambda_NEB3011")
head(lambda)

d$ErrorRate = d$NumErrors/d$Length
d$SNPErrorRate = d$NumSNPErrors / d$Length
d$IndelErrorRate = d$NumIndelErrors / d$Length


ggplot

mkPlot<-function(errorType, ylab, title) {
  fm = formula(paste(errorType,"~NumPasses+Reference"))
  q = aggregate(fm,d,mean)
  v = ggplot(q,aes_string(x="NumPasses", y = errorType, colour="Reference")) + geom_line()+ theme_bw(base_size=8)
  v = v + labs(x= "Number of Passes", y = paste(ylab, "Error Rate (per bp)"), title=title)
  v = v + scale_x_continuous(limits=c(0,100))
  v
}


a = mkPlot("ErrorRate","Total ", "Total Errors")
b = mkPlot("SNPErrorRate","SNP ", "SNP Errors")
c = mkPlot("IndelErrorRate","SNP ", "Indel Errors")



q = aggregate(ErrorRate~NumPasses+Reference, d[d$Reference=="lambda_NEB3011",], mean )
t1= rep("All",nrow(q))

q2 = aggregate(IndelErrorRate~NumPasses+Reference, d[d$Reference=="lambda_NEB3011",], mean )
names(q2)[3]<-"ErrorRate"
t2 = rep("Indel",nrow(q2))

q3 = aggregate(SNPErrorRate~NumPasses+Reference, d[d$Reference=="lambda_NEB3011",], mean )
names(q3)[3]<-"ErrorRate"
t3 = rep("SNP",nrow(q3))

res = rbind(q,q2,q3)
res$ErrorType=c(t1,t2,t3)

q = aggregate(ErrorRate~NumPasses+Reference, d[d$Reference=="lambda_NEB3011",], length )
t1= rep("All",nrow(q))

q2 = aggregate(IndelErrorRate~NumPasses+Reference, d[d$Reference=="lambda_NEB3011",], length )
names(q2)[3]<-"ErrorRate"
t2 = rep("Indel",nrow(q2))

q3 = aggregate(SNPErrorRate~NumPasses+Reference, d[d$Reference=="lambda_NEB3011",], length )
names(q3)[3]<-"ErrorRate"
t3 = rep("SNP",nrow(q3))

res2 = rbind(q,q2,q3)
res2$ErrorType=c(t1,t2,t3)


head(res)
ph = function(x) { -10*log10(x)}
res$QV = sapply(res$ErrorRate,ph)

d$JittedPos = d$NumPasses+.5 - runif(nrow(d))
ld =d[d$Reference=="lambda_NEB3011",]
n = nrow(ld)
sd = ld[sample(n, n*.05),]

sum(ld$NumPasses>9)/nrow(ld)

v = ggplot(res,aes(x=NumPasses, y = QV, colour=ErrorType))+geom_line()+geom_point()+theme_classic()+scale_x_continuous(limits=c(0,60)) 
v = v+labs(x="Number of Passes", y = "Per Base Error Rate (Phred QV Scale)",title=expression(paste("Accuracy in ",lambda)))
v = v + geom_vline(xintercept=10, colour="red") +geom_rug(aes(x=JittedPos, y = NULL, colour=NULL),data=sd)
v



setwd("/Users/nigel/git/cafe-quality/master_full/")
pdf("errors.pdf",width=11,height=8.5)
grid.arrange(a,b,c,v, ncol = 2)
dev.off()

hist(lb$PredictedCCSAccuracy,40)
qplot(lb$Length,lb$PredictedCCSAccuracy)+geom_smooth()
qplot(lb$NumSubReads,lb$PredictedCCSAccuracy)+geom_smooth()

hist(lb$ErrorRate,50)
median(lb$ErrorRate)
