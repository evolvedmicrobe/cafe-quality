library(ggplot2)
library(gridExtra)
library(grid)

setwd("/Users/nigel/git/cafe-quality/master_full/")

getld <-function(fname) {
  d= read.csv(fname)
  #remove duds
  d = d[!is.na(d$Reference) & d$Reference != "SmrtBellSequence",]
  d$ErrorRate = d$NumErrors/d$Length
  d$SNPErrorRate = d$NumSNPErrors / d$Length
  d$IndelErrorRate = d$NumIndelErrors / d$Length
  d$MeanGC = (d$SnrG + d$SnrC) / 2
  d[d$Reference=="lambda_NEB3011",]
}
d = getld("/Users/nigel/git/cafe-quality/master_full/local_read_combined_reads.csv")
d2 = getld("/Users/nigel/git/cafe-quality/master_full/C2_polymer_combined_reads.csv")

head(d2)

in1 = intersect(d$ZMW,d2$ZMW)
d2 = d2[d2$ZMW%in%in1,]
d = d[d$ZMW%in%in1,]

cd = rbind(d,d2)
cd$Reference=  factor(c(rep("Original",nrow(d)),rep("Polished",nrow(d2))))

mkPlot<-function(errorType, ylab, title) {
  fm = formula(paste(errorType,"~NumPasses+Reference"))
  q = aggregate(fm,cd,mean)
  v = ggplot(q,aes_string(x="NumPasses", y = errorType, colour="Reference")) + geom_line()+ theme_bw(base_size=8)
  v = v + labs(x= "Number of Passes", y = paste(ylab, "Error Rate (per bp)"), title=title)
  v = v + scale_x_continuous(limits=c(0,100))
  v
}


a = mkPlot("ErrorRate","Total ", "Total Errors")
b = mkPlot("SNPErrorRate","SNP ", "SNP Errors")
c = mkPlot("IndelErrorRate","SNP ", "Indel Errors")

grid.arrange(a,b,c)

q = aggregate(ErrorRate~NumPasses+Reference, cd, mean )
t1= rep("All",nrow(q))

q2 = aggregate(IndelErrorRate~NumPasses+Reference, cd, mean )
names(q2)[3]<-"ErrorRate"
t2 = rep("Indel",nrow(q2))

q3 = aggregate(SNPErrorRate~NumPasses+Reference, cd, mean )
names(q3)[3]<-"ErrorRate"
t3 = rep("SNP",nrow(q3))

res = rbind(q,q2,q3)
res$ErrorType=c(t1,t2,t3)



head(res)
ph = function(x) { -10*log10(x)}
res$QV = sapply(res$ErrorRate,ph)

res = res[res$ErrorType=="All",]

v = ggplot(res,aes(x=NumPasses, y = QV, colour=Reference))+geom_line()+geom_point()+theme_classic()+scale_x_continuous(limits=c(0,60)) 
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


## Variant analysis
dv = read.csv("/Users/nigel/git/cafe-quality/master_full/master_5ba5286_all_variants.csv")
ldv = dv[dv$Ref=="lambda_NEB3011",]
ldv2 = read.csv("/Users/nigel/git/cafe-quality/data/C2_polymer_fix_1934603_all_variants.csv" )
ldv2 = ldv2[ldv2$Ref=="lambda_NEB3011",]

cdv = rbind(ldv,ldv2)
cdv$Ref = c(rep("Original",nrow(ldv)),rep("Polished",nrow(ldv2)))
cdv = cdv[cdv$zmw%in%in1,]
cdv$homopolymerChar[cdv$homopolymerChar=="T"]="A"
cdv$homopolymerChar[cdv$homopolymerChar=="G"]="C"


errorDrop <-function(minPasses) {
  ok = cd$ZMW[cd$NumPasses>=minPasses]
  gd = cdv[cdv$zmw%in%ok,] 
  head(gd)
  res = aggregate(Pos~Ref,gd,length)
  print(res)
  res[2,2]/res[1,2]
}

s=seq(0,30)
y = sapply(s,errorDrop)
b=data.frame(NumPasses = s, pr = y)
pdf(height=4,width=4,"ImprovementsByCoverage.pdf")
ggplot(b,aes(x=NumPasses,y=pr))+geom_point()+labs(x="Minimum Number of Passes Required", y = "Percentage of Errors Removed", title ="How Polishing Effect Varies by Coverage Cutoff")+theme_bw(base_size=8)
dev.off()


errorDrop2 <-function(minGC) {
  ok = cd$ZMW[cd$MeanGC>=minGC & cd$Mean]
  gd = cdv[cdv$zmw%in%ok,] 
  head(gd)
  res = aggregate(Pos~Ref,gd,length)
  print(res)
  res[2,2]/res[1,2]
}

s=seq(0,15)
y2 = sapply(s,errorDrop2)
b=data.frame(NumPasses = s, pr = y2)

pdf(height=4,width=4,"ImprovementsByCoverage.pdf")
ggplot(b,aes(x=NumPasses,y=pr))+geom_point()+labs(x="Minimum Number of Passes Required",
                                                  y = "Percentage of Errors Removed",
                                                  title ="How Polishing Effect Varies by Coverage Cutoff")+theme_bw()
dev.off()


### CONTOUR PLOT
errorContour <-function(minGC, minPasses) {
  ok = cd$ZMW[cd$MeanGC>=(minGC - 2) & cd$MeanGC>=(minGC + 2) & cd$NumPasses>=minPasses]
  gd = cdv[cdv$zmw%in%ok,] 
  head(gd)
  res = aggregate(Pos~Ref,gd,length)
  print(res)
  res[2,2]/res[1,2]
}
k = Vectorize(errorContour)
x = seq(0,14,.5) # Min GC
y = seq(0,20) # Min Passes
z = outer(x,y,k)

pdf("ErrorsRemovedContour.pdf",width=6.5,height=5)
filled.contour(x,y,z, xlab="Mean G+C SNR",ylab="Minimum Number of Passes", main="% Errors Removed By SNR and Coverage")
dev.off()

pdf("LogScaleErrors.pdf", width=10.5, height=7.5)
par(mfrow=c(2,2))
errorCount <-function(minGC, minPasses) {
  ok = d[d$MeanGC>=(minGC - 2) & d$MeanGC>=(minGC + 2) & d$NumPasses>=minPasses,]
  sum(ok$NumIndelErrors)
}
k = Vectorize(errorCount)
x = seq(0,14,.5) # Min GC
y = seq(0,20) # Min Passes
z = outer(x,y,k)
filled.contour(x,y,z, xlab="Mean G+C SNR",
               ylab="Minimum Number of Passes", 
               main="Total Indel Errors by\nSNR and Coverage", 
               key.title=title(main="Count"))

filled.contour(x,y,log10(z), xlab="Mean G+C SNR",
               ylab="Minimum Number of Passes", 
               main="Log Total Indel Errors by\nSNR and Coverage", 
               key.title=title(main="Log10\nCount"))


readCount <-function(minGC, minPasses) {
  sum(d$Length[d$MeanGC>=(minGC - 2) & d$MeanGC>=(minGC + 2) & d$NumPasses>=minPasses & d$Length<2500])
}
k = Vectorize(readCount)
x = seq(0,14,.5) # Min GC
y = seq(0,20) # Min Passes
z = outer(x,y,k)

filled.contour(x,y,z, xlab="Mean G+C SNR",
               ylab="Minimum Number of Passes", 
               main="Total CCS Basepairs", 
               key.title=title(main="BP\nCount"))

errorRate <-function(minGC, minPasses) {
  ok = d[d$MeanGC>=(minGC - 2) & d$MeanGC>=(minGC + 2) & d$NumPasses>=minPasses & d$Length<2500,]
  mean(ok$IndelErrorRate)
}
k = Vectorize(errorRate)
x = seq(0,14,.5) # Min GC
y = seq(0,20) # Min Passes
z = outer(x,y,k)

filled.contour(x,y,log10(z), xlab="Mean G+C SNR",
               ylab="Minimum Number of Passes", 
               main="CCS Log10 Indel Error Rate")

dev.off()








ok = cd$ZMW[cd$MeanGC>=(minGC - 2) & cd$MeanGC>=(minGC + 2) & cd$NumPasses>=minPasses]
gd = cdv[cdv$zmw%in%ok,] 
head(gd)
res = aggregate(Pos~Ref+type+homopolymerLength+homopolymerChar,cdv,length)

head(res)

res2 = res[res$homopolymerChar=="C",]
head(res2)
pdf("LambdaIndelErrorsPolishedAllData.pdf",width=6,height=4)
v = ggplot(res2,aes(x=homopolymerLength,y=Pos, colour=Ref))+facet_grid( .~indeltype)+geom_point(size=1.5)+theme_bw(base_size=9)
v = v +labs(x="Homopolymer Length",y="Count of Total Errors", 
            title=expression(paste(lambda, " G/C indel errors divided by genomic context for all Reads")))
v = v + scale_colour_discrete(name="Method",labels = c("Original","Polished"))#+geom_vline(xintercept=3.5,colour="red")
v
dev.off()

or = res2[res2$Ref=="Original",]
e1d = sum(or$Pos[or$homopolymerLength>4 & or$indeltype=="Deletion"])
e1i = sum(or$Pos[or$homopolymerLength>3 & or$indeltype=="Insertion"])
e1=e1d+e1i

or2 = res2[res2$Ref=="Polished",]
e1d = sum(or2$Pos[or2$homopolymerLength>4 & or2$indeltype=="Deletion"])
e1i = sum(or2$Pos[or2$homopolymerLength>3 & or2$indeltype=="Insertion"])
e2=e1d+e1i

e2/e1

head(or)


