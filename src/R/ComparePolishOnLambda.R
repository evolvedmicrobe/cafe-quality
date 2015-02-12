library(ggplot2)
library(gridExtra)
library(grid)

setwd("/Users/nigel/git/cafe-quality/NotTracked/")

getld <-function(fname) {
  d= read.csv(fname)
  #remove duds
  d = d[!is.na(d$Reference) & d$Reference != "SmrtBellSequence",]
  d$ErrorRate = d$NumErrors/d$Length
  d$SNPErrorRate = d$NumSNPErrors / d$Length
  d$IndelErrorRate = d$NumIndelErrors / d$Length
  d$MeanGC = (d$SnrG + d$SnrC) / 2
  d$MeanSnr = .25*(d$SnrG + d$SnrC + d$SnrT + d$SnrA) 
  #d[d$Reference=="lambda_NEB3011",]
  d[d$Reference=="ALL4MER.V2.01",]
}
d = getld("master_5ba5286_combined_reads.csv")
#d2 = getld("C2_polymer_d067558_combined_reads.csv")
#d2 = getld("del_tag_combined_reads.csv")
#d3 = getld("no_del_combined_reads.csv")
d2 = getld("16s_training_combined_reads.csv")


in1 = intersect(d$ZMW,d2$ZMW)
d2 = d2[d2$ZMW%in%in1,]
d = d[d$ZMW%in%in1,]

cd = rbind(d,d2)
cd$Analysis=  factor(c(rep("Original",nrow(d)),rep("New",nrow(d2))))
cd$MeanGC = .5*(cd$SnrG + cd$SnrC)

res = aggregate(ErrorRate~NumPasses+Analysis, cd, FUN=mean)
head(res)
ph = function(x) { -10*log10(x)}
res$QV = sapply(res$ErrorRate,ph)

v = ggplot(res,aes(x=NumPasses, y = QV, colour=Analysis))+geom_line()+geom_point()+theme_classic(base_size=8)+scale_x_continuous(limits=c(0,60)) 
v = v+labs(x="Number of Passes", y = "Per Base Error Rate (Phred QV Scale)",title=expression(paste("Accuracy in ",lambda)))
v = v + geom_vline(xintercept=10, colour="red") +geom_hline(yintercept=30, colour="red")
v

flipped = dcast(res,NumPasses~Analysis,value.var="ErrorRate")
nd = t(apply(flipped[,2:3],1, function(x) x/x[2]))
head(nd)
flipped[,2:3]<-nd
fmd = melt(flipped, id.vars="NumPasses")
v2 = ggplot(fmd,aes(x=NumPasses, y = value, colour=variable))+geom_line()+geom_point()+theme_classic(base_size=8)+scale_x_continuous(limits=c(0,60)) 
v2 = v2+labs(x="Number of Passes", y = "Relative Error Rate",title=expression(paste("Relative Accuracy in ",lambda)))+scale_y_continuous(labels=percent)
v2
pdf("LambdaComparison.pdf",width=10, height  =3.5)
grid.arrange(v,v2,ncol=2)
dev.off()

pdf("errors.pdf",width=11,height=8.5)
grid.arrange(a,b,c,v, ncol = 2)
dev.off()

hist(lb$PredictedCCSAccuracy,40)
qplot(lb$Length,lb$PredictedCCSAccuracy)+geom_smooth()
qplot(lb$NumSubReads,lb$PredictedCCSAccuracy)+geom_smooth()

hist(lb$ErrorRate,50)
median(lb$ErrorRate)


## Variant analysis
dv = read.csv("master_5ba5286_all_variants.csv")
ldv = dv[dv$Ref=="lambda_NEB3011",]
ldv2 = read.csv("C2_polymer_fix_d067558_all_variants.csv" )
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
ggplot(b,aes(x=NumPasses,y=pr))+geom_point()+labs(x="Minimum Number of Passes Required", y = "Percentage of Errors Removed", title ="How Polishing Effect Varies by Coverage Cutoff")+theme_bw(base_size=8)+geom_hline(yintercept=1, color="red")
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
  1-res[2,2]/res[1,2]
}
k = Vectorize(errorContour)
x = seq(0,14,.5) # Min GC
y = seq(0,20) # Min Passes
z = outer(x,y,k)

pdf("ErrorsRemovedContourPolish2orMore.pdf",width=6.5,height=5)
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


minPasses=20
ok = cd$ZMW[cd$NumPasses>=minPasses]
gd = cdv[cdv$zmw%in%ok,] 
res = aggregate(Pos~Ref+type+homopolymerLength+homopolymerChar+indeltype,gd,length)

head(res)


res2 = res[res$homopolymerChar=="C" & res$Ref=="Original",]
head(res2)
pdf("LambdaIndelErrorsPolishedAllData.pdf",width=6,height=4)
v = ggplot(res2,aes(x=homopolymerLength,y=Pos, colour=Ref))+facet_grid( .~indeltype)+geom_point(size=3)+theme_bw(base_size=9)
v = v +labs(x="Homopolymer Length",y="Count of Total Errors", 
            title=expression(paste(lambda, " G/C indel errors divided by genomic context for all reads with >20X coverage")))
v = v + scale_colour_discrete(name="Method",labels = c("Original","Full Polish"))#+geom_vline(xintercept=3.5,colour="red")
v
dev.off()

head(res)




res3 = res2[res2$type=="Indel",]
head(res3)
aggregate(Pos~Ref,res3,sum)
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


