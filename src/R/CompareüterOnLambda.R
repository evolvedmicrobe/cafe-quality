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
  snr=d[,grep("Snr",colnames(d))]
  d$MinSNR = apply(snr,1,min)
  d$MeanSNR = apply(snr,1,mean)
  d$RoundGC = factor(round(4*floor(d$MeanSNR/4)))
  d$RoundRQ = factor(round(5*floor(100*d$PredictedRawAccuracy/5))+2.5)
  s = nrow(d)
  #d = d[d$PredictedRawAccuracy > .85,]
  e = nrow(d)
  print(e/s)
  d=d[d$Reference=="lambda_NEB3011",]
  #d=d[d$Reference=="HP.V1.02",]
  #d = d[d$Reference=="ALL4MER.V2.01",]
  return(d)
}

d = getld("master_5ba5286_combined_reads.csv")
d2 = getld("lambda_read_combined_reads.csv")
#d2 = getld("del_tag_combined_reads.csv")
#d3 = getld("no_del_combined_reads.csv")
#d2 = getld("16s_go2_combined_reads.csv")
#d3 = getld("round1_read_combined_reads.csv")

in1 = intersect(d$ZMW,d2$ZMW)
d2 = d2[d2$ZMW%in%in1,]
d = d[d$ZMW%in%in1,]

cd = rbind(d,d2)
cd$Analysis =  factor(c(rep("Original",nrow(d)),rep("Phase 1",nrow(d2))))


# Get the total error difference.
aggregate(NumErrors~Analysis, data=cd, FUN=sum)




ph = function(x) { 
  if(x==1) {0} else if(x==0) {60} else { -10*log10(x) }
}


## Make a quick graphic
pdf("LambdaerrorsWithCoverage.pdf", width=4, height=3)
#nd = d[d$Reference=="HP.V1.02",]
toUse = cd$ZMW[cd$PredictedCCSAccuracy > .99 & cd$Analysis =="Phase 1" & cd$PredictedRawAccuracy > 0.85]
nd = cd[cd$ZMW%in%toUse,]
res = aggregate(ErrorRate~NumPasses+Analysis, nd, FUN=mean)
res$QV = sapply(res$ErrorRate, ph)
ggplot(res[res$NumPasses>3,], aes(x=NumPasses, y=QV, color=Analysis))+geom_line()+
  theme_bw(base_size=8)+labs(x="Number of Passes", y ="Empirical Error Rate (Phred Scaled)", title="Quality by Coverage in Lambda\n(Raw > 0.85, CCS > .99)") 
head(nd)
dev.off()


cd$PredQV = sapply((1-cd$PredictedCCSAccuracy),ph)
cd$ActualQV = sapply(cd$ErrorRate,ph)

ggplot(cd,aes(x=PredQV, y=ActualQV, colour=Analysis)) + geom_point() + geom_smooth()

ggplot(cd, aes(x=PredictedCCSAccuracy, fill=Analysis))  + geom_density(alpha=.2) + scale_x_continuous(limits=c(0.99, 1.0))

res$QV = sapply(res$ErrorRate,ph)

#pdf("ErrorHistogram.pdf")
ggplot(d2[d2$NumPasses>8 & d2$ErrorRate>0.0,], aes(x=ErrorRate)) + geom_histogram() + theme_bw(base_size=12) + labs(x="CCS Empirical Error Rate", y = "Count", title="Histogram of Error Rate With > 8 Passes")
#dev.off()


toUse = cd$ZMW[cd$PredictedCCSAccuracy > .99 & cd$Analysis=="Original"]
toUse = cd$ZMW[cd$PredictedCCSAccuracy > .99 & cd$Analysis=="Phase 1"]

toUse = cd$ZMW[cd$PredictedRawAccuracy > 0.85 cd$PredictedCCSAccuracy > .99 & cd$Analysis=="Phase 1"]


res = aggregate(ErrorRate~Analysis, cd[cd$PredictedCCSAccuracy > .999,], FUN=length)
res = aggregate(ErrorRate~NumPasses+Analysis, cd[cd$ZMW%in%toUse,], FUN=mean)
res$QV = sapply(res$ErrorRate, ph)
head(res)


v = ggplot(res,aes(x=NumPasses, y = QV, colour=Analysis))+geom_line()+geom_point()+theme_classic(base_size=11)#+scale_x_continuous(limits=c(0,25)) 
v = v+labs(x="Number of Passes", y = "Per Base Error Rate (Phred QV Scale)")
v = v + geom_vline(xintercept=10, colour="red") +geom_hline(yintercept=30, colour="red")+scale_colour_discrete(name="Binned Raw\nQuality")
v



qplot(d2$MinSNR,d2$MeanSNR)
ggplot(d2,aes(x=MinSNR,fill=Reference))+geom_density(alpha=.2)+geom_vline(xintercept=3.5, col="red")
ggplot(d2,aes(x=MeanGC,fill=Reference))+geom_density(alpha=.2)+geom_vline(xintercept=3.5, col="red")

dq2 = aggregate(IndelErrorRate~NumPasses+RoundGC, d2, mean )
ph = function(x) { -10*log10(x)}
dq2$IndelErrorRate = ph(dq2$IndelErrorRate)
pdf("ErrorBySNR.pdf")
v = ggplot(dq2[dq2$NumPasses>1,],aes(x=NumPasses, y = IndelErrorRate, colour=RoundGC)) + geom_line()+ theme_bw(base_size=8)
v = v + labs(x= "Number of Passes", y ="Indel Error Rate (Phred Scaled)")
v = v + scale_x_continuous(limits=c(0,100))
v
dev.off()
names(q2)[3]<-"ErrorRate"
t2 = rep("Indel",nrow(q2))




#pdf("ErrorInAll4Mer.pdf")
v = ggplot(res,aes(x=NumPasses, y = QV, colour=Analysis))+geom_line()+geom_point()+theme_classic(base_size=11)+scale_x_continuous(limits=c(0,60)) 
v = v+labs(x="Number of Passes", y = "Per Base Error Rate (Phred QV Scale)",title="Mean accuracy in ALL4MERS ")
v = v + geom_vline(xintercept=10, colour="red") +geom_hline(yintercept=30, colour="red")
v
#dev.off()


pdf("ErrorInLambda.pdf")
v = ggplot(res,aes(x=NumPasses, y = QV, colour=Analysis))+geom_line()+geom_point()+theme_classic(base_size=11)+scale_x_continuous(limits=c(0,60)) 
v = v+labs(x="Number of Passes", y = "Per Base Error Rate (Phred QV Scale)",title=expression(paste("Mean accuracy in ",lambda, " with > 0.85 Raw Accuracy")))
v = v + geom_vline(xintercept=10, colour="red") +geom_hline(yintercept=30, colour="red")
v
dev.off()


res2 = aggregate(NumErrors~NumPasses+Analysis, cd, FUN=sum)
v = ggplot(res2,aes(x=NumPasses, y = NumErrors, colour=Analysis))+geom_line()+geom_point()+theme_classic(base_size=11)+scale_x_continuous(limits=c(0,60)) 
v = v+labs(x="Number of Passes", y = "Per Base Error Rate (Phred QV Scale)",title=expression(paste("Mean accuracy in ",lambda)))
v = v + scale_x_continuous(limits= c(10,20)) + scale_y_continuous(limits=c(0,1000))
v


cds = cd[cd$NumPasses==10,]
cds$QV = ph(cds$ErrorRate)
cds$QV[is.infinite(cds$QV)]=40

plot(cds$QV,cds$MeanGC)
ggplot(cds,aes(x=QV,fill=Analysis))+geom_histogram(alpha=.5,position="dodge")
ggplot(cds,aes(x=QV,colour=Analysis))+geom_freqpoly(alpha=.5)

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
dv = read.csv("master_variants_all_variants.csv")
ldv = dv[dv$Ref=="lambda_NEB3011",]
ldv2 = read.csv("lambda_all_variants.csv" )
ldv2 = ldv2[ldv2$Ref=="lambda_NEB3011",]

cdv = rbind(ldv,ldv2)
cdv$Analysis = c(rep("Original",nrow(ldv)),rep("Phase 1",nrow(ldv2)))
cdv = cdv[cdv$zmw%in%in1,]
cdv$homopolymerChar[cdv$homopolymerChar=="T"]="A"
cdv$homopolymerChar[cdv$homopolymerChar=="G"]="C"

aggregate(Ref~Analysis+type, data =cdv, FUN=length)

aggregate(Ref~Analysis+type+indelSize, data=cdv, FUN=length)



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
ggplot(b,aes(x=NumPasses,y=pr))+geom_point()+labs(x="Minimum Number of Passes Required", y = "Percentage of Errors Removed", title ="How Polishing Effect Varies by Coverage Cutoff")+theme_bw(base_size=8)


errorDrop2 <-function(minGC) {
  ok = cd$ZMW[cd$MeanGC>=minGC & cd$Mean]
  gd = cdv[cdv$zmw%in%ok,] 
  head(gd)
  res = aggregate(Pos~Ref,gd,length)
  print(res)
  res[2,2]/res[1,2]
}

s=seq(7,15)
y2 = sapply(s,errorDrop2)
b=data.frame(NumPasses = s, pr = y2)

pdf(height=4,width=4,"ImprovementsBySNR.pdf")
ggplot(b,aes(x=NumPasses,y=pr))+geom_point()+labs(x="Minimum Number of Passes Required",
                                                  y = "Percentage of Errors Removed",
                                                  title ="How Polishing Effect Varies by GC SNR Cutoff")+theme_bw()
dev.off()


### CONTOUR PLOT
errorContour <-function(minGC, minPasses) {
  ok = cd$ZMW[cd$MeanGC>=(minGC - 1) & cd$MeanGC<=(minGC + 1) & cd$NumPasses>=(minPasses-1) & cd$NumPasses<=(minPasses+1)]
  print(minGC)
  print(minPasses)
  gd = cdv[cdv$zmw%in%ok,]
  print(nrow(gd))
  
  res = aggregate(Pos~Analysis,gd,length)
  print(res)
  1-res[2,2]/res[1,2]
}
k = Vectorize(errorContour)
x = seq(7,13,.5) # Min GC
y = seq(2,20) # Min Passes
z = outer(x,y,k)

pdf("ErrorsRemovedContourPolish2orMore.pdf",width=6.5,height=5)
filled.contour(x,y,z, color.palette=rainbow, xlab="Mean G+C SNR",ylab="Minimum Number of Passes", main="% Errors Removed By SNR and Coverage")
dev.off()

pdf("LogScaleErrors.pdf", width=10.5, height=7.5)
par(mfrow=c(2,2))
errorCount <-function(minGC, minPasses) {
  ok = d[d$MeanGC>=(minGC - 2) & d$MeanGC<=(minGC + 2) & d$NumPasses>(minPasses-2) & d$NumPasses<(minPasses+2),]
  sum(ok$NumIndelErrors)
}
k = Vectorize(errorCount)
x = seq(1,14,.5) # Min GC
y = seq(1,20) # Min Passes
z = outer(x,y,k)
filled.contour(x,y,z,color.palette=rainbow, xlab="Mean G+C SNR",
               ylab="Minimum Number of Passes", 
               main="Total Indel Errors by\nSNR and Coverage", 
               key.title=title(main="Count"))

filled.contour(x,y,log10(z),color.palette=rainbow, xlab="Mean G+C SNR",
               ylab="Minimum Number of Passes", 
               main="Log Total Indel Errors by\nSNR and Coverage", 
               key.title=title(main="Log10\nCount"))


readCount <-function(minGC, minPasses) {
  sum(d$Length[d$MeanGC>=(minGC - 2) & d$MeanGC<=(minGC + 2) & d$NumPasses>=(minPasses-1) & d$NumPasses<=(minPasses+1) & d$Length<2500])
}
k = Vectorize(readCount)
x = seq(1,14,.5) # Min GC
y = seq(1,20) # Min Passes
z = outer(x,y,k)

filled.contour(x,y,z,color.palette=rainbow, xlab="Mean G+C SNR",
               ylab="Minimum Number of Passes", 
               main="Total CCS Basepairs", 
               key.title=title(main="BP\nCount"))

errorRate <-function(minGC, minPasses) {
  ok = d[d$MeanGC>=(minGC - 2) & d$MeanGC<=(minGC + 2) & d$NumPasses>=(minPasses-1) & d$NumPasses<=(minPasses+1) & d$Length<2500,]
  mean(ok$IndelErrorRate)
}
k = Vectorize(errorRate)
x = seq(3,14,.5) # Min GC
y = seq(1,20) # Min Passes
z = outer(x,y,k)

filled.contour(x,y,log10(z),color.palette=rainbow, xlab="Mean G+C SNR",
               ylab="Minimum Number of Passes", 
               main="CCS Log10 Indel Error Rate")

dev.off()


minPasses=15
ok = cd$ZMW[cd$NumPasses>=minPasses]
gd = cdv[cdv$zmw%in%ok & cdv$QV>40,] 
res = aggregate(Pos~Analysis,gd,length)

res = aggregate(Pos~Analysis+type+homopolymerLength+homopolymerChar+indeltype,gd,length)
head(res)

res2 = res[res$homopolymerChar=="A",]
res2= res
head(res2)

pdf("LambdaByContext.pdf", width=7, height = 4)
v = ggplot(res2,aes(x=homopolymerLength,y=Pos, colour=Analysis, shape=homopolymerChar))+facet_grid( .~indeltype)+geom_point(size=3)+theme_bw(base_size=8)
v = v +labs(x="Homopolymer Length",y="Count of Total Errors", 
            title=expression(paste(lambda, " Indel errors divided by genomic context (Passes > 15, Only most confidently called variants)")))
v = v + scale_colour_discrete(name="Method")#+geom_vline(xintercept=3.5,colour="red")
v
dev.off()
ggplot(cdv, aes(x=QV, fill=Analysis)) + geom_density(alpha=.3)

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


