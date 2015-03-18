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
  d=d[d$Reference!="lambda_NEB3011",]
  #d=d[d$Reference=="HP.V1.02",]
  #d = d[d$Reference=="ALL4MER.V2.01",]
  #d = d[d$Referend=="NOHP.V1.01"]
  d$Reference=factor(d$Reference)
  return(d)
}

d = getld("master_5ba5286_combined_reads.csv")
#d2 = getld("C2_polymer_d067558_combined_reads.csv")
#d2 = getld("del_tag_combined_reads.csv")
#d3 = getld("no_del_combined_reads.csv")
#d2 = getld("16s_go2_combined_reads.csv")
d2 = getld("round1_read_combined_reads.csv")

q = merge(d, d2,by="ZMW")
blacklist = q$ZMW[q$Reference.x!=q$Reference.y]
sum(as.character(q$Reference.x)!=q$Reference.y)
#pdf("GvCSnr.pdf", height=6, width =6)
#ggplot(d[sample(size=5000, x=nrow(d)),], aes(x=SnrC, y = SnrG)) + geom_hline(yintercept=mean(d$SnrG), col="red", lwd = 1.5) + geom_vline(xintercept=mean(d$SnrC), col="red", lwd=1.5) +
#  geom_point() +geom_smooth() +theme_bw() + 
#  labs(x="C channel SNR", y = "G channel SNR", title="Relationship between G and C SNR")
#dev.off()  

in1 = intersect(d$ZMW,d2$ZMW)
in1=in1[!(in1%in%blacklist)]
d2 = d2[d2$ZMW%in%in1,]
d = d[d$ZMW%in%in1,]

cd = rbind(d,d2)
cd$Analysis =  factor(c(rep("Original",nrow(d)),rep("Phase 1",nrow(d2))))
ph = function(x) { 
  if(x==1) {0} else if(x==0) {60} else { -10*log10(x) }
}
cd$PredQV = sapply((1-cd$PredictedCCSAccuracy),ph)
cd$ActualQV = sapply(cd$ErrorRate,ph)

## Make a quick graphic
pdf("errorsWithCoverage.pdf", width=8, height=6)
nd = d[d$Reference=="HP.V1.02",]
res = aggregate(ErrorRate~NumPasses, nd, FUN=mean)
res$QV = sapply(res$ErrorRate, ph)
ggplot(res[res$NumPasses>3,], aes(x=NumPasses, y=QV))+geom_line()+
  theme_bw(base_size=12)+labs(x="Number of Passes", y ="Empirical Error Rate (Phred Scaled)", title="Quality by Coverage in HP.V1.02") 
head(nd)
dev.off()

# Get the total error difference.
aggregate(NumErrors~Analysis+Reference, data=cd, FUN=sum)

aggregate(NumErrors~Analysis+Reference, data=cd[cd$NumPasses>20,], FUN=sum)

aggregate(NumErrors~Analysis+Reference, data=cd[cd$NumPasses>20 & cd$PredictedRawAccuracy > .84,], FUN=sum)


ggplot(cd,aes(x=PredQV, y=ActualQV, colour=Analysis)) + geom_point() + geom_smooth()
ggplot(cd, aes(x=PredictedCCSAccuracy, fill=Analysis))  + geom_density(alpha=.2) + scale_x_continuous(limits=c(0.99, 1.0))

res$QV = sapply(res$ErrorRate,ph)

#pdf("ErrorHistogram.pdf")
cutoff = 20
ggplot(d2[d2$NumPasses>cutoff & d2$ErrorRate>0.0,], aes(x=ErrorRate)) + geom_histogram() + theme_bw(base_size=12) + labs(x="CCS Empirical Error Rate", y = "Count", title=paste("Histogram of Error Rate With >", cutoff, "Passes"))
#dev.off()


toUse = cd$ZMW[cd$PredictedCCSAccuracy > .99 && cd$Analysis=="Original"]
toUse = cd$ZMW[cd$PredictedRawAccuracy > 0.85 & cd$Reference=="NOHP.V1.01"]


res = aggregate(ErrorRate~NumPasses+Analysis+Reference, cd[cd$PredictedCCSAccuracy > .99,], FUN=mean)
res = aggregate(ErrorRate~NumPasses+Analysis+Reference, cd[cd$ZMW%in%toUse,], FUN=mean)
#res = aggregate(ErrorRate~NumPasses+Analysis+Reference, cd[cd$ZMW%in%toUse,], FUN=length)
res$QV = sapply(res$ErrorRate, ph)
head(res)

pdf("NOHPbyCoverage.pdf")
v = ggplot(res,aes(x=NumPasses, y = QV, color=Analysis, fill=Analysis))+geom_line()+geom_point()+theme_bw(base_size=12)+scale_x_continuous(limits=c(0,60)) 
v = v+labs(x="Number of Passes", y = "Per Base Error Rate (Phred QV Scale)", title="NOHP.V1 Accuracy by Coverage (Raw Accuracy > 0.85)")
v = v + scale_color_manual(name="Analysis", values=c("blue","red")) + scale_fill_manual(name="Analysis", values=c("blue","red"))  + scale_y_continuous(limits=c(12,60)) +geom_smooth()
v
dev.off()

quantile(d)

geom_vline(xintercept=10, colour="red") +geom_hline(yintercept=30, colour="red")

v = ggplot(res[res$Reference=="ALL4MER.V2.01",],aes(x=NumPasses, y = QV, color=Analysis))+geom_line()+geom_point()+theme_classic(base_size=11)+scale_x_continuous(limits=c(0,25)) 
v = v+labs(x="Number of Passes", y = "Per Base Error Rate (Phred QV Scale)")
v = v +scale_colour_discrete(name="Method") + scale_y_continuous(limits=c(15,45))
v


org = 10^(-28.66406/10)
new = 10^-(42.50566/10)
new/org


qplot(log10(1-d2$PredictedCCSAccuracy), log10(d2$ErrorRate))
qplot(d2$MinSNR,d2$MeanSNR)
ggplot(d2,aes(x=MinSNR,fill=Reference))+geom_density(alpha=.2)+geom_vline(xintercept=3.5, col="red")
ggplot(d2,aes(x=MeanGC,fill=Reference))+geom_density(alpha=.2)+geom_vline(xintercept=3.5, col="red")

dq2 = aggregate(IndelErrorRate~NumPasses+RoundGC, d2, mean )
ph = function(x) { -10*log10(x)}
dq2$IndelErrorRate = ph(dq2$IndelErrorRate)
v = ggplot(dq2[dq2$NumPasses>1,],aes(x=NumPasses, y = IndelErrorRate, colour=RoundGC)) + geom_line()+ theme_bw(base_size=8)
v = v + labs(x= "Number of Passes", y ="Indel Error Rate (Phred Scaled)")
v = v + scale_x_continuous(limits=c(0,100))
v




#pdf("ErrorInAll4Mer.pdf")
v = ggplot(res,aes(x=NumPasses, y = QV, colour=Analysis))+geom_line()+geom_point()+theme_classic(base_size=11)+scale_x_continuous(limits=c(0,60)) 
v = v+labs(x="Number of Passes", y = "Per Base Error Rate (Phred QV Scale)",title="Mean accuracy in ALL4MERS ")
v = v + geom_vline(xintercept=10, colour="red") +geom_hline(yintercept=30, colour="red")
v
#dev.off()



v = ggplot(res,aes(x=NumPasses, y = QV, colour=Analysis))+geom_line()+geom_point()+theme_classic(base_size=11)+scale_x_continuous(limits=c(0,60)) 
v = v+labs(x="Number of Passes", y = "Per Base Error Rate (Phred QV Scale)",title=expression(paste("Mean accuracy in ",lambda, " with > 0.85 Raw Accuracy")))
v = v + geom_vline(xintercept=10, colour="red") +geom_hline(yintercept=30, colour="red")
v


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
ct = "ALL4MER.V2.01"
ct = "HP.V1.02"
ct = "NOHP.V1.01"
dv = read.csv("master_5ba5286_all_variants.csv")
ldv = dv[dv$Ref==ct,]
ldv2 = read.csv("round1b_all_variants.csv" )
ldv2 = ldv2[ldv2$Ref==ct,]

cut1 = quantile(d$PredictedCCSAccuracy, .2)
cut2 = quantile(d2$PredictedCCSAccuracy, .2)
#ldv= ldv[ldv$zmw%in%(d$ZMW[d$PredictedCCSAccuracy >cut1]),]
#ldv2= ldv2[ldv2$zmw%in%(d2$ZMW[d2$PredictedCCSAccuracy >cut2]),]

cdv = rbind(ldv,ldv2)
cdv$Analysis = c(rep("Original",nrow(ldv)),rep("Phase 1",nrow(ldv2)))
valid = cd$ZMW[cd$NumPasses > 10 & cd$PredictedRawAccuracy > 0.85]
cdv = cdv[cdv$zmw%in%valid,]
cdv$homopolymerChar[cdv$homopolymerChar=="T"]="A"
cdv$homopolymerChar[cdv$homopolymerChar=="G"]="C"

ldv2 = ldv2[ldv2$zmw%in%valid,]
aggregate(Ref~QV+type, data = ldv2, FUN=length)
aggregate(Ref~type, data = ldv2, FUN=length)


errors = aggregate(Ref~Analysis+type, data =cdv, FUN=length)
#pdf("ErrorsInALL4Mers.pdf", width=8, height=6)
ggplot(errors,aes(x=Analysis, y =Ref, fill=Analysis)) + geom_bar(stat="identity") +facet_grid(.~type) +
  scale_fill_manual(values=c("blue","red")) +theme_bw(base_size=14) +labs(y="Error Count", title="ALL4MERS Performance Comparison (Top 80% of data)")
#dev.off()

spots = aggregate(Ref~Analysis+type+indelSize+indeltype+homopolymerChar, data=cdv, FUN=length) 
ggplot(spots, aes(x=indelSize, y=Ref, fill=Analysis))+geom_bar(stat="identity", position="dodge") + 
  labs(x="Template Position", y="Error Count", title="NOHP Errors by Position (Raw Accuracy > 0.85, Coverage > 10)") +
  theme_bw(base_size=14) + scale_fill_manual(values=c("blue", "red")) + facet_wrap(indeltype~homopolymerChar)




errors2 = aggregate(Ref~Analysis+Pos, data =cdv, FUN=length)
#pdf("ErrorPositionDistribution.pdf",width=11, height=6)
ggplot(errors2, aes(x=Pos, y=Ref, fill=Analysis))+geom_bar(stat="identity", position="dodge") + 
  labs(x="Template Position", y="Error Count", title="ALL4MERS Errors by Position (Raw Accuracy > 0.85, Coverage > 15)") +
  theme_bw(base_size=14) + scale_fill_manual(values=c("blue", "red"))
#dev.off()

bp= aggregate(indelSize~Pos+Analysis, data = cdv, FUN=length)

# compare original and current
ggplot(bp, aes(x=Pos, y=indelSize, fill=Analysis)) + geom_bar(stat="identity", position="dodge")

#Now do that for over 120X
toUse = cd$ZMW[cd$NumPasses >=20 & cd$PredictedRawAccuracy > 0.85 & cd$PredictedCCSAccuracy > .999 & cd$Analysis=="Phase 1"]
bp2= aggregate(indelSize~Pos+Analysis, data = cdv[cdv$zmw%in%toUse,], FUN=length)
bp2= bp2[bp2$Pos%in%c(29,44,75,95),]
bp2$BP = factor(rep(c("G", "G", "A", "A"),2))
bp2$Length = factor(rep(c("5", "10", "10", "5"), 2))
bp2$Region = factor(bp2$Length:bp2$BP, levels=c("5:G", "10:G", "5:A", "10:A"))
#pdf("ImprovementsWithPhase1.pdf", width=9, height=6)
ggplot(bp2, aes(x=Region, y=indelSize, fill=Analysis)) + geom_bar(stat="identity", position="dodge") +
  labs(x="HP Type (Length:BP)", y = "Error Counts", title="Homopolymer error counts from CCS reads with > 20 Passes") + theme_bw(base_size=14) + scale_fill_manual(values=c( "blue","red"))
#dev.off()




v = ggplot(res2,aes(x=homopolymerLength,y=Pos, colour=Analysis, shape=homopolymerChar))+facet_grid( .~indeltype)+geom_point(size=3)+theme_bw(base_size=9)
v = v +labs(x="Homopolymer Length",y="Count of Total Errors", 
            title=expression(paste(lambda, " G/C indel errors divided by genomic context for all Reads")))
v = v + scale_colour_discrete(name="Method")#+geom_vline(xintercept=3.5,colour="red")
v

head(res)


# Let's compare for the 

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


