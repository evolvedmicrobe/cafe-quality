library(ggplot2)
library(gridExtra)
library(grid)

setwd("/Users/nigel/git/cafe-quality/NotTracked/iqvComparison/")

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
  #d = d[d$Referend=="NOHP.V1.01"]
  d$Reference=factor(d$Reference)
  return(d)
}

#d = getld("BaumWelch_0dcd9a6_combined_reads.csv")
d = getld("master_5ba5286_combined_reads.csv")
#d2 = getld("C2_polymer_d067558_combined_reads.csv")
#d2 = getld("del_tag_combined_reads.csv")
#d3 = getld("no_del_combined_reads.csv")
#d2 = getld("16s_go2_combined_reads.csv")
#d2 = getld("lambdaParams_28fbe97_combined_reads.csv")
d2 = getld("lambdaParams_0abe1b2_combined_reads.csv")
d3 = getld("iqvlambda_cbacb93_combined_reads.csv")
#d4 = getld("iqvlambda_uniform_combined_reads.csv")
#d5 = getld("iqvlambda_oldparams_combined_reads.csv")
#d4 = getld("mergeQV_937384c_combined_reads.csv")
d4 = getld("pulsewidth_79b8946_combined_reads.csv")
d5 = getld("mergeQV_3b6e9b4_combined_reads.csv")
head(d5)

d5$Accepted = with(d5, NumberAcceptedMutations / NumberTriedMutations)

pdf("percentageAccepted.pdf")
ggplot(d5, aes(x=Accepted)) + geom_density(fill="blue") + scale_x_continuous(label=percent) + theme_bw() +labs(x="Percentage of accepted mutations out of all tested", title="Lambda 2K Example")
dev.off()


d6 = getld("mergeQVDelTag_c0535b6_combined_reads.csv")
d4$Method="New"
d5$Method = "Old"

gd = rbind(d4,d5)
aggregate(ProcessingTime~Method, gd, sum)

ggplot(gd, aes(x=I(ProcessingTime), fill=Method)) + 
  geom_density(alpha = .75) + theme_bw(base_size=12) 
ggplot(gd, aes(x=ProcessingTime, y=Length, color=Method)) + geom_point(alpha = .1)
ggplot(gd, aes(x=NumPasses, y=ProcessingTime, color=Method)) + geom_point(alpha = .1)
ggplot(gd, aes(x=NumberTriedMutations, y=ProcessingTime, color=Method)) + geom_point(alpha = .1)



pdf("RuntimeImprovements.pdf")
ggplot(gd, aes(x=I(log10(ProcessingTime)), fill=Method)) + 
  geom_density(alpha = .75) + theme_bw(base_size=12) +
  labs(x="Log10 Processing Time (seconds)", title="Lambda-2K per ZMW Processing Time Comparison")
dev.off()


q = merge(d, d2, by="ZMW")
q2 = merge(q, d3, by="ZMW")
blacklist = q2$ZMW[q2$Reference.x!=q2$Reference.y | q2$Reference.x!=q2$Reference]
sum(as.character(q$Reference.x)!=q$Reference.y)
#pdf("GvCSnr.pdf", height=6, width =6)
#ggplot(d[sample(size=5000, x=nrow(d)),], aes(x=SnrC, y = SnrG)) + geom_hline(yintercept=mean(d$SnrG), col="red", lwd = 1.5) + geom_vline(xintercept=mean(d$SnrC), col="red", lwd=1.5) +
#  geom_point() +geom_smooth() +theme_bw() + 
#  labs(x="C channel SNR", y = "G channel SNR", title="Relationship between G and C SNR")
#dev.off()  

in1 = intersect(d$ZMW,d2$ZMW)
in1 = intersect(in1, d3$ZMW)
in1=in1[!(in1%in%blacklist)]
d = d[d$ZMW%in%in1,]
d2 = d2[d2$ZMW%in%in1,]
d3 = d3[d3$ZMW%in%in1,]
d4 = d4[d4$ZMW%in%in1,]
d5 = d5[d5$ZMW%in%in1,]
d6 = d6[d6$ZMW%in%in1,]



hd = merge(d3, d2, by=c("ZMW")) # uniform = x, new = y
ggplot(hd[hd$NumPasses.x>20,], aes(x=NumErrors.x, y = NumErrors.y)) + geom_point()
ggplot(hd[hd$NumPasses.x>20,], aes(x=NumErrors.x-NumErrors.y)) + geom_density()
onedif = hd[hd$NumPasses.x>20 & ((hd$NumErrors.x - hd$NumErrors.y) == 1),]

head(onedif)


d2 = d2[,colnames(d2)[colnames(d2)%in%colnames(d)]]
d3 = d3[,colnames(d3)[colnames(d3)%in%colnames(d)]]
d4 = d4[,colnames(d4)[colnames(d4)%in%colnames(d)]]
d5 = d5[,colnames(d5)[colnames(d5)%in%colnames(d)]]
d6 = d6[,colnames(d6)[colnames(d6)%in%colnames(d)]]




cd = rbind(d,d2, d3, d4, d5, d6)
cd$Analysis =  factor(c(rep("Original CCS",nrow(d)),
                        rep("Unite-EM",nrow(d2)), 
                        rep("IQV New", nrow(d3)),
                        rep("Pulse-Width New", nrow(d4)),
                        rep("Merge-QV New", nrow(d5)),
                        rep("Merge-QV DT", nrow(d6))
                        ))
                        #rep("IQV-Uniform", nrow(d4)),
                        #rep("IQV-Old Params", nrow(d5))))
ph = function(x) { 
  if(x==1) {0} else if(x==0) {60} else { -10*log10(x) }
}
cd$PredQV = sapply((1-cd$PredictedCCSAccuracy),ph)
cd$ActualQV = sapply(cd$ErrorRate,ph)

aggregate(ActualQV~Analysis, cd, function(x) sum(x>=30))
aggregate(ActualQV~Analysis + NumPasses, cd, function(x) sum(x>=30) / length(x))
aggregate(ActualQV~Analysis, cd[cd$NumPasses<50 & cd$NumPasses>10,], function(x) sum(x>=30) / length(x))
aggregate(ActualQV~Analysis + NumPasses, cd, median)


ggplot(cd[cd$Analysis!="IQV",], aes(y=ActualQV, x=PredQV, color = Analysis)) + geom_smooth() + geom_abline()


## Make a quick graphic
#pdf("errorsWithCoverageLambdaMergePulseWidth.pdf", width=8, height=6)
#nd = d[d$Reference=="HP.V1.02",]
#cd = cd[cd$Reference=="lambda_NEB3011",]
toUse = cd$ZMW[cd$Analysis =="Original CCS"]
#toUse = cd$ZMW[cd$Analysis=="Pulse-Width New" & cd$PredictedCCSAccuracy >= .99]
toUse = cd$ZMW[cd$Analysis=="Unite-EM" & cd$PredictedCCSAccuracy >= .999]
toUse = cd$ZMW[cd$Analysis=="Merge-QV DT" & cd$PredictedCCSAccuracy >= .999]
#toUse = cd$ZMW[cd$Analysis =="Merge-QV New"]
toUse = cd$ZMW[cd$Analysis =="Original CCS" & cd$PredictedCCSAccuracy >= .999]
toUse = cd$ZMW[cd$Analysis =="Original CCS" & cd$ActualQV >= 20]


nd = cd[cd$ZMW%in%toUse,]

res = aggregate(ErrorRate~NumPasses+Analysis+Reference, nd, FUN=mean)
res$QV = sapply(res$ErrorRate, ph)
pdf("MergeDTCoverage.pdf", height=8, width=10)
ggplot(res[res$NumPasses>3,], aes(x=NumPasses, y=QV, color=Analysis)) + geom_smooth(fill="white", alpha=.1) + theme_bw(base_size=14)  + labs(x="Number of Passes", y ="Empirical Error Rate (Phred Scaled)", title="Quality by Coverage (Pred. QV > 20 in original)")  + facet_wrap(~Reference, scale="free")
dev.off()
res = aggregate(ErrorRate~NumPasses+Analysis+Reference, nd, FUN=median)


res$QV = sapply(res$ErrorRate, ph)

pdf("ZoomedIn.pdf", width=10, height=8)
ggplot(res[res$NumPasses>3,], aes(x=NumPasses, y=QV, color=Analysis) )+geom_line()+ geom_point() + theme_bw(base_size=14)  + labs(x="Number of Passes", y ="Empirical Error Rate (Phred Scaled)", title="Quality by Coverage in Lambda (All Data)") + facet_wrap(~Reference, scale="free")+ scale_x_continuous(limits=c(5,15)) + scale_y_continuous(limits=c(25,40)) 
 head(nd)
dev.off()



#pdf("QV_Distribution.pdf")
ggplot(nd[nd$NumPasses>15,], aes(x=ActualQV, fill=Analysis)) + geom_histogram(position="dodge") + theme_bw(base_size=14) + labs(x="Empirical QV distribution (NP > 15)")

#dev.off()


#pdf("errorsWithCoverageHQHP.pdf", width=8, height=6)
#nd = d[d$Reference=="HP.V1.02",]
#cd = cd[cd$Reference=="lambda_NEB3011",]
toUse = cd$ZMW[cd$PredictedCCSAccuracy >= .999 & cd$Analysis =="IQV-Model"]
toUse = cd$ZMW[cd$PredictedCCSAccuracy >= .999 & cd$Analysis =="Unite-EM"]
nd = cd[cd$ZMW%in%toUse,]
res = aggregate(ErrorRate~NumPasses+Analysis, nd, FUN=mean)
res$QV = sapply(res$ErrorRate, ph)
ggplot(res[res$NumPasses>3,], aes(x=NumPasses, y=QV, color=Analysis))+geom_line()+ geom_point() +
  theme_bw(base_size=14)  + labs(x="Number of Passes", y ="Empirical Error Rate (Phred Scaled)", title="Quality by Coverage in Lambda (Pred QV>=30)") 
head(nd)
#dev.off()

head(nd)
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

#pdf("NOHPbyCoverage.pdf")
v = ggplot(res,aes(x=NumPasses, y = QV, color=Analysis, fill=Analysis))+geom_line()+geom_point()+theme_bw(base_size=12)+scale_x_continuous(limits=c(0,60)) 
v = v+labs(x="Number of Passes", y = "Per Base Error Rate (Phred QV Scale)", title="NOHP.V1 Accuracy by Coverage (Raw Accuracy > 0.85)")
v = v + scale_color_manual(name="Analysis", values=c("blue","red", "green")) + scale_fill_manual(name="Analysis", values=c("blue","red", "green"))  + scale_y_continuous(limits=c(12,60)) +geom_smooth()
v
#dev.off()

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
#ct = "NOHP.V1.01"
#ct = "lambda_NEB3011"
dv = read.csv("master_variants_all_variants.csv")
ldv = dv[dv$Ref==ct,]
#ldv2 = read.csv("lambdaParams_0abe1b2_all_variants.csv" )
ldv2 = read.csv("mergeQV_3b6e9b4_all_variants.csv")
ldv2 = ldv2[ldv2$Ref==ct,]
#ldv3 = read.csv("lambdaParams_0abe1b2_all_variants.csv")
ldv3 = read.csv("mergeQVDelTag_c0535b6_all_variants.csv")
ldv3 = ldv3[ldv3$Ref==ct,]


cut1 = quantile(d$PredictedCCSAccuracy, .2)
cut2 = quantile(d2$PredictedCCSAccuracy, .2)
#ldv= ldv[ldv$zmw%in%(d$ZMW[d$PredictedCCSAccuracy >cut1]),]
#ldv2= ldv2[ldv2$zmw%in%(d2$ZMW[d2$PredictedCCSAccuracy >cut2]),]

cdv = rbind(ldv,ldv2, ldv3)
cdv$Analysis = factor(c(rep("Original",nrow(ldv)),rep("PW-MODEL",nrow(ldv2)), rep("Merge-DelTag", nrow(ldv3))))
valid = cd$ZMW[cd$NumPasses > 20 ]
#valid = cd$ZMW[SnrC > 9]
cdv = cdv[cdv$zmw%in%valid,]


#cdv = cdv[cdv$zmw%in%valid & cdv$QV > 40,]
cdv$homopolymerChar[cdv$homopolymerChar=="T"]="A"
cdv$homopolymerChar[cdv$homopolymerChar=="G"]="C"


aggregate(Ref~QV+Analysis, data = cdv, FUN=length)
aggregate(Ref~type+Analysis, data = cdv, FUN=length)


errors = aggregate(Ref~Analysis+type, data =cdv, FUN=length)
#pdf("ErrorsInLambda_15.pdf", width=8, height=6)
ggplot(errors,aes(x=Analysis, y =Ref, fill=Analysis)) + geom_bar(stat="identity") +facet_grid(.~type) +theme_bw(base_size=14) +labs(y="Error Count", title="Lambda Performance Comparison (NumP > 15)")
#dev.off()


spots = aggregate(Ref~Analysis+type+indelSize+indeltype+homopolymerChar, data=cdv, FUN=length) 
ggplot(spots, aes(x=indelSize, y=Ref, fill=Analysis))+geom_bar(stat="identity", position="dodge") + 
  labs(x="Template Position", y="Error Count", title="Lambda Indel Errors by Type( (NP > 15)") +
  theme_bw(base_size=14)+ facet_wrap(indeltype~homopolymerChar)

#pdf("ErrorTypes.pdf", width=8, height=6)
spots = aggregate(Ref~Analysis+homopolymerLength+indeltype+homopolymerChar, data=cdv, FUN=length) 
ggplot(spots, aes(x=homopolymerLength, y=I(Ref), fill=Analysis))+geom_bar(stat="identity", position="dodge") + 
  labs(x="HP Length", y="Error Count", title="HP Errors by Position (Coverage > 15)") +
  theme_bw(base_size=14)  + facet_wrap(indeltype~homopolymerChar)


spots = aggregate(Ref~Analysis+homopolymerLength, data=cdv, FUN=length) 
ggplot(spots, aes(x=homopolymerLength, y=Ref, fill=Analysis))+geom_bar(stat="identity", position="dodge") + 
  labs(x="HP Length", y="Error Count", title="HP Errors by Position (Coverage > 15)") +
  theme_bw(base_size=14)

#dev.off()

spots = aggregate(Ref~Analysis+homopolymerLength+indeltype+homopolymerChar, data=cdv, FUN=length) 
ggplot(spots, aes(x=homopolymerLength, y=Ref, fill=Analysis))+geom_bar(stat="identity", position="dodge") + 
  labs(x="HP Length", y="Error Count", title="HP Errors by Position (Coverage > 15)") +
  theme_bw(base_size=14)  + facet_wrap(indeltype~homopolymerChar)



spots = aggregate(Ref~Analysis+homopolymerLength+indeltype+homopolymerChar, data=cdv[cdv$QV>40,], FUN=length) 
ggplot(spots, aes(x=homopolymerLength, y=Ref, fill=Analysis))+geom_bar(stat="identity", position="dodge") + 
  labs(x="HP Length", y="Error Count", title="Lambbda Errors by Position (Coverage > 15)") +
  theme_bw(base_size=14)  + facet_wrap(indeltype~homopolymerChar)


errors2 = aggregate(Ref~Analysis+Pos, data =cdv, FUN=length)
#pdf("ErrorPositionDistribution.pdf",width=11, height=6)
ggplot(errors2, aes(x=Pos, y=Ref, fill=Analysis))+geom_bar(stat="identity", position="dodge") + 
  labs(x="Template Position", y="Error Count", title="Lambda Errors by Position (Coverage > 20)") +
  theme_bw(base_size=14) + scale_fill_manual(values=c("blue", "red"))
#dev.off()

ggplot(cdv, aes(x=Pos, fill=Analysis))+geom_histogram(position="dodge", binwidth=500) + 
  labs(x="Template Position", y="Error Count", title="Lambda Errors by Position (Coverage > 20)") +
  theme_bw(base_size=14) + scale_fill_manual(values=c("blue", "red"))

#Get positions
b = errors2[errors2$Analysis=="Original",]
b[order(b$Ref),]

bp= aggregate(indelSize~Pos+Analysis, data = cdv, FUN=length)

# compare original and current
ggplot(bp, aes(x=Pos, y=indelSize, fill=Analysis)) + geom_bar(stat="identity", position="dodge")

#Now do that for over 120X
toUse = cd$ZMW[cd$NumPasses >=20 & cd$PredictedRawAccuracy > 0.85 & cd$PredictedCCSAccuracy > .999 & cd$Analysis=="Phase 1"]
aggregate(indelSize~Analysis+QV, data = cdv[cdv$zmw%in%toUse,], FUN=length)
bp2= aggregate(indelSize~Pos+Analysis, data = cdv[cdv$zmw%in%toUse,], FUN=length)
bp2= bp2[bp2$Pos%in%c(29,44,75,95),]
bp2$BP = factor(rep(c("G", "G", "A", "A"),2))
bp2$Length = factor(rep(c("5", "10", "10", "5"), 2))
bp2$Region = factor(bp2$Length:bp2$BP, levels=c("5:G", "10:G", "5:A", "10:A"))
#pdf("ImprovementsWithPhase1.pdf", width=9, height=6)
ggplot(bp2, aes(x=Region, y=indelSize, fill=Analysis)) + geom_bar(stat="identity", position="dodge") +
  labs(x="HP Type (Length:BP)", y = "Error Counts", title="Homopolymer error counts from CCS reads with > 20 Passes") + theme_bw(base_size=14) + scale_fill_manual(values=c( "blue","red"))
#dev.off()


ggplot(cdv[cdv$zmw%in%toUse,], aes(x=QV, fill=Analysis)) + geom_density(alpha=.5)

aggregate(QV~Analysis, cdv, function(x) sum(x<=30) / length(x) )

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


