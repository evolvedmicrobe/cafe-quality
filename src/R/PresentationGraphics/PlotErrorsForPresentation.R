library(ggplot2)
#setwd("/Users/nigel/git/cafe-quality/data")
setwd("/Users/nigel/git/cafe-quality/NotTracked/master_full")
#Section to plot the variant data.
#d = read.csv("variants_nochange.csv")
d = read.csv("master_5ba5286_all_variants.csv")
head(d)
#how many errors are G indels?
ng = sum(d$type=="Indel" & d$homopolymerChar=="G" & d$indeltype == "Deletion" & d$homopolymerLength > 1)
ng/nrow(d)
ng = sum(d$type=="Indel" & d$indeltype == "Deletion")
ng/nrow(d)
ng = sum(d$type=="Indel")
ng/nrow(d)

cnts = aggregate(zmw~Ref+type,data=d,FUN=length)
#pdf()
head(cnts)
cnts = cnts[!(as.character(cnts$Ref)%in%c("lambda_NEB3011","SmrtBellSequence")),]

pdf("ErrorsByReference.pdf",width=10,height=7)
ggplot(data = cnts,aes(x=type,y=zmw, fill=type)) + geom_bar(stat="identity") + facet_grid( . ~ Ref) +
   labs(x="Error Type", y = "Count", title = "Errors by Reference") + theme_bw(base_size=14)+scale_fill_manual(name="", guide=FALSE, values=c("blue","black"))
dev.off()

hp = d[d$Ref=="HP.V1.02",]
cnts = aggregate(zmw~Pos+type, data=hp, FUN=length)
cnts = cnts[order(-cnts$zmw),]
head(cnts)
cnts$zmw= cnts$zmw/sum(cnts$zmw)
library(scales)
pdf("ErrorsByLocHP.pdf",width=10.5, height=5)
ggplot(data = cnts, aes(x=Pos,y=zmw, fill=type)) + geom_bar(stat="identity") + scale_y_continuous(labels = percent_format()) + scale_fill_discrete(name="Error Type")+
  labs(x="Position of Error", y = "Percentage of Total Errors", title = "HP.V1.02 Errors by Template Location") + theme_bw(base_size=14)
dev.off()

## Make plot of insertion/deletion frequeincy
cntType = aggregate(zmw~indeltype+Pos, data=hp, FUN=length)
rd = cntType[cntType$Pos%in%c(29,44,75,95),]
rd$Pos2=factor(rd$Pos,labels=c("5bp-G","10bp-G","5bp-T","10bp-T"), levels=c(29,44,95,75))
rd$Mutation = rd$indeltype:factor(rd$Pos2)
rd$Mutation = factor(rd$Mutation, levels=levels(rd$Mutation)[c(1,5,2,6,3,7,4,8)])
pdf("ErrorsForHPCounts.pdf",height=6,width=13)
ggplot(rd, aes(x=Mutation,y=zmw, fill=indeltype))+geom_bar(stat="identity")+theme_bw()+
  labs(x="Mutation:HP-Type",y="Error Counts", title="Counts of Indel Error Types for HP.V1.02 Homopolymers")+scale_fill_discrete(name="Indel Type")
dev.off()

# Make plot of relative sizes for insertions/deletions
rd = hp[hp$Pos%in%c(29,44,75) & hp$indeltype=="Deletion",]
rd$Pos = factor(rd$Pos)
ggplot(rd, aes(x=indelSize)) +
  geom_histogram(bin=1)+facet_grid( . ~ Pos) + theme_bw() +
  labs(x="Deletion Size", y= "Count")




cnts
15490/(15490+1325)
48834/(48834+684)
5120/(5120+603)

#rd = read.csv("ccs_nochange.csv")
rd = read.csv("master_5ba5286_combined_reads.csv")
rd = rd[rd$Reference!="SmrtBellSequence" & rd$Reference!="lambda_NEB3011",]
rd$Reference = factor(rd$Reference)
rd$ErrorRate = (rd$NumErrors/rd$Length)
rd$SNPErrorRate = (rd$NumSNPErrors / rd$Length)
rd$IndelErrorRate = (rd$NumIndelErrors / rd$Length)
tmp = rd[,c("Reference","IndelErrorRate","SNPErrorRate","NumPasses")]
vd = melt(tmp, id =c("Reference","NumPasses") )
head(vd)
fd = aggregate(value~Reference+NumPasses+variable, vd, mean)
head(fd)
fd$value=-10*log10(fd$value)
pdf("ErrorsByRefAndType.pdf", width=10, height=7)
ggplot(fd,aes(x=NumPasses , y = value, colour=Reference, shape = variable, linetype = variable)) +
  geom_line()+geom_point()+theme_bw(base_size=14)+ scale_shape_discrete(name="Error Type", labels=c("Indel", "SNP")) +
  labs(x="Number of Subreads",y = "Mean Empirical Error Rate (Phred QV Scaled)", title="Errors by Reference/Subreads")+
  scale_x_continuous(limits=c(0,75)) + scale_linetype_discrete(name="Error Type", labels=c("Indel", "SNP"))
dev.off()

head(rd)
rd$MeanGC = (rd$SnrG + rd$SnrC) / 2

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


ggplot(rd,aes(x=PredictedCCSAccuracy,y=NumPasses, colour=Reference))+geom_point()+geom_smooth()+scale_x_continuous(limits=c(.85,1.0))



merrors = aggregate(ErrorRate~NumSubReads+Reference,data=rd[rd$NumSubReads > 2,],FUN=mean)
merrors$QV = -10 * log10(merrors$ErrorRate)


pdf("ErrorsByRef.pdf", width=10, height=7)
ggplot(merrors,aes(x=NumSubReads , y = QV, group=Reference, colour = Reference)) +
  geom_line()+geom_point()+theme_bw(base_size=14)+
  labs(x="Number of Subreads",y = "Mean Empirical Error Rate (Phred QV Scaled)", title="Errors by Reference/Subreads")+ scale_x_continuous(limits=c(0,75))
dev.off()


hp = rd[rd$Reference=="HP.V1.02",]
h2 = hp[hp$NumSubReads==100 ,]
mean(h2$NumErrors)
head(h2)

hist(hp$NumSubReads,450,col="blue",main="Subread distribution for HP Reference", xlab = "Number of SubReads")
tops =aggregate(ZMW~NumSubReads,data=hp,FUN=length)
max(tops$ZMW)

d136 = hp[hp$NumSubReads==136,]
hist(d136$NumErrors)

ggplot(data=rd, aes(x=NumSubReads)) + 
  geom_histogram(binwidth=5) + 
  facet_grid(.~Reference) + theme_bw() +
  title("Distributions of Subread Counts") + labs(x="Subread count")

#hp = rd[rd$Reference=="ALL4MER.V2.01",]
hp = rd[rd$Reference=="HP.V1.02",]

aggregate(ZMW~NumSubReads,data=hp,FUN=length)

hist(hp$NumSubReads,500)

head(hp)
ggplot(hp, aes(y=NumIndelErrors,x=NumSubReads)) + 
  geom_jitter() + stat_smooth() + labs(x="Sub Read Count",y="Mean Indel Error Count"
