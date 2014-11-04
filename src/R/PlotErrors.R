library(ggplot2)
setwd("/Users/nigel/git/cafe-quality/data")

#Section to plot the variant data.
d = read.csv("variants.csv")
head(d)
#how many errors are G indels?
ng = sum(d$type=="Indel" & d$homopolymerChar=="G" & d$indeltype == "Deletion" & d$homopolymerLength > 1)
ng/nrow(d)
ng = sum(d$type=="Indel" & d$indeltype == "Deletion")
ng/nrow(d)
ng = sum(d$type=="Indel")
ng/nrow(d)

cnts = aggregate(zmw~Ref+type,data=d,FUN=length)
ggplot(data = cnts,aes(x=type,y=zmw)) + geom_bar(stat="identity") + facet_grid( . ~ Ref) +
   labs(x="Error Type", y = "Count", title = "Errors by Reference") + theme_bw()


hp = d[d$Ref=="HP.V1.02",]
cnts = aggregate(zmw~Pos+type, data=hp, FUN=length)
cnts = cnts[order(-cnts$zmw),]
head(cnts)
ggplot(data = cnts, aes(x=Pos,y=zmw)) + geom_bar(stat="identity") +
  labs(x="Position of Error on Refernce", y = "Count", title = "Error Count by Reference Location") + theme_bw()

## Make plot of insertion/deletion frequeincy
cntType = aggregate(zmw~indeltype+Pos, data=hp, FUN=length)
rd = cntType[cntType$Pos%in%c(29,44,75),]
rd$Mutation = rd$indeltype:factor(rd$Pos)
ggplot(rd, aes(x=Mutation,y=zmw))+geom_bar(stat="identity")+theme_bw()+labs(x="Mutation:Position",y="Counts")

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

rd = read.csv("ccs.csv")

rd = rd[rd$Reference!="SmrtBellSequence",]
rd$Reference = factor(rd$Reference)
merrors = aggregate(NumErrors~NumSubReads+Reference,data=rd[rd$NumSubReads > 2,],FUN=mean)



ggplot(merrors,aes(x=NumSubReads , y = NumErrors, group=Reference, colour = Reference)) +
  geom_line()+geom_point()+theme_classic(base_size=14)+
  labs(x="Number of Sub Reads",y = "Mean Number of Errors", title="Errors by Reference/Subreads")

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
