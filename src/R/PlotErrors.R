library(ggplot2)

setwd("/Users/nigel/git/cafe-quality/data")
d = read.csv("variants.csv")
head(d)

cnts = aggregate(zmw~Ref+type,data=d,FUN=length)

ggplot(data = cnts,aes(x=type,y=zmw)) + geom_bar(stat="identity") + facet_grid( . ~ Ref) +
   labs(x="Error Type", y = "Count", title = "Errors by Reference") + theme_bw()

cnts
15490/(15490+1325)
48834/(48834+684)
5120/(5120+603)

rd = read.csv("ccs.csv")
head(rd)
hist(rd$NumSubReads,450,col="blue")

ggplot(data=rd, aes(x=NumSubReads)) + 
  geom_histogram(binwidth=5) + 
  facet_grid(.~Reference) + theme_bw() +
  title("Distributions of Subread Counts") + labs(x="Subread count")

#hp = rd[rd$Reference=="ALL4MER.V2.01",]
hp = rd[rd$Reference=="HP.V1.02",]

head(hp)
ggplot(hp, aes(y=NumIndelErrors,x=NumSubReads)) + 
  geom_jitter() + stat_smooth() + labs(x="Sub Read Count",y="Indel Error Count") + title("Errors in HP Reference")
