library(ggplot2)
library(car)
# First to load the data. 
setwd("/Users/nigel/git/cafe-quality/data/")

d = read.csv("homopolymerDeepDiveDiagnostics.csv")
bad = d$BaseCalls=="NULL"
percMissing = sum(bad)/nrow(d) # less than 5%
percMissing
cb = function (x) sum(x == "NULL")/ length(x)
# Create a set of only good data
gd = d[!bad,]


head(d)
gd$FullDif=gd$FullLengthCorrectSumProductScore - gd$FullLengthIncorrectSumProductScore
bd=gd[!gd$Correct,]
sum(bd$FullDif>0)/nrow(bd)
gd$Dif = gd$NoErrorSumProductScore - gd$OneDeletionSumProductScore

gd$Correct = gd$ConsensusIndelSize=="0"
cate <- function(x) {
  if(x=="SNP") {return("SNP")}
  x = as.numeric(as.character(x))
  res ="none"
  if(x < -2) {res="Below -2"}
  for(i in seq(-2,2,1)) {
    if(x==i)
    {
      res = as.character(i)
      break
    }
  }
  if(x > 2) {
    res = "Above 2"
  }
  return(res)
}
gd$HPSizeGroup = factor(sapply(gd$HomopolymerLength,cate), levels= c("Below -2","-2","-1","0","1","2","Above 2","SNP"))
smp = .025
gds = gd[sample.int(nrow(gd),nrow(gd)*smp),]


b = aggregate(FullDif~Zmw+Correct,gd,sum)
bad = b[b$Correct==FALSE,]
sum(bad$FullDif>0)/nrow(bad)
pdf("ScoreDifSumProduct.pdf")
ggplot(b,aes(x=FullDif,fill=Correct))+geom_density(alpha=.5)+geom_vline(xintercept=0,color="red")+theme_bw(base_size=16)+labs(x="Score Difference Between Correct and Incorrect Sequence", y = "Density")
dev.off()

ggplot(d,aes(x=Dif, y=SummedPulseWidthForHP, colour=HPSizeGroup))+geom_point()

ggplot(gds,aes(x=Dif,y=FullDif,colour=Correct))+geom_point()
