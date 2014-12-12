setwd('/Users/nigel/CCS_P6_C4/')
library(ggplot2)
library(gridExtra)

d = read.csv("combined_read_report.csv")
aggregate(ZMW~Algorithm,d,length)

dv = read.csv("combined_variant_report.csv")
aggregate(zmw~Algorithm,dv,length)


levels(d$Algorithm) <- c("Polish HP w/ C2", "P6-C4 SumProduct", "Current")


#plot by depth
mkPlot <- function(group) {
  nd = d[d$Reference==group,]
  res = aggregate(NumErrors~NumSubReads+Algorithm+Reference,nd,mean)
  v = ggplot(res,aes(x=NumSubReads,y= NumErrors, color=Algorithm))+geom_line()+theme_classic(base_size=16)
  v = v + labs(x= "Number of Sub Reads", y="Mean Errors Per Read",title=group)
  v
}
a=mkPlot("ALL4MER.V2.01")
mkPlot("HP.V1.02")
b =mkPlot("NOHP.V1.01") 

grid.arrange(a,b, ncol=2)

mkPlot2 <- function(group) {
  nd = d[d$Reference==group,]
  res = aggregate(NumErrors~NumSubReads+Algorithm+Reference,nd,length)
  v = ggplot(res,aes(x=NumSubReads,y= NumErrors, color=Algorithm))+geom_line()+theme_classic(base_size=16)
  v = v + labs(x= "Number of Sub Reads", y="Number of Reads",title=group)
  v
}
mkPlot2("ALL4MER.V2.01")
mkPlot2("HP.V1.02")
mkPlot2("NOHP.V1.01") 


reportDif <- function(group, cutoff=10,minReads=20 ) {
  nd  = d[d$Reference==group & d$NumSubReads > minReads,]
  c4 = nd[nd$Algorithm=="C2",]
  old = nd[nd$Algorithm=="master",]
  cd = merge(old,c4,by="ZMW")
  cd = cd[cd$NumErrors.y<cutoff & cd$NumErrors.x<cutoff  ,]
  v= ggplot(cd, aes(x=NumErrors.x, y = NumErrors.y))+geom_jitter(alpha=.2)+geom_abline(intercept=0,slope=1)+labs(x="Error Count Reads, Current Viterbi", y="Error Count, C2 HP SumProduct Polish",title=group)+theme_bw()
  dif = cd$NumErrors.y-cd$NumErrors.x
  improved = sum(dif<0)/length(dif)
  worsened = sum(dif>0)/length(dif)
  percLess = 1 - sum(cd$NumErrors.y)/sum(cd$NumErrors.x)
  print(group)
  print(paste("Total Original Errors = ", sum(cd$NumErrors.x)))
  print(paste("Percentage Error Reduction = ", percLess))
  print(paste("Percentage Reads Improved  = ", improved))
  print(paste("Percentage Reads Worsened  = ", worsened))
  v
}

reportDif("ALL4MER.V2.01",100,20)
  reportDif("ALL4MER.V2.01",4,50)


reportDif("HP.V1.02",200,1)
reportDif("HP.V1.02",4,50)

reportDif("NOHP.V1.01",200)
reportDif("NOHP.V1.01",4)


sd = dv[dv$Ref=="HP.V1.02" & (dv$homopolymerLength==5 | dv$homopolymerLength==10),]
levels(sd$Algorithm) <- c("Polish HP w/ C2", "P6-C4 SumProduct", "Current")
res = aggregate(Pos~Algorithm+homopolymerLength,sd,length)
cnts = aggregate(zmw~Algorithm,dv,length)
dif = rep(cnts[,2],2)
res$ReadCount = dif
res$ErrorRate = res$Pos / res$ReadCount
res
