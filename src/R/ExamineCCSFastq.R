# Data derived from InvestigateFastQ.py 
library(ggplot2)
library(gridExtra)
library(grid)
library(wesanderson)


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
  #d=d[d$Reference!="lambda_NEB3011",]
  d=d[d$Reference=="HP.V1.02",]
  #d = d[d$Reference=="ALL4MER.V2.01",]
  return(d)
}

nd = getld("master_5ba5286_combined_reads.csv")
d = read.csv("qv_compare.csv")

d$Size= factor(d$Size)
d$Type = d$BP:d$Size

aggregate(ZMW~BP+Size+Correct, d, FUN=length)

ld = merge(d, nd, by="ZMW")
lds = ld[ld$NumPasses>20,]

head(lds)
#ggplot(ld, aes(x=Correct, y= QV, color=Size, shape=BP))+geom_point() +geom_jitter()
pdf("QV_Calibration.pdf")
ggplot(lds, aes(x=Correct, y= QV, fill=Correct, color=Correct))+geom_violin() + facet_grid(BP~Size) +
  theme_bw(base_size=12) + labs(title="Violin plots of QV values for different HP calls\nfor reads with > 20 passes") +
  scale_fill_manual(values=wes_palette("Cavalcanti")) +scale_color_manual(values=wes_palette("Cavalcanti"))
dev.off()

hist(d$QV,500)
max(d$QV)

te = d[d$Type=="G:10" & d$QV==42,]
aggregate(ZMW~Correct, te, FUN=length)
aggregate(ZMW~Correct+Type, d, FUN=length)




