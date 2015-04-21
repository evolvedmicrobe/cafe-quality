library(ggplot2)
library(scales)
library(gridExtra)
library(mlogit)
library(nnet)
library(reshape2)
setwd('/Users/nigel/git/cafe-quality/NotTracked/')
d = read.csv("qv_compare.csv")
d2 = read.csv("qv_compare_master.csv")
q = read.csv("master_5ba5286_combined_reads.csv")
b = q$ZMW[q$NumPasses > 20]
ins = intersect(intersect(d$ZMW, d2$ZMW), b)
ins = b
cd = rbind(d,d2)
cd$Analysis=factor(c(rep("Phase1", nrow(d)), rep("Original", nrow(d2))))
cd = cd[cd$ZMW%in%ins,]


ggplot(cd, aes(x=QV, fill=Analysis))+geom_density()+facet_wrap(Correct~.)

b = aggregate(ZMW~Correct+Analysis+Size+BP, data = cd[cd$QV ==42,], FUN=length)
b
