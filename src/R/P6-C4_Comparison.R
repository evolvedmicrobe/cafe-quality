### This file is designed to look at 5 bp deletions.  
# It is designed to compare how the parameter settings for the C2 chemistry compares to the P6-C4 Chemistry
library(ggplot2)
library(car)
# First to load the data. 
setwd("/Users/nigel/git/cafe-quality/data/")
#sd = read.csv("homopolymerDeepDive3.csv")
d = read.csv("homopolymerDeepDiveDiagnostics5bpP6go4.csv")
plot(d$HomopolymerLengthFromAlignment, d$HomopolymerLengthFromCounting)

qd = d[d$HomopolymerLengthFromAlignment !="SNP" & d$HomopolymerLengthFromCounting > -10,]
head(qd$HomopolymerLengthFromAlignment)
qd$hpaln = as.numeric(as.character(qd$HomopolymerLengthFromAlignment))
qplot(qd$hpaln, qd$HomopolymerLengthFromCounting)+geom_jitter(alpha=.2)+geom_smooth()
## Now to count how many subreads have null bases, this occurs if any of the following occurs
# 1 - The subread could not align to the reference
# 2 - The subread aligned, but the alignment did not overlap with the homopolymer13/59

bad = d$BaseCalls=="NULL"
percMissing = sum(bad)/nrow(d) # less than 5%
percMissing
cb = function (x) sum(x == "NULL")/ length(x)


# Create a set of only good data
gd = d[!bad,]
gd$ConsensusIndelSize = factor(gd$ConsensusIndelSize)
gd$Zmw = factor(gd$Zmw)
gd$Correct = gd$ConsensusIndelSize==0
gd$ViterbiDifferenceP6 = gd$NoErrorViterbiScoreP6 - gd$OneDeletionErrorViterbiScoreP6
gd$ViterbiDifferenceC2 = gd$NoErrorViterbiScoreC2 - gd$OneDeletionErrorViterbiScoreC2
gd$spdifP6 = gd$NoErrorSumProductScoreP6 - gd$OneDeletionSumProductScoreP6
gd$spdifC2 = gd$NoErrorSumProductScoreC2 - gd$OneDeletionSumProductScoreC2




plot(gd$NoErrorViterbiScoreC2, gd$OneDeletionErrorViterbiScoreC2)

hist(gd$ViterbiDifferenceC2)
bat6 = aggregate(ViterbiDifferenceP6 ~ Zmw + Correct, gd, sum)
bat2 = aggregate(ViterbiDifferenceC2 ~ Zmw + Correct, gd, sum)

b6 = aggregate(spdifP6 ~ Zmw + Correct, gd, sum)
b2 = aggregate(spdifC2 ~ Zmw + Correct, gd, sum)

sum(gd$ViterbiDifferenceP6 > 0) /nrow(gd)
sum(gd$ViterbiDifferenceC2 > 0) /nrow(gd)
sum(gd$spdifP6 > 0) /nrow(gd)
sum(gd$spdifC2 > 0) /nrow(gd)





head(bat2)
vd = merge(b6,b2)
nd = merge(bat2,bat6)
nd=merge(vd,nd)
head(nd)
ggplot(nd,aes(x=ViterbiDifferenceC2, y = ViterbiDifferenceP6, colour=Correct )) + geom_point()+
  geom_abline(slope=1,xintercept=0) + labs(x= "Viterbi Score Difference C2", y = "Viterbi Score Difference P6") + theme_bw(base_size=16) +
  geom_vline(xintercept=0, color="red", lwd=2) 

ggplot(nd,aes(x=spdifC2, y = spdifP6, colour=Correct )) + geom_point() + 
  labs(x= "Sum Product Score Difference C2", y = "Sum Product Score Difference P6") + theme_bw(base_size=16) +
  geom_vline(xintercept=0, color="red", lwd=2) 

sum(nd$ViterbiDifferenceP6 > 0) /nrow(nd)
sum(nd$ViterbiDifferenceC2 > 0) /nrow(nd)
sum(nd$spdifP6 > 0) /nrow(nd)
sum(nd$spdifC2 > 0) /nrow(nd)

cur = (1 - sum(nd$ViterbiDifferenceP6 > 0) /nrow(nd)) 
fut = (1 - sum(nd$spdifC2 > 0) /nrow(nd))
1 - fut/cur
