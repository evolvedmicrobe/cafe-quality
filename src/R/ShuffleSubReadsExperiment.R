### This file is designed to look at 5 bp deletions.  
# It analyzes a data file of covariates produced by the FSharp file ExamineHP.fs
# This looks at all reads with 25 bp size of the HP reference template with > 126 subreads
# them to see if 
library(ggplot2)
library(car)
library(RColorBrewer)
# First to load the data. 
setwd("/Users/nigel/git/cafe-quality/data/")
d = read.csv("homopolymerDeepDive5bpLong.csv")
head(d)
nrow(d)


#Merge in SNR Data?
snrD = read.csv("/Users/nigel/git/cafe-quality/NotTracked/master_full/master_5ba5286_combined_reads.csv")
snrD = snrD[snrD$Reference=="HP.V1.02",]
head(snrD)
colnames(snrD)[colnames(snrD)=="ZMW"]<-"Zmw"


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
gd$ViterbiDifference = gd$NoErrorViterbiScore - gd$OneDeletionErrorViterbiScore
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
res = rep(NA,nrow(gd))
for (i in 1:nrow(gd)) {
  if(gd$ReverseComplementedOriginally[i]=="True") {c = 'C'} else {c = 'G'}
  s = gd$BaseCalls[i]
  s = as.character(s)
  s2 <- gsub(c,"",as.character(s))
  res[i]= (nchar(s) - nchar(s2))
}
gd$NumC = res
aggregate(NumC~Correct, gd, mean)
# FINISH DATA LOAD



percneg <- function(x) {
  x = as.numeric(as.character(x[x!="SNP"]))
  sum(x < 0) / length(x)
}
exists<- function(x) {
  x = as.numeric(as.character(x[x!="SNP"]))
  length(x)
}
resneg=aggregate(HomopolymerLength~Zmw+Correct,data=gd,FUN=percneg)
resExist = aggregate(HomopolymerLength~Zmw+Correct,data=gd,FUN=exists)
ggplot(resneg, aes(x = HomopolymerLength, fill=Correct))+geom_density(alpha=.8) + theme_bw(base_size=16)+labs(x = "Percentage of Subreads With Deletion", title="Density of % Subreads with HP Deletion by whether CCS is Correct")

head(resExist)

tmp = merge(resneg,resExist,by=c("Zmw","Correct"))
colnames(tmp)[3:4] <- c("PD","SubReadCount")
ggplot(tmp,aes(x=SubReadCount,y=PD,colour=Correct))+geom_jitter()+labs(y="Percentage of Subreads With Deletion",x="Usable Subread Count")+theme_bw(base_size=16)+geom_smooth()


hist(tmp$SubReadCount)


n=110
tmp2 = tmp[tmp$SubReadCount>=n,]


tmp2$Correct=factor(tmp2$Correct)

y = glm(Correct~PD+I(PD^2),data=tmp2,family="binomial")
summary(y)
#Get decision point
co = coefficients(y)
crit_point = -co[1]/co[2]
a=co[3]
b=co[2]
c=co[1]
crit_point = (-b - sqrt(b^2 - 4*a*c)) / (2*a)


pdf("DecisionPoint.pdf",width=10, height=8)
ggplot(tmp2, aes(x = PD, fill=Correct))+
  geom_density(alpha=.8) + theme_bw(base_size=16)+
  labs(x = "Percentage of subreads with a deletion", title="Density of % subreads from a 5 bp HP with a deletion\n divided by whether CCS is correct", y="Density") +
  geom_vline(xintercept=crit_point,color="red", size =2, linetype="dashed", name="Decision Point")+scale_fill_discrete(name="CCS correctly called")
dev.off()

#Verify that the critical function is at a minimum
numWrong <-function(x) {
  del = tmp2$PD < x
  right_del = sum(del & (tmp2$Correct=="TRUE"))
  right_bad = sum((!del) & (tmp2$Correct=="FALSE"))
  return(1 - (right_del+right_bad)/nrow(tmp2))
}
cp = seq(.3,.5,.01)
yy=sapply(cp,numWrong)
plot(cp,yy)



#Get simulated accuracy
getRes <-  function() {
  #First sample a read
  n = sample(tmp2$SubReadCount,1)
  freqs = sample(tmp2$PD, n)
  reads = rbinom(n,1,freqs)
  res = sum(reads)/n
  #return whether called correctly
  return (res<crit_point)
}
results = replicate(2000,getRes())
percMixed = sum(results)/length(results) # percent correct when mixing
percActual = sum(tmp2$Correct) / nrow(tmp2) # percent correct in CCS
percRule= sum(tmp2$PD<crit_point)/ nrow(tmp2) # percent correct by simple rule


exp = data.frame(Group=factor(c("CCS", "CCS Reads w/\nRule", "CCS Reads w/\nRule + Mixing")), pc = c(percActual, percRule, percMixed))
#exp$Group=factor(exp$Group,levels=levels(exp$Group)[c(1,3,2)])
exp$Error = 1 - exp$pc
exp$QV = -10*log10(exp$Error)
exp
pdf("MixingExperiment.pdf",width=)
ggplot(exp, aes(x=Group,y=QV, fill=Group) )+geom_bar(stat="identity") + 
  theme_classic(base_size=16) + labs(y="Empirical Error Rate for 5 bp G HP (Phred Scaled)", x="") +
  scale_fill_manual(guide=FALSE, values=brewer.pal(3,"Set1")) + 
  geom_text(aes(x=Group, y=QV, ymax=QV, label=sprintf("%1.2f", QV)), position = position_dodge(width=1), vjust= -1)
dev.off()


#most frequent seems to be 120
n=120
tmp2 = tmp[tmp$SubReadCount==n,]
mle = mean(tmp2$PD)
stdev = sqrt(n*mle*(1-mle))
getP <-function() { sum(runif(n)<mle)}#/120}
res=replicate(12000,getP())
toP =data.frame(group=c(rep("Real",nrow(tmp2)),rep("Simulated",length(res))), freq = c(tmp2$PD,res/120))
pdf("DecisionPointWhenSimulated.pdf",width=10, height=8)
ggplot(toP, aes(freq, ..density..,fill=group)) +geom_density(alpha=.6) +
  theme_bw(base_size=20)+labs(x="Percentage of subreads with a deletion", y="Density", title="Comparison of real and simulated 120X data")+ 
  geom_vline(xintercept=crit_point,color="red", size =2, linetype="dashed", name="Decision Point")+scale_fill_manual(name="Data", values=brewer.pal(3,"Set1")[1:2])
dev.off()
#geom_freqpoly(binwidth=1/120)+



tmp2 = tmp[tmp$SubReadCount>110,]
combined = merge(tmp2,snrD)
combined$MeanSNR = .5*(combined$SnrG+combined$SnrC)
head(combined)
pdf("DeletionRateBySNR.pdf",width=10, height=8)
ggplot(combined,aes(x=MeanSNR,y=PD)) + geom_point() + geom_smooth() + labs(y="Percentage of subreads with a deletion", x=" Mean G + C SNR",
                                                                           title="Deletion Rate by SNR") +theme_bw(base_size=16) +
  geom_hline(yintercept=crit_point,color="red", size =2, linetype="dashed", name="Decision Point")      
dev.off()











