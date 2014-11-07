### This file is designed to look at 5 bp deletions.  
# It analyzes a data file of covariates produced by the FSharp file ExamineHP.fs
# This looks at all reads with 25 bp size of the HP reference template with > 126 subreads
# them to see if 
library(ggplot2)
library(car)
# First to load the data. 
setwd("/Users/nigel/git/cafe-quality/data/")
#sd = read.csv("homopolymerDeepDive3.csv")
d = read.csv("homopolymerDeepDive5bpLong.csv")

head(d)
nrow(d)

## Now to count how many subreads have null bases, this occurs if any of the following occurs
# 1 - The subread could not align to the reference
# 2 - The subread aligned, but the alignment did not overlap with the homopolymer
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
b=summary(gd$HPSizeGroup)
b= data.frame(HPSize=factor(names(b), levels= c("Below -2","-2","-1","0","1","2","Above 2","SNP")) ,Count=b)
ggplot(b,aes(x=HPSize,y=Count))+geom_bar(stat="identity")+theme_classic(base_size=16)

v=aggregate(Zmw~HPSizeGroup+Correct,gd,FUN=length)
ggplot(v,aes(x=HPSizeGroup,y=Zmw,group=Correct,fill=Correct))+geom_bar(stat="identity",position=position_dodge())+theme_classic(base_size=16)+labs(x="Indel Size Relative to Correct Read")

#count number of errors
t2 = aggregate(ConsensusIndelSize~Zmw+Correct,gd,FUN=length)
sum(t2$Correct)/nrow(t2)

# Get how many subreads didn't align by the 
cnts = aggregate(BaseCalls~Zmw+ConsensusIndelSize, data = gd, FUN=cb)
head(cnts)

#Verify our result matches the expected outcome
res=aggregate(ViterbiDifference~Zmw+Correct+ConsensusIndelSize,data=gd,FUN=sum)
ggplot(res, aes(x = ViterbiDifference))+geom_histogram( aes(y=..density..)) + theme_bw(base_size=16)+facet_grid(.~Correct)+labs(title="Histograms of Summed Difference in Subread Viterbi Scores by CCS being Correct")
# How many am I off by
rtemp = res[res$Correct==FALSE,]
summary(rtemp$ConsensusIndelSize)
sum(rtemp$ViterbiDifference>0) / nrow(rtemp)


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
ggplot(resneg, aes(x = HomopolymerLength, fill=Correct))+geom_density(alpha=.8) + theme_bw(base_size=16)+labs(x = "Percentage of Reads With Deletion", title="Density of % Subreads with HP Deletion by whether CCS is Correct")
tmp = merge(resneg,resExist,by=c("Zmw","Correct"))
colnames(tmp)[3:4] <- c("PD","SubReadCount")
ggplot(tmp,aes(x=SubReadCount,y=PD,colour=Correct))+geom_jitter()+labs(y="Percentage of Reads With Deletion",x="Usable Subread Count")+theme_bw(base_size=16)+geom_smooth()

q=aggregate(Zmw~SubReadCount,tmp,FUN=length)
sum(q$Zmw)

#how many would be fixed by simple voting?
wrong = tmp$Correct==FALSE
sum(tmp[wrong,"PD"]<.5) / length(wrong)


#most frequent seems to be 120
n=120
tmp2 = tmp[tmp$SubReadCount==n,]
mle = mean(tmp2$PD)
stdev = sqrt(n*mle*(1-mle))
getP <-function() { sum(runif(n)<mle)}#/120}
res=replicate(12000,getP())
toP =data.frame(group=c(rep("Real",nrow(tmp2)),rep("Simulated",length(res))), freq = c(tmp2$PD,res/120))
ggplot(toP, aes(freq, ..density.., colour = group,fill=group)) +geom_density(alpha=.2) +   geom_freqpoly(binwidth=1/120)+theme_bw(base_size=20)+labs(x="Percent Errors in SubReads")


#make qqplot of comparisons
res2 = rnorm(12000,120*mle,sqrt(120*mle*(1-mle)))
qqnorm((n*tmp2$PD-mle*n)/stdev)
abline(0,1)
ggplot()
qqnorm((res2-mle*n)/stdev)
abline(0,1)
hist(tmp2$PD,50)

# how does the plus one minus one rule stack up?
rule = rep(-999,nrow(gd))
sub = gd$Correct==TRUE
vals = gd[!sub,]
cntinc <- function(x) {
  r1 = sum(x == "-2")
  r2 = sum(x == "0")
  r1/(r2+r1)  
}
cntcor <- function(x) {
  r1 = sum(x == "-1")
  r2 = sum(x == "1")
  r1/(r2+r1)  
}
set1 = aggregate(HPSizeGroup~Zmw,data=vals,FUN=cntinc)
set1$Correct = rep(FALSE,nrow(set1))
vals2 = gd[sub,]
set2 = aggregate(HPSizeGroup~Zmw,data=vals2,FUN=cntcor)
set2$Correct = rep(TRUE,nrow(set2))
tmp = rbind(set1,set2)
ggplot(tmp,aes(x=HPSizeGroup,colour=Correct))+geom_histogram(aes(y=..density..))+facet_grid(.~Correct)+labs(x="Nigel's `Your Wrong` Score")
ggplot(tmp,aes(x=HPSizeGroup,fill=Correct))+geom_density(aes(y=..density..))+labs(x="Read Bias Score", title = "Read Bias by Length")+theme_bw(base_size=12)

#Count Fixes
sum(set1$HPSizeGroup > .5)/nrow(set1)

# END RULE CHECK




cnts = aggregate(BaseCalls~Zmw+ConsensusIndelSize, data = gds, FUN=cb)
head(cnts)


#Create a sampled set of good data
smp = .2
gds = gd[sample.int(nrow(gd),nrow(gd)*smp),]


ggplot(gds, aes(x = ViterbiDifference))+geom_histogram( aes(y=..density..)) + theme_bw(base_size=16)+facet_grid(.~HPSizeGroup)+labs(title="Histograms of Difference in Subread Viterbi Score by Subread HP Length")
ggplot(gds, aes(x = NoErrorViterbiScore, y = OneDeletionErrorViterbiScore, colour=HPSizeGroup))+geom_point(alpha=.85) + theme_bw(base_size=16)+geom_abline(intercept=0,slope=1)+labs(x="Correct Sequence Score",y="Incorrect Sequence Score",title = "Likelihood of Templates For Subreads")
ggplot(gds, aes(x = NoErrorViterbiScore, y = OneDeletionErrorViterbiScore, colour=Correct))+geom_point(alpha=.85) + theme_bw(base_size=16)+geom_abline(intercept=0,slope=1)+labs(x="Correct Sequence Score",y="Incorrect Sequence Score",title = "Likelihood of Templates For Subreads")


ggplot(gds, aes(x = NoErrorViterbiScore, y = OneDeletionErrorViterbiScore, colour=HomopolymerLength))+geom_point(alpha=.85) + theme_bw(base_size=16)+geom_abline(intercept=0,slope=1)+labs(x="Correct Sequence Score",y="Incorrect Sequence Score",title = "Likelihood of Templates For Subreads")




# Verify hplength variable
len <- function(x) {
  v = as.character(x)
  nchar(v)
}
l = sapply(gds$BaseCalls,len)
plot(gds$HPSectionLength,l)
hist(gds$HPSectionLength,50)

# Count up most frequent Haplotypes
cnts = aggregate( Zmw ~ BaseCalls, data = gds, FUN = length)
cnts = cnts[order(-cnts$Zmw),]
cnts$Zmw = cnts$Zmw/nrow(gds)
colnames(cnts)[2]<- "PerecentageSubReads"
cnts[1:10,]

rd = gds[gds$BaseCalls%in%cnts$BaseCalls[1:10],]
rd$BaseCalls= factor(rd$BaseCalls)
ggplot(rd,aes(x=NoErrorViterbiScore ,y = OneDeletionErrorViterbiScore, colour=BaseCalls,shape=HomopolymerLength))+geom_point()+theme_bw(base_size=16)+geom_abline(intercept=0,slope=1)


 ### RUN T-SNE
library(tsne)
tosne = rd[,!(colnames(rd)%in%c("BaseCalls", "Zmw","Correct"))]
tosne$HomopolymerLength=as.numeric(as.character(tosne$HomopolymerLength))
tosne$ConsensusIndelSize=as.numeric(as.character(tosne$ConsensusIndelSize))
tosne$ReverseComplementedOriginally=as.numeric(tosne$ReverseComplementedOriginally)
summary(tosne)
res=tsne(tosne)



plot(1:nrow(cnts),cumsum(cnts$PerecentageSubReads),xlab= "SubRead Sequence (Ordered by Frequency)",ylab="Cumulative Frequency")



toR = gds[gds$NoErrorViterbiScore > -40 & gds$HPSectionLength < 25,!(colnames(gds)%in%c("BaseCalls","Zmw","ViterbiDifference","OneDeletionErrorViterbiScore","NoErrorSumProductScore","OneDeletionSumProductScore"))]
y=lm(NoErrorViterbiScore~.,data=toR)
summary(y)
b=step(y)
summary(b)
y=lm(NoErrorViterbiScore~HPSectionLength+I(HPSectionLength^2)+HomopolymerLength+ConsensusIndelSize,data=toR)
summary(y)
crPlots(y)
lenRange = seq(0,82,1)
ylen = coefficients(y)[1]+coefficients(y)[2]*lenRange
mn = data.frame(HPSec=lenRange,Score= ylen)
p = ggplot(toR,aes(x=HPSectionLength, y = NoErrorViterbiScore)) + geom_point()
p = p + geom_line(data=mn,aes(x=HPSec,y=Score))   
p
lenRange,ylen,pch=18, xlab= "HP Section Length", ylab = "Mean Predicted Score")

y=lm(NoErrorViterbiScore~HPSectionLength+Mean_InsertionQV+HPSectionLength:Mean_InsertionQV+I(Mean_InsertionQV^2),data=toR)
summary(y)
y=lm(NoErrorViterbiScore~HPSectionLength+Mean_InsertionQV+Max_InsertionQV,data=toR)
summary(y)
plotplot(y)

ggplot(gds, aes(x = ViterbiDifference, group=Correct))+geom_histogram(aes(y=..density..)) + theme_bw()+facet_grid(.~ReverseComplementedOriginallyTrue)

ggplot(gds, aes(x = ViterbiDifference, group=Correct))+geom_histogram(aes(y=..density..)) + theme_bw()+facet_grid(.~ReverseComplementedOriginallyTrue)




gds$SubRight = gds$ViterbiDifference > 0
y = glm(SubRight~ReverseComplementedOriginally+, data=gds, family="binomial")
summary(y)

bad = c("BaseCalls","Zmw","ViterbiDifference","NoErrorViterbiScore","OneDeletionErrorViterbiScore","NoErrorSumProductScore","OneDeletionSumProductScore")
nd = gds[,!(colnames(gds)%in%bad)]
red = glm(SubRight~1,data=nd, family="binomial")
full = glm(SubRight~., data=nd, family="binomial")
summary(full)

y=step(red,scope=list(lower=red,upper=full))


y = glm(SubRight ~ HPSectionLength + Mean_DeletionQV + ReverseComplementedOriginally +  Mean_MergeQV + Mean_InsertionQV + Min_MergeQV + Mean_SubstitutionQV +  Max_InsertionQV + Min_InsertionQV + Max_DeletionTag + Min_IpdInFrames + SpikeMergeQVCount + Min_DeletionQV + Min_DeletionTag + RQ + SubReadNumber + Max_PulseWidthInFrames + HomopolymerLength, data=nd,family="binomial")
summary(y)
y = glm(SubRight ~ HPSectionLength, data=nd,family="binomial")
summary(y)

hist(gds$NoErrorViterbiScore - gds$OneDeletionErrorViterbiScore,500)


plot(gds$ViterbiDifference, gds$HPSectionLength)

aggregate(HPSectionLength~SubRight, gds, median)

hist(gdnrow(gds)

     
# Model based accuracy predictions
getAccuracy <- function(cov) {
  reps = 1000
  getSamp <- function() { sum(gds$ViterbiDifference[sample.int(nrow(gds),cov)]) }
  res = replicate(reps,getSamp())
  return(sum(res>0)/reps)
}
cov = seq(1,150,5)
y=sapply(cov, FUN=getAccuracy)
plot(cov,y)
