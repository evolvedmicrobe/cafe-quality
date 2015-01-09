### This file is designed to look at relatiionships between QV values 
# It analyzes a data file of covariates produced by the FSharp file ExamineQVValues.fs
# This looks at a sample of bases from reads with 25 bp size of the HP reference template 
library(ggplot2)
library(gridExtra)
library(reshape2)
# First to load the data. 
setwd("/Users/nigel/git/cafe-quality/data/")
d = read.csv("QVvalues.csv")
head(d)


nrow(d)

dns = sample(1:nrow(d),1e4)
nd = d[dns,]
head(nd)
nd$Zmw <- NULL
nd$MergeQV.1 <- NULL
plot(nd)

cols = c("MergeQV", "DeletionQV", "DeletionTag", "SubsQv", "InsQV")
pg <- function(name) {
  p=ggplot(nd[nd$DeletionTag!=78,],aes_string(x=name))+geom_histogram(binwidth=1)
  p=p+geom_density()+theme_classic()+labs(y="Count")
  p
}
g=lapply(cols, pg)
do.call(grid.arrange,g)


pgb <- function(name) {
p=ggplot(nd[nd$DeletionTag!=78,],aes_string(x=name))+geom_histogram(aes(y=..density..),colour="black", fill="white",binwidth=1)
p=p+geom_density()+theme_classic() +facet_grid(.~Base)
p
}
g=lapply(cols[2:length(cols)], pgb)
do.call(grid.arrange,g)

sum(nd$DeletionTag==78)/nrow(nd)

# Now to look at recalibration
useChemP6=TRUE

if(useChemP6) {
  
  ### Retrained P6-C4 parameters
  v_lab = "- P6-C4 Chemistry"
  v_name = "qv_violin_p6c4.pdf"
  matchVal = ### Retrained P6-C4 parameters
    v_lab = "- P6-C4 Chemistry"
  v_name = "qv_violin_p6c4.pdf"
  matchVal = 0.139078334
    
  delNVal = -1.45442641 
  recal = data.frame( Mismatch =c (-1.85595274, -0.0312577896), 
                      Branch = c(-0.426779211, -0.0509828515),
                      DeletionWithTag  = c(0.0239052698, -0.0412806273),
                      Nce = c(-0.165034905,-0.0826156288),
                      Merge_A =c ( -0.656411886, -0.0198618043),
                      Merge_C = c( 0.0718991458, -0.0777479857),
                      Merge_G = c(0.281680286, -0.0718968958),
                      Merge_T = c(-0.717717707, -0.149290472))
  
  
  #Original ones below 
  if(TRUE) {
    #p6_c4 parameters
    v_lab = "- P6-C4 Chemistry"
    v_name = "qv_violin_p6c4.pdf"
    matchVal = .494617
    delNVal = -3.07623 
    recal = data.frame( Mismatch =c (-3.45427,-0.0199875), 
                        Branch = c(-0.28307, -0.121895),
                       DeletionWithTag  = c(-0.0134004, -0.0844901),
                       Nce = c(-0.0213797,-0.175615),
                       Merge_A =c ( -1.95056, -0.0357092),
                      Merge_C = c(-0.0101704, -0.0878416),
                      Merge_G = c(-0.225816, -0.0759892),
                      Merge_T = c(-0.344726,-0.0570072))
  
  }
} else {
  

## c2 Chemistry
  v_lab = "- C2 Chemistry"
  v_name = "qv_violin_c2.pdf"
  delNVal   =  -2.740747389
  matchVal  =   0.560464158
  recal = data.frame( Mismatch    = c(-3.24367743, -0.066470768),
                         Branch	=   c(-1.215244837, -0.101115918),                       
                         DeletionWithTag	=c ( 0.049440612, -0.093502954),
                         Nce =	c ( -0.251601059, -0.199185423),
                         Merge_A	= c( -1.510719452, -0.035472633),
                         Merge_C	=  c( 0.133876836,   -0.161696672),
                         Merge_G	=  c(0.513768983, -0.181189888),
                         Merge_T	=  c( -1.316929942, -0.287860148))
}

rd = nd
head(rd)
rd$DeletionQV = nd$DeletionQV*recal$DeletionWithTag[2] + recal$DeletionWithTag[1]
d=qplot(nd$DeletionQV,rd$DeletionQV,xlab="Original Score", ylab = "Recalibrated Score")+geom_jitter()+theme_classic()+labs(title="DeletionQV")+geom_abline(slope = -.1*log(10),colour="red")
d
rd$SubsQv = nd$SubsQv*recal$Mismatch[2] + recal$Mismatch[1]
s=qplot(nd$SubsQv,rd$SubsQv,xlab="Original Score", ylab = "Recalibrated Score")+geom_jitter()+theme_classic()+labs(title="SubQV")+geom_abline(slope = -.1*log(10),colour="red")
s
rd$InsQvBranch = nd$InsQV*recal$Branch[2] + recal$Branch[1]
i1=qplot(nd$InsQV,rd$InsQvBranch,xlab="Original Score", ylab = "Recalibrated Score")+geom_jitter()+theme_classic()+labs(title="InsQV - Branch Recalibration")+geom_abline(slope = -.1*log(10),colour="red")

rd$InsQvNce = nd$InsQV*recal$Nce[2] + recal$Nce[1]
i2=qplot(nd$InsQV,rd$InsQvNce,xlab="Original Score", ylab = "Recalibrated Score")+geom_jitter()+theme_classic()+labs(title="InsQV - Nce Recalibration")+geom_abline(slope = -.1*log(10),colour="red")

grid.arrange(d,s,i1,i2)

qplot(rd$InsQvBranch,rd$InsQvNce)+geom_abline(intercept=0,slope=1)+labs(x="Branch Recalibrated Insertion", y = "NceRecalibrated Insertion",title= "InsQV Comparison")+theme_classic(base_size=14)


#now for merge scores
bases = levels(rd$Base)
for(bp in bases)
{
  ro = rd$Base==bp
  cos = recal[,paste("Merge_",bp,sep="")]
  rd$MergeQV[ro] = cos[2]*rd$MergeQV[ro]+cos[1]
}
head(rd)
name = "MergeQV"
p=ggplot(nd,aes_string(x=name))+geom_histogram(aes(y=..density..),colour="black", fill="white",binwidth=1)
p=p+geom_density()+theme_classic() +facet_grid(.~Base) + labs(title="Original")

q=ggplot(rd,aes_string(x=name))+geom_histogram(aes(y=..density..),colour="black", fill="white",binwidth=0.1)
q=q+geom_density()+theme_classic() +facet_grid(.~Base) + labs ( title = "Recalibrated")

grid.arrange(p,q)
nd$RecalMergeQV = rd$MergeQV
v=ggplot(nd, aes(x=MergeQV, y=RecalMergeQV, group=Base, colour = Base ))+geom_point(size = 2)+theme_classic(base_size=16)+labs(x = "Original Merge QV", y= "Recalibrated Merge QV")
v=v+geom_abline(slope = -.1*log(10),colour="red")
v

# Now for violin plot.v
vd = rd
vd$DeletionQV[vd$DeletionTag==78]=NA
#now to make A,C,G,T distributions
bases = levels(rd$Base)
for(bp in bases) {
  nm =paste("Merge_",bp,sep="")
  vd[,nm] = rep(NA,nrow(vd))
  ro = vd$Base==bp
  vd[ro,nm]=rd$MergeQV[ro]
}
vd2 = vd[,c(3,5,7,8,9,10,11,12)]
vd2$row = 1:nrow(vd2)
vd3 = melt(vd2,id="row")
n=1000

#fakeData=data.frame(row=(nrow(vd3):(nrow(vd3)+n-1)),variable=rep("Match",n),value=rep(matchVal,n)+rnorm(0,0.01))
#vd4=rbind(vd3,fakeData)
#head(vd4)
#trimm the tons of 27 values
tooMany = (1:nrow(vd3))[vd3$variable=="SubsQv" & vd3$value==median(vd$SubsQv)]
rm = sample(tooMany,length(tooMany)*.95)
vd4=vd3[-rm,]
low_y = min(min(vd4$value,na.rm=TRUE),delNVal)
high_y = max(max(vd4$value,na.rm=TRUE),matchVal)
pdf(paste("../doc/",v_name,sep=""),width=6,height=4.45)
ggplot(vd4,aes(variable,value,fill=variable))+geom_violin(show_guide=FALSE)+theme_classic(base_size=16)+
  labs(y="Recalibrated Score",x="",title=paste("Distribution of Quality Scores at Bases",v_lab))+
  geom_hline(yintercept=matchVal,colour="red")+ scale_y_continuous(limits=c(low_y, high_y)) +
  geom_hline(yintercept=delNVal,colour="blue") +
  geom_hline(yintercept=delNVal+matchVal,colour="blue",linetype="dotted")
dev.off()

b= aggregate(value~variable,vd4,median)

#Calculate case 1, no tag deletion
matchVal*3+delNVal
matchVal*3+b[b[,1]=="InsQvBranch",2]

#Case case 2, tag deletion
matchVal*3+b[b[,1]=="DeletionQV",2]

# Case 3, merge deletion
matchVal*2+b[b[,1]=="Merge_G",2]
matchVal*2+b[b[,1]=="Merge_A",2]


# No Tag path
matchVal*3 + (matchVal*3+b[b[,1]=="InsQvBranch",2]) #3 bp average score for two reads
(matchVal*3+delNVal) + matchVal*4  # 4bp average

# Tag score path
matchVal*3 + (matchVal*3+b[b[,1]=="InsQvBranch",2]) # 3 bp score
(matchVal*3+b[b[,1]=="DeletionQV",2]) + matchVal*4  # 4 bp score 

# Merge score path
matchVal*3 + (matchVal*3+b[b[,1]=="InsQvBranch",2]) # 3 bp score
(matchVal*2+b[b[,1]=="Merge_G",2]) + matchVal*4  # 4 bp score 



# intercept / slope
