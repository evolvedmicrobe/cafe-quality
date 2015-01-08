""" THIS SCRIPT LOADS DATA PRODUCED BY THE F# SCRIPT ExamineDelTags.fs."
setwd("/Users/nigel/git/cafe-quality/data/")
d = read.csv("DelTagInvestigation.csv")
#Reduce data down to useful things
d = d[d$AssignedReferenceName!="NULL",colnames(d)%in%c("AssignedReferenceName",
                        "ReverseComplementedOriginally",
                        "BaseCalls",
                        "RQ",
                        "OriginalSubReadLength",
                        "Zmw",
                        "AlignedLength",
                        "CountDelTags",
                        "CountCorrectDelTags")]
head(d)

d$AssignedReferenceName = factor(d$AssignedReferenceName)
levels(d$AssignedReferenceName)



hist(d$CountDelTags,50)

#####################################
d$RQF = factor(d$RQ)
d$Rate = d$CountDelTags / d$AlignedLength
rateByValue = aggregate(Rate~RQF,d,mean)
rateByValue$RQ = as.numeric(as.character(rateByValue$RQF))
pdf("TagRateByRQ.pdf",width=3.5, height=3.5)
ggplot(rateByValue, aes(x=RQ, y=Rate))+geom_point()+geom_smooth() + labs(x="Read Quality", y = "Percentage of Bases with a DelTag")+theme_bw(base_size=10)+scale_y_continuous(labels=percent)
dev.off()

#####################################
d$CorrectRate = d$CountCorrectDelTags / d$CountDelTags
hist(d$CorrectRate,50)
plot(d$CountDelTags,d$CorrectRate, col=d$AssignedReferenceName)

x=aggregate(CountCorrectDelTags~RQF+AssignedReferenceName,d,sum)
y=aggregate(CountDelTags~RQF+AssignedReferenceName,d,sum)
z= merge(x,y)
z$RateCorrect = z$CountCorrectDelTags / z$CountDelTags
z$RQ = as.numeric(as.character(z$RQF))
pdf("TagCorrectRateByRQ.pdf",width=3.5, height=2.5)
ggplot(z, aes(x=RQ, y=RateCorrect))+geom_point()+geom_smooth() + labs(x="Read Quality", y = "Percentage of 'Correct' DelTags ")+theme_bw(base_size=10)+scale_y_continuous(labels=percent)
dev.off()

