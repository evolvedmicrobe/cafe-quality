library("ggplot2")
d = read.csv("/Users/nigel/git/cafe-quality/data/ratios.csv", header=FALSE)
d$BP=factor(d$BP)
colnames(d)<-c("RefSeq","OriginallyRevComped","Zmw","Start","HPLength","BP","Ratio","Insertion","Deletion","Subreads","Alns")
head(d)

plot(d$Subreads,d$Alns)
d$HPLength=factor(d$HPLength)
ggplot(d,aes(x=HPLength,y=Ratio)) + geom_point()+ geom_boxplot()

ggplot(d,aes(x=Insertion,y=Deletion,colour=HPLength)) + geom_point()
plot(d$Ratio,(d$Deletion/(d$Deletion+d$Insertion)))

vd = d[d$HPLength < 6 & d$Alns>50,]
vd$HPLength=factor(vd$HPLength)
ggplot(vd,aes(x=Ratio,fill=HPLength)) + geom_density(alpha = .4)+geom_vline(xintercept=.333)
ggplot(d,aes(x=Subreads,y=Deletion))+geom_smooth()+geom_jitter()

ggplot(vd,aes(x=Alns,y=Ratio,colour=HPLength))+geom_smooth()+geom_jitter()+geom_vline(xintercept=65)

ggplot(vd,aes(x=Subreads,y=Ratio,colour=HPLength))+geom_smooth()+geom_jitter()+geom_vline(xintercept=65)
ggplot(d,aes(x=Subreads,y=Ratio,colour=FailsRatioTest))+geom_smooth()+geom_jitter()+geom_vline(xintercept=65)
ggplot(d,aes(x=Subreads,y=Ratio,colour=RefSeq))+geom_smooth()+geom_jitter()+geom_vline(xintercept=65)


aggregate(FailsRatioTest~RefSeq,d,FUN=sum)
aggregate(FailsRatioTest~RefSeq,d,FUN=length)



hist(d$Start,120)
d$FailsRatioTest = d$Ratio <=.34 & (d$Insertion + d$Deletion)/d$Alns > .1 & (d$Insertion+d$Deletion>=4)
d$FailsRatioTest2 = d$Ratio <=.34 & (d$Insertion + d$Deletion)/d$Alns > .1 

table(d[,grep("FailsRatio",colnames(d))])



ggplot(d,aes(x=Start,y=HPLength,shape=RefSeq,colour=BP))+geom_jitter(size=2)+geom_vline(xintercept=25)+geom_vline(xintercept=40)+labs(x="HP Start on CCS Sequence",y="Jittered CCS HP Length")+theme_bw(base_size=16)
ggplot(d,aes(x=Start,y=HPLength,shape=RefSeq,colour=FailsRatioTest))+geom_jitter()+geom_vline(xintercept=25)+geom_vline(xintercept=45)+labs(x="HP Start on CCS Sequence",y="Jittered CCS HP Length")

hd = d[d$RefSeq=="HP.V1.02" & d$Alns>1 & d$FailsRatioTest == FALSE & d$FailsRatioTest2 == TRUE ,]
hd$BP=factor(hd$BP)
ggplot(hd,aes(x=Start,y=HPLength,colour=FailsRatioTest,shape=BP))+geom_jitter()+geom_vline(xintercept=25)+geom_vline(xintercept=40)+labs(x="HP Start on CCS Sequence",y="Jittered CCS HP Length")+theme_bw(base_size=16)
ggplot(hd,aes(x=Start,y=HPLength,shape=RefSeq,colour=FailsRatioTest))+geom_point()+geom_vline(xintercept=25)+geom_vline(xintercept=40)+labs(x="HP Start on CCS Sequence",y="Jittered CCS HP Length")+theme_bw(base_size=16)

ld = hd[hd$Start>110,]
hist(ld$Ratio,200)
head(ld)
plot(ld$Ratio,ld$FailsRatioTest)
plot(ld$Insertion,ld$Deletion)
abline(0,1)
plot(ld$Ratio,ld$Deletion)
sum(ld$Insertion)/sum(ld$Insertion+ld$Deletion)
sum(ld$Insertion)/sum(ld$Alns)
