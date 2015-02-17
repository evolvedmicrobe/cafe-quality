library(ggplot2)
library(gridExtra)
library(grid)

setwd("/Users/nigel/git/cafe-quality/src/ConstantModelOptimizer/bin/Debug/")
t = read.csv("TrueParameters2.csv")
d = read.csv("Parameters2.csv")

makeEffectiveParameters <- function(d) {
  merges = grep("(A|C|G|T|N)(A|C|G|T)\\.Merge",colnames(d))
  darks = grep("(A|C|G|T|N)(A|C|G|T)\\.Dark",colnames(d))
  d[,darks] = d[,merges]+d[,darks]
  d[,-merges]
}

do = d
to = t 
d = makeEffectiveParameters(d)
t = makeEffectiveParameters(t)
d$Iteration = 1:nrow(d)
head(d)

h = 5
w = 5
pdf("LikelihoodConvergence.pdf", height=h, width = w)
ggplot(d, aes(x = Iteration, y = Likelihood))+geom_point()+geom_line()+theme_bw(base_size=14) + labs(x= "Iteration", y = "Log Likelihood", title = "Likelihood Convergence")
dev.off()

head(t)

maxdif <- function(x) {
  v = abs(x[2:(length(x)-1)] - to[1,])
  print(v)
  v = as.numeric(v)
  print(mean(v[v!=0]))
  mean(v[v!=0])
}
difs = apply(do,1,maxdif)
dd = data.frame(Iteration=1:length(difs), Dif = difs)

maxdif <- function(x) {
  v = abs(x[2:(length(x)-1)] - t[1,])
  v = as.numeric(v)
  mean(v[v!=0])
}
difs2 = apply(d,1,maxdif)
dd2 = data.frame(Iteration=1:length(difs2), Dif = difs2)

dd = rbind(dd,dd2)
dd$ParameterSet = c(rep("Original", nrow(dd2)), rep("Effective", nrow(dd2)))




pdf("MaximumConvergence.pdf",height=h, width = w+2)
ggplot(dd,aes(x=Iteration, y= Dif, colour=ParameterSet))+geom_point()+geom_line()+labs(x="Iteration", y = expression(paste("Mean |", theta, ' - ', hat(theta), '|' )), title="Parameter Convergence")+theme_bw(base_size=14)+scale_y_continuous(limits=c(0,.08))
dev.off()

dif = data.frame(t(d[nrow(d),2:(ncol(d)-1)] - t[1,]))
dif$Rate = rownames(dif)
getType <-function(x) { strsplit(x, "\\.")[[1]][2]}
h=sapply(dif$Rate,getType)
h[1]="Miscall"
dif$Type=factor(h)
dif$Rate = factor(dif$Rate)
colnames(dif)[1] = "Error"
pdf("Differences.pdf", height=6.5, width=10)
ggplot(dif,aes(x=Rate, y = Error, fill=Type))+geom_bar(stat="identity")+ theme_bw(base_size=10)+theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x= "Model Parameter", y = expression(paste( hat(theta), ' - ', theta )), title="Parameter Differences ")
dev.off()
