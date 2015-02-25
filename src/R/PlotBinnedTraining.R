library(ggplot2)
library(scales)
library(gridExtra)
library(mlogit)
library(nnet)
library(reshape2)
setwd('/Users/nigel/git/cafe-quality/NotTracked/')
d= read.csv('Combined.csv')
head(d)
hist(d$Count)

ms = grep(".Merge",colnames(d))
ds = grep(".Dark", colnames(d))
d[,ds] = d[,ms] + d[ds]
d = d[, -ms]

getData <-function(bp) {
  cd = d[d$BP==bp,]
  ctx1 = paste("N", bp, sep="")
  ctx2 = paste(bp, bp, sep="")
  rcol1 = c(1:3, grep(ctx1, colnames(cd)))
  rcol2 = c(1:3, grep(ctx2, colnames(cd)))
  
  ncd = cd[,rcol1]
  cnames = c("BP","SNR","Count","Match","Branch","Dark","Stick")
  colnames(ncd) = cnames
  hcd = cd[,rcol2]
  
  colnames(hcd) = cnames
  nd = rbind(ncd,hcd)
  nd$Context = c(rep(ctx1,nrow(ncd)),rep(ctx2, nrow(hcd)))
  nd$Homopolymer = c(rep(FALSE,nrow(ncd)),rep(TRUE, nrow(hcd)))
  return(nd)
}
a = getData("A")
g = getData("G")
c = getData("C")
t = getData("T")
nd = rbind(a,c,g,t)
nd$HP = factor(nd$HP)
mkplot <- function(y) {
  v = ggplot(nd, aes_string(x="SNR", y = y, color="BP", linetype="Homopolymer"))+geom_point()+geom_line()+theme_bw(base_size=8) + labs(x="Channel SNR")
  return(v)
}
pdf("BinnedSnrEstimates.pdf", width=12, height = 8)
a = mkplot(colnames(nd)[4])
b = mkplot(colnames(nd)[5])
c = mkplot(colnames(nd)[6])
dg = mkplot(colnames(nd)[7])
grid.arrange(a,b,c,dg)
dev.off()


# Now to simulate and fit with a multinomial model
a = getData("C")
moves = c("Match","Branch","Dark","Stick")
cd = a[a$Context=="CC",]
results = c()
contexts = levels(factor(nd$Context))
ctx = "NA"
for(i in 1:nrow(cd)) {
  n=5000
  outcomes = sample(moves,n, prob=cd[i,4:7], replace=TRUE)
  td = data.frame(Context=factor(rep(ctx,n), levels= contexts), Outcome = outcomes, SNR=rep(cd$SNR[i],n))
  results = rbind(results,td)
}
mod<-multinom(Outcome~SNR+I(SNR^2)+I(SNR^3), results)
summary(mod)
exp(coef(mod))

rng = range(results$SNR)
tmpSNR = seq(rng[1], rng[2], .02)
toy = data.frame(SNR=rep(tmpSNR,4),Outcome=rep(moves,each=length(tmpSNR)))
hd = cbind(toy,predict(mod, newdata=toy, type="probs"))
mhd = melt(hd[,-2], id.vars=c("SNR"), value.name="probability")
ggplot(mhd,aes(x=SNR,y=probability,colour=variable))+geom_line()
table(results$SNR,results$Outcome)
#b = mlogit.data(results, shape="wide", choice="Outcome")
#m <- mlogit(Outcome ~ SNR / 0 ,b)
#summary(m)
#head(results)

