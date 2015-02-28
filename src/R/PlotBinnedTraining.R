library(ggplot2)
library(scales)
library(gridExtra)
library(mlogit)
library(nnet)
library(reshape2)
setwd('/Users/nigel/git/cafe-quality/NotTracked/')
d = read.csv('Combined.csv')
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
fitModel <- function(ctx) {
  bp = substr(ctx,2,2)
  print(bp)
  a = getData(bp)
  moves = c("Match","Branch","Dark","Stick")
  cd = a[a$Context==ctx,]
  results = c()
  contexts = levels(factor(nd$Context))
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
  pl=ggplot(mhd,aes(x=SNR,y=probability,colour=variable))+geom_line()+theme_bw(base_size=12)+labs(title=ctx)
  tb=table(results$SNR,results$Outcome)
  list(tbl=tb, plot=pl, model=mod, ctx = ctx)
}

ctxs = colnames(d)[grep(".Match",colnames(d))]
ctxs = sapply(ctxs, function(x) {strsplit(x, "\\.")[[1]][1]})
pdf("results.pdf")
for(ctx in ctxs) {
  b = fitModel(ctx)
  plot(b$pl)
  write.csv(coef(b$model), file=paste(b$ctx,".params"))
}
dev.off()
#Verify the observed probabilities make sens
coef(b$model)
snr = seq(6.5,12)
snr = cbind(rep(1,length(snr)),snr, snr^2, snr^3)
co =coef(b$model)
r = co%*%t(snr)
probs = apply(r,2, function(x) exp(x) / (1 + sum(exp(x))))
apply(probs,2, function(x) { 1-sum(x)})
#b = mlogit.data(results, shape="wide", choice="Outcome")
#m <- mlogit(Outcome ~ SNR / 0 ,b)
#summary(m)
#head(results)

