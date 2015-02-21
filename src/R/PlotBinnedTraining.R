library(ggplot2)
library(scales)

setwd('/Users/nigel/git/cafe-quality/NotTracked/')
d= read.csv('Combined.csv')
head(d)
hist(d$Count)

ms = grep(".Merge",colnames(d))
ds = grep(".Dark", colnames(d))ds

d[,ds] = d[,ms] + d[ds]
d = d[, -ms]

gd = d[d$BP=="G",]
head(dg)
ggplot(gd, aes(x=SNR, y = GG.Stick))+geom_point()+geom_line()

bp = "A"
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
  nd$HP = c(rep(FALSE,nrow(ncd)),rep(TRUE, nrow(hcd)))
  
  return(nd)
}
a = getData("A")
g = getData("G")
c = getData("C")
t = getData("T")
nd = rbind(a,c,g,t)
nd$HP = factor(nd$HP)
mkplot <- function(y) {
  v = ggplot(nd, aes_string(x="SNR", y = y, color="Context", linetype="HP"))+geom_point()+geom_line()+theme_bw(base_size=20)
  return(v)
}
mkplot(colnames(nd)[4])
mkplot(colnames(nd)[5])


do.call(grid.arrange,b)
