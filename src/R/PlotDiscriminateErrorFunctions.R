


x = seq(-13.8, 4.5, .01)
y = -log(1+exp(x))
plot(x,y)


y2 = 1 / (1 + exp(-x)) 
plot(x, y2)
n = length(x)
y3 = rep(1,n)
y3[x<0] = 0
xs = c(x,x,x)
ys = c(-y,y2,y3)
group = c(rep("Log-Sigmoid (CCS)", n), rep("Logistic (Quiver/MCE)", n), rep("Empirical Error Rate",n))
d = data.frame(xs,ys,group)

ggplot(d, aes(x=xs, y = ys, colour=group))+geom_line(lwd=1.2)+ theme_bw(base_size=16)+labs(x="Mutations Score Difference", y = "Function Value", title="Scoring Schemes")+scale_colour_discrete(name="Scoring Function\n(per mutation)")
