setwd("/Users/nigel/git/cafe-quality/doc")
# from http://www.statext.com/practice/MultipleRegressionLinear01.php
d = read.csv("hamilton.csv")

png("Marginal.png",width=10, height=3, units="in", res=120)
par(mfrow=c(1,3))
p1 = lm(y~x1, d)
summary(p1)
plot(predict(p1), d$y, xlab="Predicted", ylab="Observed", main="Model for X1", pch=16) 
abline(coefficients(p1))

p2 = lm(y~x2, d)
summary(p2)
plot(predict(p2), d$y, xlab="Predicted", ylab="Observed", main="Model for X2", pch=16) 
abline(c(0,1))

plot((predict(p1)+predict(p2))/2, d$y,xlab="Predicted", ylab="Observed", main="Combination of Models for X1 and X2", pch=16)
abline(c(0,1))
dev.off()

pred = (predict(p1)+predict(p2))/2
act = d$y
1 - sum((pred-act)^2) / sum((act- mean(act))^2)


# now joint
png("Joint.png",width=3.333, height=4, units="in", res=120)
p = lm(y~., d)
summary(p)
plot(predict(p), d$y,xlab="Predicted", ylab="Observed", main="Model with X1 and X2", pch=16)
abline(c(0,1))
dev.off()




= read.csv("Data.csv")


head(d)
y = lm(Sale~Land.Value, data = d)
summary(y)
y = lm(Sale~Improvements, data = d)
summary(y)
y = lm(Sale~Area, data = d)
summary(y)

y = lm(Sale~., data = d[,-1])
summary(y)


