## Program to examine data in file
#Change to directory with data
setwd("/Users/nigel/git/cafe-quality/doc")
# Load data
d = read.csv("hamilton.csv")

# Fit a model just on the x1 variable
p1 = lm(y~x1, d)
# Assess model fit
summary(p1)
# Plot that model
plot(predict(p1), d$y, xlab="Predicted", ylab="Observed", main="Model for X1", pch=16) 
abline(c(0,1)) # this plots the 1:1 line for obs/predicted

# Now fit a model just on the x2 variable
p2 = lm(y~x2, d)
summary(p2)
# Plot that model
plot(predict(p2), d$y, xlab="Predicted", ylab="Observed", main="Model for X2", pch=16) 
abline(c(0,1))

# Now let's combine the two separate models
plot((predict(p1)+predict(p2))/2, d$y,xlab="Predicted", ylab="Observed", main="Combination of Models for X1 and X2", pch=16)
abline(c(0,1))

# Get the R2 for that "combined" model
pred = (predict(p1)+predict(p2))/2
act = d$y
1 - sum((pred-act)^2) / sum((act- mean(act))^2)


# Now for the proper joint model.
p = lm(y~., d)
summary(p)
# Plot it
plot(predict(p), d$y,xlab="Predicted", ylab="Observed", main="Model with X1 and X2", pch=16)
abline(c(0,1))
