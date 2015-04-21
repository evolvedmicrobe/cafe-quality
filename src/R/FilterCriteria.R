library(ggplot2)
library(scales)
# Generate simulated data
accuracies =rep(seq(.99, 1.0, .001), 200)
lengths = rnorm(length(accuracies), 1000,50)
empirical = runif(length(accuracies), accuracies*.85,1)
data = data.frame(Acc = accuracies, Length= lengths, Empirical = empirical)


# BELOW TAKES A STRUCTURE LIME THE SIMULATED DATA AND MAKES A PLOT
d = read.csv("/Users/nigel/git/cafe-quality/NotTracked/master_5ba5286_combined_reads.csv")
d = d[!is.na(d$Reference),]
d$Empirical = 1 - d$NumErrors/d$Length
d$Empirical[is.na(d$Empirical)] = 0.999
d$Acc=d$PredictedCCSAccuracy
data = d
#data = d[ d$Acc>.99 & d$Reference=="lambda_NEB3011" ,]
#Function for confidence intervals
getCI = function(lengths, predictedAccuracy) {
  #simulate from each read 5 times
  rep = 5
  rlength=rep(lengths, rep)
  lambda = rlength * (1-predictedAccuracy)
  n = length(rlength)
  simAccuracy = rpois(n, lambda) / rlength
  #now get 95% range
  1-quantile(simAccuracy, c(0.05,0.95))
}
getRow <- function(predictedAccuracy) { c(predictedAccuracy, getCI(data$Length[data$Acc == predictedAccuracy], predictedAccuracy))}

#make data to plot
tp = data.frame(t(sapply(unique(data$Acc), getRow)))
colnames(tp) <- c("Accuracy", "low", "high")


ggplot(data,aes(x=Acc, y=Empirical, color=Reference)) + geom_point() + theme_bw() + 
  scale_y_continuous(labels = percent) + scale_x_continuous(labels=percent) +
  geom_smooth(data=tp,aes(x=Accuracy,ymin=low, ymax=high, y=Accuracy), stat="identity", color="blue") +
  labs(x="Predicted Accuracy", y="Actual Accuracy", title="Accuracy by CCS Read Score (Scatter)") + scale_y_continuous(limits=c(.85,1)) +

