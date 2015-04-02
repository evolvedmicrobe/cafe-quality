library(ggplot2)
library(gridExtra)
library(grid)
library(unitem)
library(nnet)
library(gtools)
library(reshape)

setwd("/Users/nigel/git/cafe-quality/NotTracked/")
#setwd("/Users/nigel/git/libunitem/EM_Temp")

set.seed(6178) #Set seed for reproducible data

#### LOAD DATA

if(!exists("samps")) {
  load("job219390.rda")
  samps = sample(data, 2000)
}
######  END DATA LOAD 


# Filter out large descrepencies
# The idea is that if a read/template pair differ by more than 20% in length, they are not
# suitable for fitting.
noGood <- function(x) {
  nm = nchar(x$model)
  no = nchar(x$outcome)
  dif = abs(1 - nm/no)
  dif < .2
}
samps = Filter(noGood, samps)





### FORMULA - Simply change this line to test a new model
fmla = formula(" ~ nextIQV + nextDQV +  I(nextIQV^2) + nextDQV + nextDQV:nextMQV +  nextMQV + currMQV + snrC + I(snrC^2) + snrG + I(snrG^2)")

### END FORMULA

bps = c("A", "C", "G", "T")


# Function to remove the last BP of a read/template pair
cutLast <- function(x) {
  substr(x, 0, nchar(x)-1)
}

# Function to convert loaded data into objects for modeling
currentPos <-1
convert <- function(x) {
  currentPos <<- currentPos+1
  #print(x$model)
  #print(x$outcome)
  covars = data.frame(x$t)
  # convert numerics to numerics (as factors by default)
  for(j in 3:ncol(covars)) {
    covars[,j] = as.numeric(as.character(covars[,j]))
  }
  covars$CTX = covars$currBP:covars$nextBP
  x$covars = covars
  x$t <- NULL
  x$model = cutLast(x$model)
  x$outcome = cutLast(x$outcome)
  m = model.frame(fmla, covars)
  Terms <- attr(m, "terms")
  X2 <- model.matrix(Terms, m, contrasts, na.action=na.fail) # this is the covariation matrix
  x$covars = X2
  x$pc <- rdirichlet(nrow(X2), 40*c(.85, .05, .05, .05)) # initialize the transition probabilities to random data.
  return(x)
}

n_data = lapply(samps, convert)


# Make covariance matrix and a list of pseudo counts to start things off
covars = do.call(rbind, sapply(n_data, function(x) x$covars)) # Call concatenates matrix
head(covars)
covars = as.matrix(covars)
eps = 0.01

# Softmax conversion function
toProb <- function(x) {
  b = c(1, exp(x))
  b / sum(b)
}

estimated <- function(read, fit) {
  # Predict new transition probabilities
  # Note, this simple thing below should have worked, had to roll my own.
  #predict(fit, type="probs", newdata = read$covars) # New transition parameters
  co = coef(fit) # get coefficients
  x = t(read$covars) # Get covariates
  o = t(co%*%x) # Matrix multiply
  read$t = t(apply(o, 1, toProb)) # Soft max conversion 
  return(read)
}


fit <- function(read, eps) {
  # Fit a new model and update the pseudo counts
  read$r = FitUnitEM(read$t, read$model, read$outcome, eps )
  read$pc = read$r$pc
  return (read)
}


######   MAIN E-M LOOP  ########
tot = sum(sapply(n_data, function(x) nrow(x$covars)))
ptm <- proc.time()
ll = -Inf 
nll = -Inf
end_dif = 1e-5 # End relative tolerance
c_wts = NULL # The last regression parameters
while ( is.infinite(ll) | abs(1 - ll/nll) > end_dif) {
  print(ll)
  ll = nll # set the old likelihood equal to the new.
  # Concatenate all the read/template pseudo counts into one matrix
  outcomes = do.call(rbind, sapply(n_data, function(x) x$pc))
  # M step, fit parameters
  if( is.null(c_wts)) { # Fit it, using last parameters
    cfit = multinom(outcomes~covars-1) # Multinomial regression step, Need to better specify covariates in future, just using BP for now.
  } else {
    ## A note on this:
    # The multinom package is from Brian Ripley's Neural Network package.
    # It fits the "multinomial regression" by making the problem equivalent to a 
    # neural network with no middle layers, just inputs that directly feed into the output 
    # network.  It's a bit unclear how the target function works, but it seems to be working as 
    # the likelihood is always going up and the softcounts are akin to the functions somm=2 method.
    # The backpropogation uses BFGS and so is likely not as efficient as Fisher Scoring.
    # I set the decay = 0 as we do not scale the covariates
    # Wts is the start parameters.  
    # 
    # In the future it might be worth trying out either
    #     1 - Using the mlogit package with weights on the observations for faster/guaranteed-global convergence
    #     2 - Trying out a more complete Neural Network here, the code is already there for it!
    cfit = multinom(outcomes~covars-1, Wts = c_wts, decay=0) # Multinomial regression step, Need to better specify covariates in future, just using BP for now.
    
  }
  c_wts = cfit$wts
  # update parameters
  n_data = lapply(n_data, function(x) estimated(x, cfit))
  # update the pseudo counts
  n_data = lapply(n_data, function(x) fit(x, eps))
  # E step, update parameters
  num = sapply(n_data, function(x) x$r$sn)
  denom = sapply(n_data, function(x) x$r$sd)
  eps = sum(num) / sum(denom)  
  nlls = sapply(n_data, function(x) x$r$ll)
  nll = sum(nlls)  
  if (nll < ll) {
    stop("Algorithmic Error!")
  }
}

# Return the model
print(cfit)
p = length(coef(cfit))
aic = 2 * p - 2 * nll
print(fmla)
print(paste("LL =", nll, " P =", p, " AIC =", aic))



