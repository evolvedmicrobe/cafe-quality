library(unitem)
library(nnet)
set.seed(6178) #Set seed for reproducible data
bps = c("A", "C", "G", "T")
eps = 0.01

# Function to create a fake template
generateTemplate <- function(size) { 
  # Make a fake model sequence
  samp = sample(bps, size, replace=TRUE)
  read = paste(samp, sep="", collapse="")
  # fake SNR values at each point
  snr = sample(seq(1,size), size) + rnorm(size, 0, 2)  
  # Just hard coding transition parameters for now
  Probs = matrix(rep(c(.8,.05,.05,.1), size), ncol=4, byrow=TRUE)
  list(read = read, covar = data.frame(BP=factor(samp), SNR=snr), Probs = Probs)
}

sampleMatch <- function(bp) {
  if (runif(1) >eps) {return(bp)} else {
    nbp = bp
    while(nbp == bp) {
      nbp = sample(bps, 1)
    }
    return(nbp)
  }
}

sampleStick <-function(bp) {
  nbp = bp
  while(nbp == bp) {
    nbp = sample(bps, 1)
  }
  return(nbp)
}

simulateOutcome <- function(d) {
  n = nrow(d$covar)
  model = d$covar$BP
  prs = d$Probs
  i=1
  res = sampleMatch(as.character(model[1]))
  while(i < n) 
  {
    probs = prs[i,]
    outcome = sample(1:4, 1, prob = probs)
    cbp = as.character(model[i])
    if (outcome == 1) {
      res = c(res, sampleMatch(cbp))
      i = i + 1
    } else if (outcome == 2) { # branch
      res = c(res, cbp)
    } else if (outcome==3) {
      res = c(res, sampleStick(cbp))
    }
    else { i = i + 1} # deletion    
  }
  return(paste(res, sep="", collapse=""))
  
}

# IGNORE THIS FOR NOW
generateParameterss <- function(d) {
  # Match, stick, branch
  #m = model.frame(~ BP + SNR + I(SNR^2), t$covar)
  #Terms <- attr(m, "terms")
  #X2 <- model.matrix(Terms, m, contrasts)
}



# Generate one fake read/template pair
d = generateTemplate(500)
m = d$read
o = simulateOutcome(d)

# Now fit it with EM
ll = -Inf 
nll = 0
end_dif = 1e-3 # End tolerance
t = d$Probs # Transition probabilities
s = eps # 
# E-M main loop
while (abs(1 - ll/nll) > end_dif) {
  print(ll)
  ll = nll
  r = FitUnitEM(t, m, o, s) 
  fit = multinom(r$pc~d$covar$BP) # Multinomial regression step, Need to better specify covariates in future, just using BP for now.
  t = predict(fit, type="probs") # New transition parameters
  s = r$sn / r$sd
  nll = r$ll  
}

