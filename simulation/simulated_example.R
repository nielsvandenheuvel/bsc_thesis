### Preabmle #############################
# This script provides a simulation example for the HBART methodology.
# It produces evaluation and data plots for the thesis. The DGP is derived
# from Pratola et al. (2020).
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Sets working directory to current file directory
### Resources #############################
# remove.packages("rbart") # Remove if previously ran bankrupt.R
install.packages("rbart")
library(rbart)
library(energy)
source("../qinsamp.R")

# Set if workspace exists already
workspace_path <- ""
# Set if you want to save the workspace after completion
workspace_save_path <- ""

if (workspace_path == "") {
  ### Settings #############################
  # MCMC
  nburn <- 1000 # Number of draws discarded
  niter <- 2000 # Number of iterations for the MCMC after burn-in
  nthreads <- 5 # Number of OpenMP threads to use in computation
  
  # Simulation
  set.seed(27)
  ntrain <- 500 # Number of simulation points for training
  npred <- 500 # Number of simulation points for prediction
  ncov <- 1 # Number of covariates
  
  # (H)BART parameters
  ntree <- 200 # Number of trees in the mean model
  ncut <- 1000 # Number of cutpoints for each predictor
  k <- 5 # Prior hyperparameter for mean model
  
  ### Simulate Data #############################
  # Training data
  xtrain = matrix(sort(runif(ntrain*ncov)), ncol = ncov) # Uniform predictors
  ytrain = 4*xtrain[,1]^2 + 0.2*exp(2*xtrain[,1])*rnorm(ntrain)
  
  # Prediction data
  xpred = matrix(sort(runif(npred*ncov)), ncol = ncov) # Uniform predictors
  ypred = 4*xpred[,1]^2 + 0.2*exp(2*xpred[,1])*rnorm(npred)
  
  ### Model Fitting #############################
  # Heteroscedastic BART
  hbart_train = rbart(xtrain, ytrain, xpred, ntree = ntree, nskip = nburn, ndpost = niter, numcut = ncut, k = k, tc = nthreads)
  hbart_pred = predict.rbart(hbart_train, x.test = xpred)
  
  # Homoscedastic BART
  pbd <- c(0.7,0.0) # Sets probability of death/birth for variance trees to 0.
  bart_train = rbart(xtrain, ytrain, xpred, ntree = ntree, ntreeh = 1, pbd = pbd, nskip = nburn, ndpost = niter, numcut = ncut, k = 2, tc = nthreads)
  bart_pred = predict.rbart(bart_train, x.test = xpred)
  
  ### Quality Metrics #############################
  hbart_eval <- list(rmse = NaN,
                     estat = list(val = NaN, pdraw = NaN, qnt = NaN))
  bart_eval <- list(rmse = NaN,
                     estat = list(val = NaN, pdraw = NaN, qnt = NaN))
  
  # Root mean squared error
  hbart_eval$rmse <- sqrt(mean((ypred - hbart_pred$mmean)^2))
  bart_eval$rmse <- sqrt(mean((ypred - bart_pred$mmean)^2))
  
  # e-statistic
  hbart_eval$estat$pdraw <- hbart_pred$mdraws + hbart_pred$sdraws*matrix(rnorm(niter*npred), nrow = niter)
  bart_eval$estat$pdraw <- bart_pred$mdraws + bart_pred$sdraws * matrix(rnorm(niter*npred), nrow = niter)
  hbart_eval$estat$qnt <- qsamp(ypred, hbart_eval$estat$pdraw)
  bart_eval$estat$qnt <- qsamp(ypred, bart_eval$estat$pdraw)
  set.seed(99)
  hbart_eval$estat$val <- as.numeric(edist(matrix(c(hbart_eval$estat$qnt, runif(1000)), ncol = 1), c(npred, 1000)))
  set.seed(99)
  bart_eval$estat$val <- as.numeric(edist(matrix(c(bart_eval$estat$qnt, runif(1000)), ncol = 1), c(npred, 1000)))
  
  # Save workspace
  if (workspace_save_path != "") {
    save(file = workspace_save_path)
  }
} else {
  load(workspace_path)
}

### H-Evidence Plot #############################
qnt <- c(0.05, 0.95) # The lower- and upper quantiles for the interval of the posterior
oid <- order(hbart_pred$smean) # Indices in order of variance posterior mean
post_interval <- apply(hbart_pred$sdraws, 2, quantile, probs = qnt)
# Initialize plot panel
plot(range(hbart_pred$smean),
     range(post_interval),
     type = "n",
     xlab = expression(hat(s)(x)),
     ylab = "s(x) posterior",
     cex.lab = 1.5,
     cex.axis = 1.3)
# Plot the intervals for each posterior mean of the variance
for (i in 1:npred) {
  lines(rep(hbart_pred$smean[oid[i]],2), post_interval[,oid[i]],col = "gray")
}
abline(h = bart_pred$smean[1], col = "black", lwd = 2) # Plot BART variance
# Plot confidence bounds for BART variance
abline(h = quantile(bart_pred$sdraws,0.05), col = "black", lty = 2)
abline(h = quantile(bart_pred$sdraws,0.95), col = "black", lty = 2)

### Predictive QQ-Plots #############################
# Heteroscedastic BART
par(mfrow = c(1,2), oma = c(0,4,0,0), mar = c(5,2,2,1), cex.lab = 1.3)
qqplot(hbart_eval$estat$qnt, runif(10000), col = "grey", xlab = "", ylab = "")
abline(0, 1, col = "black", lwd = 2) # Add 45 degree ray
mtext("Uniform Draws", adj = .5, line = 3, side = 2, cex = 1.3)
mtext("Data Quantiles", adj = 1.5, line = 3, side = 1, cex = 1.3)
title("HBART")

# Homoscedastic BART
qqplot(bart_eval$estat$qnt, runif(10000), col = "grey", xlab = "", ylab = "")
abline(0, 1, col = "black", lwd = 2) # Add 45 degree ray
title("BART")


### Data Plots #############################
# The DGP can be expressed as y = fy + sy*Z, where Z~N(0,1).
fytrain <- 4*xtrain[,1]^2
sytrain <- 0.2*exp(2*xtrain[,1])
fypred <- 4*xpred[,1]^2
sypred <- 0.2*exp(2*xpred[,1])

# Plot of true prediction data
plot(xpred, ypred, col = "grey", xlab = "X", ylab = "Y")
lines(xpred, fypred, col = "black", lwd = 3)
lines(xpred, fypred + 2*sypred, col = "black", lwd = 1, lty = 2)
lines(xpred, fypred - 2*sypred, col = "black", lwd = 1, lty = 2)

# Plot of true data and estimates
plot(xpred, ypred, col = "grey", pch = 1, xlab = "X", ylab = "Y")
# True
lines(xpred, fypred, col = "grey", lwd = 3, lty = 2) #true f
lines(xpred, fypred + 2*sypred, col = "grey", lwd = 3,lty = 2)
lines(xpred, fypred - 2*sypred, col = "grey", lwd = 3,lty = 2)
# Estimated
lines(xpred, hbart_pred$mmean, col = "black", lwd = 2, lty = 1) #estimate of f
lines(xpred, hbart_pred$mmean + 2*hbart_pred$smean, col = "black", lwd = 2, lty = 1)
lines(xpred, hbart_pred$mmean - 2*hbart_pred$smean, col = "black", lwd = 2, lty = 1)

# Mean plot
rga <- range(c(hbart_pred$mmean, hbart_pred$m.lower, hbart_pred$m.upper))
plot(range(xpred), rga, xlab = "X", ylab = "f(X)", type = "n")
lines(xpred, fypred, col = "gray", lwd = 2, lty = 2) 
lines(xpred, hbart_pred$mmean, col = "black", lwd = 2, lty = 1) 
lines(xpred, hbart_pred$m.lower, lwd = 2, lty = 4, col = "black")
lines(xpred, hbart_pred$m.upper, lwd = 2, lty = 4, col = "black")

# Std. plot
rga <- range(c(hbart_pred$smean, hbart_pred$s.lower, hbart_pred$s.upper))
plot(range(xpred), rga, xlab = "x", ylab = "s(x)", type = "n")
lines(xpred, sypred,col = "gray", lwd = 2, lty = 2) # True
lines(xpred, hbart_pred$smean, col = "black", lwd = 2, lty = 1) # Estimated
# Credible intervals
lines(xpred, hbart_pred$s.lower, lwd = 2, lty = 4, col = "black") 
lines(xpred, hbart_pred$s.upper, lwd = 2, lty = 4, col = "black")

### MCMC Convergence Plots #############################
# Mean traceplot
par(mfrow = c(1, 2))
xplot <- c(0.02, 0.42, 0.63, 0.79, 0.91)
cols <- c("azure4", "aquamarine4", "darkgoldenrod", "blueviolet", "coral3")

# Mean traceplots
plot(x = NaN, ylim = c(-0.3, 4), xlim = c(0, 2000), xlab = "", cex.axis = 1.3)
for (i in 1:length(xplot)) {
  row <- which.min(abs(xtrain - xplot[i]))
  lines(hbart_pred$mdraws[,row], type = "l", col = cols[i], lwd = 1)
  lines(lowess(hbart_pred$mdraws[,row]), col = "black", lwd = 3)
}
mtext("Iteration", adj = 1.3, line = 3, side = 1, cex = 1.5)
title("Mean")

# Variance traceplots
plot(x = NaN, ylim = c(0.2, 1.5), xlim = c(0, 2000), xlab = "", cex.axis = 1.3)
for (i in 1:length(xplot)) {
  row <- which.min(abs(xtrain - xplot[i]))
  lines(hbart_pred$sdraws[,row], type = "l", col = cols[i], lwd = 1)
  lines(lowess(hbart_pred$sdraws[,row]), col = "black", lwd = 3)
}
title("Standard Deviation")
