### Preamble #############################
# In this script I run the analysis for the fish catch production example.
# It generates some evaluation plots for the thesis including H-evidence
# and predictive qq-plots.

### Resources #############################
library(rbart)
library(energy)
source("../qinsamp.R")

### Data #############################
file_name <- "Fishery-Bart.csv" # Set correct file name
df <- read.csv(file_name) # Read the file into a dataframe

set.seed(99)
idxtrain <- sample(1:6806, 4084) # Random sample for training rows
x <- rbartModelMatrix(df[,!names(df) %in% c("y")]) # Load covariates into matrix

# Training data
xtrain <- x[idxtrain,]
ytrain <- df$y[idxtrain]

# Prediction Data
xpred <- x[-idxtrain,]
ypred <- df$y[-idxtrain]

### Settings #############################
# MCMC
nburn <- 1000 # Number of draws discarded for burn-in
niter <- 2000 # Number of iterations for the MCMC after burn-in
nthreads <- 5 # Number of OpenMP threads to use in computation

# Dimensions
ntrain <- NROW(ytrain) # Number of simulation points for training
npred <- NROW(ypred) # Number of simulation points for prediction
ncov <- NCOL(xtrain) # Number of covariates

# (H)BART parameters
ntree <- 200 # Number of trees in the mean model
ncut <- 1000 # Number of cutpoints for each predictor
k <- 5 # Prior hyperparameter for mean model

# Variance prior (3.4; Pratola et al)
nu <- 3
lambda <- sqrt((sd(ytrain)^2*qchisq(0.1, nu))/nu)

### Model Fitting #############################
# Heteroscedastic BART
hbart_train <- rbart(xtrain, ytrain, xpred, ntree = ntree, nskip = nburn, ndpost = niter, numcut = ncut, k = k, tc = nthreads, overallsd = lambda, overallnu = nu)
hbart_pred <- predict.rbart(hbart_train, x.test = xpred)

# Homoscedastic BART
pbd <- c(0.7,0.0) # Sets probability of death/birth for variance trees to 0.
bart_train <- rbart(xtrain, ytrain, xpred, ntree = ntree, ntreeh = 1, pbd = pbd, nskip = nburn, ndpost = niter, numcut = ncut, k = 2, tc = nthreads, overallsd = lambda, overallnu = nu)
bart_pred <- predict.rbart(bart_train, x.test = xpred)

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
hbart_eval$estat$val <- as.numeric(edist(matrix(c(hbart_eval$estat$qnt, runif(10000)), ncol = 1), c(npred, 10000)))
set.seed(99)
bart_eval$estat$val <- as.numeric(edist(matrix(c(bart_eval$estat$qnt, runif(10000)), ncol = 1), c(npred, 10000)))

### Data Plots #############################
par(mar = c(4, 5, 1, 1), cex.lab = 1.3)
boxplot(y ~ x10, data = df, ylab = "Live-Weight of Fish Caught (kg)", xlab = "Mesh Size (mm)", outline = FALSE)
boxplot(df$y[df[,1] == 1], df$y[df[,2] == 1], df$y[df[,3] == 1], df$y[df[,4] == 1],
        names = c("3L", "3M", "3N", "3O"),
        xlab = "Zone",
        ylab = "Live-Weight of Fish Caught (kg)",
        outline = FALSE)
par(mfrow = c(1,4), oma = c(2,4,0,0), mar = c(5,2,2,1))
plot(y ~ x10, data = df[df[,1] == 1,], ylab = "", xlab = "", col = alpha("black", 0.1), cex.axis = 1.5, cex.main = 1.5)
title("Zone 3L")
mtext("Live-Weight of Fish Caught", adj = .5, line = 3, side = 2)
plot(y ~ x10, data = df[df[,2] == 1,], ylab = "", xlab = "", col = alpha("black", 0.1), cex.axis = 1.5, cex.main = 1.5)
title("Zone 3M")
mtext("Mesh Size (mm)", adj = 2.5, line = 3, side = 1)
plot(y ~ x10, data = df[df[,3] == 1,], ylab = "", xlab = "", col = alpha("black", 0.1), cex.axis = 1.5, cex.main = 1.5)
title("Zone 3N")
plot(y ~ x10, data = df[df[,4] == 1,], ylab = "", xlab = "", col = alpha("black", 0.1), cex.axis = 1.5, cex.main = 1.5)
title("Zone 3O")

dev.off()

### H-Evidence Plot #############################
par(mar = c(4, 5, 1, 1))
qnt <- c(0.05, 0.95) # The lower- and upper quantiles for the interval of the posterior
oid <- order(hbart_pred$smean) # Indices in order of variance posterior mean
post_interval <- apply(hbart_pred$sdraws, 2, quantile, probs = qnt)
# Initialize plot panel
plot(range(hbart_pred$smean),
     range(post_interval),
     type = "n",
     xlab = expression(hat(s)(x)),
     ylab = "s(x) posterior",
     cex.axis = 1.3,
     cex.lab = 1.5)
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
par(mfrow = c(1,2), oma = c(0,4,0,0), mar = c(5,2,2,1))
qqplot(hbart_eval$estat$qnt, runif(10000), col = "grey", cex.axis = 1.3, xlab = "", ylab = "")
abline(0, 1, col = "black", lwd = 2) # Add 45 degree ray
mtext("Uniform Draws", adj = .5, line = 3, side = 2, cex = 1.5)
mtext("Data Quantiles", adj = 1.7, line = 3, side = 1, cex = 1.5)
title("HBART")

# Homoscedastic BART
qqplot(bart_eval$estat$qnt, runif(10000), col = "grey", cex.axis = 1.3, xlab = "", ylab = "")
abline(0, 1, col = "black", lwd = 2) # Add 45 degree ray
title("BART")
