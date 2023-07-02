### Preamble #############################
# In this script I perform the main analysis for the extension of the thesis.
# The data used here comes from Preprocessing.m and the output is analysed by
# Postprocessing.m. This script also includes some boxplots and an H-evidence
# plot.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Sets working directory to current file directory

### Resources #############################
remove.packages("rbart") # Make sure the existing rbart package is unloaded
# Comment this out if package was already built
setwd("../surv_rbart")
# install.packages(devtools) # Make sure this is installed
library(devtools)
devtools::install() # Build and install the package
library(rbart)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Set working directory back again

### Data #############################
file_name <- "Bankruptcy-Bart.csv" # Set correct file name
df <- read.table(file_name, sep = ",", header = TRUE) # Read the file into a dataframe

# Convert time values
df[, "t"] <- df[, "t"] %% 100
df <- df[complete.cases(df), ]
df <- df[df[ , "t"] < 21,]

# Extract key values
keys <- unique(df[, "gvkey"])
bkeys <- unique(df[df[, "y"] == 1, "gvkey"])
keytrain <- unique(df[df[, "train"] == 1, "gvkey"])
keypred <- setdiff(keys, keytrain)
idxtrain <- which(df[, "gvkey"] %in% keytrain) # Random sample for training rows
idxpred <- which(df[, "gvkey"] %in% keypred) # Random sample for pred rows
x <- rbartModelMatrix(df[,!names(df) %in% c("gvkey", "train", "y", "status")]) # Load covariates into matrix

# Exploratory Analysis
par(mfrow = c(2,4), oma = c(2,2,0,0), mar = c(2,2,2,1))
boxplot(pe_inc ~ status, data = df, outline = FALSE, xlab = "", names = c("Viable", "Bankrupt"), cex.axis = 1.3) 
title("Price/Earnings")
boxplot(bm ~ status, data = df, outline = FALSE, xlab = "", names = c("Viable", "Bankrupt"), cex.axis = 1.3)
title("Book/Market")
boxplot(cash_lt ~ status, data = df, outline = FALSE, xlab = "", names = c("Viable", "Bankrupt"), cex.axis = 1.3)
title("Cash Balance/Total Liabilities")
boxplot(cash_debt ~ status, data = df, outline = FALSE, xlab = "", names = c("Viable", "Bankrupt"), cex.axis = 1.3)
title("Cash Flow/Total Debt")
boxplot(ps ~ status, data = df, outline = FALSE, xlab = "", names = c("Viable", "Bankrupt"), cex.axis = 1.3)
title("Price/Sales")
boxplot(lt_debt ~ status, data = df, outline = FALSE, xlab = "", names = c("Viable", "Bankrupt"), cex.axis = 1.3)
title("Long-term Debt/Total Liabilities")
boxplot(roe ~ status, data = df, outline = FALSE, xlab = "", names = c("Viable", "Bankrupt"), cex.axis = 1.3)
title("Return on Equity")
boxplot(short_debt ~ status, data = df, outline = FALSE, xlab = "", names = c("Viable", "Bankrupt"), cex.axis = 1.3)
title("Short-term Debt/Total Debt")

# Training data
xtrain <- x[idxtrain,]
ytrain <- df$y[idxtrain]

# Prediction Data
xpred <- x[idxpred,]
ypred <- df$y[idxpred]

### Settings #############################
# MCMC
nburn <- 500 # Number of draws discarded for burn-in
niter <- 2000 # Number of iterations for the MCMC after burn-in
nthreads <- 15 # Number of OpenMP threads to use in computation

# Simulation
ntrain <- NROW(ytrain) # Number of simulation points for training
npred <- NROW(ypred) # Number of simulation points for prediction
ncov <- NCOL(xtrain) # Number of covariates

# (H)BART parameters
ntree <- 200 # Number of trees in the mean model
ncut <- 1000 # Number of cutpoints for each predictor
k <- 1 # Prior hyperparameter for mean model

# Variance prior (3.4; Pratola et al)
nu <- 3
lambda <- sqrt(qchisq(0.1, nu)/nu)

### Model Fitting #############################
# Heteroscedastic BART
hbart_train <- rbart(xtrain, ytrain, xpred, ntree = ntree, nskip = nburn, ndpost = niter, numcut = ncut, k = k, tc = nthreads, overallsd = lambda, overallnu = nu)
hbart_pred <- predict.rbart(hbart_train, x.test = xpred)

# Homoscedastic BART
pbd <- c(0.7,0.0) # Sets probability of death/birth for variance trees to 0.
bart_train <- rbart(xtrain, ytrain, xpred, ntree = ntree, ntreeh = 1, pbd = pbd, nskip = nburn, ndpost = niter, numcut = ncut, k = k, tc = nthreads, overallsd = lambda, overallnu = nu)
bart_pred <- predict.rbart(bart_train, x.test = xpred)
plot(pnorm(hbart_pred$mmean/hbart_pred$smean))
### Survival Analysis #############################
# Initialize data frames
heter_surv <- data.frame(matrix(nrow = length(hbart_pred$mmean), ncol = 5))
homo_surv <- data.frame(matrix(nrow = length(hbart_pred$mmean), ncol = 5))
colnames(heter_surv) <- c("gvkey", "t", "y", "p", "S")
colnames(homo_surv) <- c("gvkey", "t", "y", "p", "S")

# Set the company keys
heter_surv[,"gvkey"] <- df[df[,"gvkey"] %in% keypred, "gvkey"]
homo_surv[,"gvkey"] <- df[df[,"gvkey"] %in% keypred, "gvkey"]

# Loop over companies
for (key in keypred) {
  # Obtain rows corresponding to company key
  keyrows <- which(heter_surv[,"gvkey"] == key)
  
  # Set time
  heter_surv[keyrows, "t"] <- df[df[,"gvkey"] == key, "t"]
  homo_surv[keyrows, "t"] <- df[df[,"gvkey"] == key, "t"]
  
  # Set bankruptcy status
  heter_surv[keyrows, "y"] <- df[df[, "gvkey"] == key, "y"]
  homo_surv[keyrows, "y"] <- df[df[, "gvkey"] == key, "y"]
  
  # Set hazard rate
  heter_surv[keyrows, "p"] <- pnorm(hbart_pred$mmean[keyrows]/hbart_pred$smean[keyrows])
  homo_surv[keyrows, "p"] <- pnorm(bart_pred$mmean[keyrows]/bart_pred$smean[keyrows])
  
  # Set survival rate
  for (i in keyrows) {
    heter_surv[i,"S"] <- prod(1 - heter_surv[keyrows[1]:i,"p"])
    homo_surv[i,"S"] <- prod(1 - homo_surv[keyrows[1]:i,"p"])
  }
}
plot(heter_surv[heter_surv$gvkey == keypred[10], "S"], ylim = c(0, 1))

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

# I already wrote these, but feel free to overwrite.
# write.csv(heter_surv, "heter_surv.csv")
# write.csv(homo_surv, "homo_surv.csv")
