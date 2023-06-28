# Alcohol Example
This folder contains the code for the replication of the individual alcohol consumption example in Pratola, Chipman, and McCulloh (2020). It models the number of alcohol consumptions of the last two weeks.

## Data
The Alcohol-Bart.csv file contains the data for this example. It consists of demographics, the outcome variable, and an indicator for if the individual received a physician's advice to limit alcohol consumption.

## Code
This section contains only important code snippets; the code is supplied with comprehensive comments for a more detailed overview.

### Parameters
As this is a replication, the same parameters will be used as in Pratola et al. (2020). These are:

| (Hyper)Paramater | Value | Description |
| --- | --- | --- |
| `nskip` | $1,000$ | The number of draws discarded for convergence of the Gibbs sampler |
| `ndpost` | $2,000$ | The number of draws used as posterior distribution |
| `ntree` | $200$  | The number of trees in the mean model | 
| `ntreeh` | $40$ | The number of trees in the variance model |
| `ncut` | $1,000$ | The number of cutpoints for each predictor |
| `k` | $\kappa = 2$ | Prior hyperparameter for the mean model |
| `overallnu` | $\nu=3$ | Degrees of freedom for the prior of the variance model* |
| `overallsd` | $\lambda=\sqrt{Q_{0.1}(f_{\chi^2(\nu)})/\nu}$ | Scale parameter for the prior of the variance model* |
| `pbd` | $(0.7, 0.0)$ | Probability of birth/death for homoscedastic model |
