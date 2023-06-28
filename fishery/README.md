# Fishery Example
This folder contains the code for the replication of the fish catching example in Pratola, Chipman, and McCulloh (2020). It models the live-weight of fish caught by fishery boats on specific days.

## Data
The Fishery-Bart.csv file contains the data for this example. It consists of characterstics of the boat and the outcome variable.

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
| `overallsd` | $\lambda=\sqrt{sd(y_{i}:i\in I_{train})^2Q_{0.1}(f_{\chi^2(\nu)})/\nu}$ | Scale parameter for the prior of the variance model* |
| `pbd` | $(0.7, 0.0)$ | Probability of birth/death for homoscedastic model |

Note that $I_{train}$ denotes the set of indices used in estimating the model.

$*$ These quantities are still altered in the R packages, as in Pratola et al. (2020) section 3.4.

### Evaluation
Pratola et al. (2020) use several evaluation methods, including root mean squared error (RMSE), $e$-statistics, qq-plots, and H-evidence plots. For an explanation of these metrics and graphics, please see the original paper. These evaluation metrics/graphics are created as:
| Metric/Graphic | Description |
| --- | --- |
| RMSE | $\sqrt{\sum_{i\in I_{test}}(\widehat{y_{i}}-y_{i})^2}$ |
| $e$-statistic | Package `energy` |
| qq-plot | Plots the quantiles from the posterior against 10,000 random uniform draws |
| H-evidence | Obtains quantiles from the posterior and sorts these on the value of the mean |

Note that $I_{test}$ denotes the set of indices used in estimating the model.
