# Bankruptcy Example
This folder contains the code for the extension of SHBART to bankruptcy prediction.
## Data


## Code
This section contains only important code snippets; the code is supplied with comprehensive comments for a more detailed overview.

### Parameters
In the extension, the following settings are used:

| (Hyper)Paramater | Value | Description |
| --- | --- | --- |
| `nskip` | $500$ | The number of draws discarded for convergence of the Gibbs sampler |
| `ndpost` | $2,000$ | The number of draws used as posterior distribution |
| `nthreads` | $15$ | The number of OpenMP threads to use in computation | 
| `ntree` | $200$  | The number of trees in the mean model | 
| `ntreeh` | $40$ | The number of trees in the variance model |
| `ncut` | $1,000$ | The number of cutpoints for each predictor |
| `k` | $\kappa = 1$ | Prior hyperparameter for the mean model |
| `overallnu` | $\nu=3$ | Degrees of freedom for the prior of the variance model* |
| `overallsd` | $\lambda=\sqrt{Q_{0.1}(f_{\chi^2(\nu)})/\nu}$ | Scale parameter for the prior of the variance model* |
| `pbd` | $(0.7, 0.0)$ | Probability of birth/death for homoscedastic model |

Note that $I_{train}$ denotes the set of indices used in estimating the model.

$*$ These quantities are still altered in the R packages, as in Pratola et al. (2020) section 3.4.

### Evaluation
Several evaluation statistics and graphics are employed in this example. These include Brier scores, H-evidence plot, and ROC curve. These quantities/graphics are created as:
| Metric/Graphic | Description |
| --- | --- |
| Brier score | $BS_{j}=\sqrt{\sum_{i\in I_{test}}(\widehat{y_{ij}}-y_{ij})^2}$ |
| H-evidence | Obtains quantiles from the posterior and sorts these on the value of the mean |
| ROC curve | Please see the thesis (Section 3.6) for a comprehensive explanation. |

Note that $I_{test}$ denotes the set of indices used in estimating the model. Moreover, these quantities/graphics are created in the file `Postprocessing.m`.

---

## Elaborations

### Preprocessing.m
**Helper Functions**

I created some helper function for in-line functional programming. The `iff` function works as an in-line if statement.
```matlab
iff = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();
```
The curly function allows us to access created cell arrays in-line.
```matlab
curly = @(x, varargin) x{varargin{:}};
```

**Loading Data**

For the bankruptcy data, I only need the year in which the firm went bankrupt and the corresponding `gvkey`, which is called `GvkeyBefore` in the bankrutpcy dataset.
```matlab
% Drop variables
bankruptcyData(:, 1:71) = [];
bankruptcyData(:, 73:128) = [];
bankruptcyData(:, 130:end) = [];

% Coarsen DateFiled to YearFiled
bankruptcyData.YearFiled = year(bankruptcyData.DateFiled);
bankruptcyData.DateFiled = [];
```

**Cross-sectional Format**

In order to limit the sample, I drop all the observations before the year 2000.
```matlab
ratiosData(ratiosData{:, "public_date"} < "2000-01-01", :) = [];
```

Then, I can select subsamples for estimation and prediction. I do this by randomly sampling from the keys of the viable firms and add the bankrupt firms to this.
```matlab
rng(0);
gvkeySample = [randsample(gvkeyLive, 1000); gvkeyBankrupt];
```

Before I lose the time-series structure of the covariates, I plot the course of 8 covariates over years until the last observation.
!!!!!!

Now, I can format the data. I loop over all firms and find the corresponding rows in the `ratiosData` table. 
```matlab
ratiosDataRows = find(ratiosData{:, "gvkey"} == gvkeySample(i));
```
Then, I add the `gvkey` and `status` variable, and add the mean of each covariate to the `struct`.
```matlab
for t = 1:length(T)
    if sum(ratiosData.public_date(ratiosDataRows) == T(t)) >= 1
        ratios(t, :) = mean(ratiosData{ratiosDataRows, selectedCovariates}, 1, 'omitnan');
    end
end
```
I can now add the data to the `struct` as a `table` object.
```matlab
data(i).table = array2table([T, y, ratios], "VariableNames", ["t", "y", selectedCovariates]);
```

Finally, I refine the data by deleting viable firms with less than 9 years data and bankrupt firms with less than 4 years data. The bias that may be introduced by this is discussed in the conclusion.
```matlab
filter = arrayfun(@(x) size(x.table, 1) >= 10 | (x.status == 1 & size(x.table, 1) >= 4), data);
filteredData = data(filter);
```

**Missing Data**

To get a brief overview of where the missing data is located, I create a table that contains the fraction of observations missing for each covariate.
```matlab
fracMissing = array2table(arrayfun(@(x) sum(arrayfun(@(i) sum(isnan(i.table{:, x})), filteredData)), 1:70)/sum(arrayfun(@(i) size(i.table, 1), filteredData)));
fracMissing.Properties.VariableNames = filteredData(1).table.Properties.VariableNames;
```
I then add targets for variables that are not informative, likely to be subject to measurement error, or have more than 4% missing data (or a combination of these) to an array which will be used when writing the data back to a `.csv` file.
```matlab
% These variables are often 0 and I expect them to carry little explanatory
% power.
dropMissingVariables = ["rd_sale", "staff_sale", "adv_sale"];

% Add variables that have more than 4% missing data
dropMissingVariables = [dropMissingVariables, fracMissing.Properties.VariableNames(fracMissing{:, :} > 0.04)];
```

I also get rid of bankrupt firms that have some instance of missing data at the year of bankruptcy, as otherwise `complete.cases` in R gets rid of these observations, converting the firm to a viable one.
```matlab
missingBankruptFirms = arrayfun(@(i) iff(i.status == 1, sum(isnan(i.table{size(i.table, 1), ~ismember(i.table.Properties.VariableNames, dropMissingVariables)})), true, 0), filteredData) > 0;
filteredData(missingBankruptFirms, :) = [];
```

**Estimation/Prediction Sample**

To easily identify the observations that go into the estimation or the prediction, I add an identifier that is `1` if an observation is in the estimation sample and `0` if it is in the prediction sample.
```matlab
filteredData(i).train = ismember(i, iTrain);
```
If an observation is in the prediction sample, I add a prediction data table that has complete data for the whole period, as if the firm was viable (in case it is not). This allows me to extract probabilities for the whole period and evaluate the predictions at each point in time.
```matlab
if ~filteredData(i).train
    T = (filteredData(i).table.t(1):2022)';
    covData = repmat(filteredData(i).table{1, 3:end}, length(T), 1);
    filteredData(i).pred = array2table([T, covData], "VariableNames", ["t", selectedCovariates]);
end
```

**Write Back**

At last, I write back the `filteredData` to a table `dataArray`, such that it can be written to a `.csv` file.
```matlab
writetable(dataArray, "Bankruptcy-Bart.csv");
```

---

### Postprocessing.m

**Loading Data**

To allow us to compare the estimations with the actual results, I need to load in the prepared data set `dataArray`. To match the data as used for estimation in R, I drop uncomplete cases and data points after 2020.
```matlab
dataArray = dataArray(~sum(isnan(dataArray{:, :}), 2) > 0, :);
dataArray = dataArray(dataArray.t < 2021, :);
```

**Cross-sectional Format**

Next, I format the data into a struct. This creates an entry for each firm and stores the `gvkey`, the identifier `status`, and the data for HBART `heter` and BART `homo` with the time variable `t` and the outcome variable `y`.
```matlab
data(numFirms) = struct();
for i = 1:numFirms
    dataRows = find(hbartSurv{:, "gvkey"} == gvkey(i));
    
    data(i).gvkey = gvkey(i);
    data(i).status = ismember(gvkey(i), gvkeyBankrupt);
    data(i).heter = hbartSurv(dataRows, 3:6);
    data(i).homo = bartSurv(dataRows, 3:6);
    
    data(i).heter.t = data(i).heter.t + 2000;
    data(i).homo.t = data(i).homo.t + 2000;
    
    data(i).heter.y = dataArray.y(find(dataArray.gvkey == gvkey(i)));
    data(i).homo.y = dataArray.y(find(dataArray.gvkey == gvkey(i)));
end
```

This allows us to make boxplots for some covariates, in order to compare the values for viable and bankrupt firms. The first thing I need to do is make an array that identifies the boxplots for bankrupt (`1`) and viable (`2`) firms. 
```matlab
g = [repmat(1,length(unique(dataArray{ismember(dataArray.gvkey, gvkeyBankrupt), test(i).ratio})),1); repmat(2,length(unique(dataArray{~ismember(dataArray.gvkey, gvkeyBankrupt), test(i).ratio})),1)];
```

Then, I concatenate the data for each of these groups together.
```matlab
x = [unique(dataArray{ismember(dataArray.gvkey, gvkeyBankrupt), test(i).ratio}); unique(dataArray{~ismember(dataArray.gvkey, gvkeyBankrupt), test(i).ratio})];
```

**Time-series Format**

Hereafter, I format the data into a time-series structure, so that I can easily access the data corresponding to each year. Some important lines in this code are the following, which set the status at time `t` to `-1` if the firm is already bankrupt for at least one year.
```matlab
if t == 1
    status = arrayfun(@(x) iff(ismember(panel(t).t, x.heter.t), x.heter.y(x.heter.t == panel(t).t), true, @() -1), data)';
else
    status = arrayfun(@(x) iff(ismember(panel(t).t, x.heter.t), x.heter.y(x.heter.t == panel(t).t), true, @() -1), data)';
    status = arrayfun(@(i) iff(panel(t - 1).data.status(i), -1, true, status(i)), 1:length(status))';
end
```

**ROC curve**

Plotting the ROC curve is not trivial. First, I sort the data on the probability of defaulting and filter out the firms that are already bankrupt.
```matlab
[~, I] = sort(panel(t).data.heter_p);
binData = panel(t).data(I, :);
binData = binData(binData.status >= 0, :);
```

Then, I create 100 bins and compute the number of failed, high-risk firm as a fraction of the total number of failed firms. 
```matlab
[idBin, bins] = discretize(binData.heter_p, 100);
for i = 1:length(bins)
    heter_fracFailed(i, t) = sum(binData.status(idBin <= i))/numBankrupt;
end
```

This is done identically for the results using the BART model.

**Survival Functions**

I have plotted survival functions by simply using the median of the posterior distributions of the induced survival functions. Moreover, I have added shaded regions for 50% and 90% interquantile ranges. This was done using the fill function, which takes in an array that outlines the shaded area. So, the array contains first the lower quantiles from left to right, and then the upper quantiles from right to left, exactly tracing the outline of the shape I want to plot.
```matlab
fill([1:21, fliplr(1:21)], [arrayfun(@(t) quantile(t.data.heter_S(ismember(t.data.gvkey, gvkeyBankrupt)), 0.25), panel), fliplr(arrayfun(@(t) quantile(t.data.heter_S(ismember(t.data.gvkey, gvkeyBankrupt)), 0.75), panel))], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
fill([1:21, fliplr(1:21)], [arrayfun(@(t) quantile(t.data.heter_S(ismember(t.data.gvkey, gvkeyBankrupt)), 0.05), panel), fliplr(arrayfun(@(t) quantile(t.data.heter_S(ismember(t.data.gvkey, gvkeyBankrupt)), 0.95), panel))], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
```

**Brier Scores**

The Brier scores are the equivalent to the root mean squared error for continuous outcomes for probabilistic forecasts.
```matlab
brierScore.heter = arrayfun(@(t) mean((t.data.heter_p(t.data.status >= 0) - t.data.status(t.data.status >= 0)).^2), panel);
brierScore.homo = arrayfun(@(t) mean((t.data.homo_p(t.data.status >= 0) - t.data.status(t.data.status >= 0)).^2), panel);
```
