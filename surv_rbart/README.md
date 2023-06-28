# rbart Package
The package created for the Pratola, Chipman, and McCulloch (2020) is adapted here for a probit framework. Note that when using `rbart` one should be cautious. For the bankruptcy prediction example you should build this package. Make sure `devtools` are installed and attached by using
```
install.packages("devtools")
library(devtools)
```
Moreover, verify that any previous installation of the package are removed by running
```
remove.packages("rbart") 
```
Then, build the package using
```
cd <root dir>
devtools::check()
devtools::install()
```
Note that the `devtools::check()` is completely optional, but is a recommended practice. When using `rbart` for the other examples, again ensure that the previous installation is removed and install it like any other CRAN package
```
remove.packages("rbart")
install.packages("rbart")
library("rbart")
```

## Differences
To allow for dichotomous data and using the probit framework, I have altered the `rbart` slightly. Note that I have written this section as if it was produced by `git getdiff`. However, I have altered the notation for the lines (between `@@ @@`). The first pair of numbers corresponds to deletion of lines (if present) and the second to the addition of lines (if present). The first number in each pair corresponds to the number of lines added or removed and the second corresponds to the line number of the original file in case of removal and in the new file in case of addition. For example, `@@ -1,10 +1,10 @@` means that a line at linenumber 10 got removed in the original file and replaced by a new line in the new file.

### R/rbart.R
In the file that is being called in R I have made a few small changes. First of all, I have removed the line that centered the outcome data because this is nonsensical in the context of dichotomous data.
```diff
@@ -1,18 @@
- y.train = y.train - fmean
```
Secondly, I have altered the value for `tau`, to match that which was the recommended setting in the probit application of Chipman, George, McCulloch (2010). That is,

```diff
@@ -2,31 +1,31 @@
- rgy = range(y.train)
- tau = (rgy[2] - rgy[1])/(2 * sqrt(m) * k)
+ tau =  2.25/(sqrt(m)*k)
```

### src/cpsambrt.cpp
For the Gibbs sampler I need to draw from the truncated normal distribution. For this I have borrowed a function from the `bart` package written by McCulloch, Sparapani, and Gramacy. 
```diff
@@ +1,46 @@
+ #include "rtnorm.h"
```

I added the extra step for the Gibbs sampler. This consists of drawing initial values for the latent variables.
```diff
@@ +10,245 @@
+ //initial latent variable draws
+ double *iz = new double[n];
+ //loop over obs
+ for(size_t j=0;j<n;j++) {
+     if (y[j]==0) {
+     iz[j] = -rtnorm(0, 1, gen);
+     } else {
+     iz[j] = rtnorm(0, 1, gen);
+     }
+ }
```

To actually use this data in the analysis, I added `iz` to the data used in constructing the trees.
```diff
@@ -1,260 +1,260 @@
- di.n=n;di.p=p,di.x = x; di.y = iz; di.tc=tc;
+ di.n=n;di.p=p,di.x = x; di.y = y; di.tc=tc;
```

Then, the code does some adaptations to the trees, in which we start to incorporate drawing the latent variables in the Gibbs sampler.
```diff
@@ +10,354 @@
+ //loop over obs
+       for(size_t j=0;j<n;j++) {
+         //see probit transformation
+         if (y[j]==0) {
+           di.y[j] = -rtnorm(-ambm.f(j), psbm.f(j), gen); //mirror distribution
+         } else {
+           di.y[j] = rtnorm(ambm.f(j), psbm.f(j), gen);
+         }
+       }
+       if(i % printevery==0) COUT << "At iteration " << i << " of " << nadapt << " for adapt" << endl;
```

This snippet is repeated at line `372` and `394`.

### src/rtnorm.h and src/rtnorm.cpp
This is the implementation of the truncated normal distribution. It is copied from the `bart` package by McCulloch, Sparapani, and Gramacy. These files were originally not in the package.
