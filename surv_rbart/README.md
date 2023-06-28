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
