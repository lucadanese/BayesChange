
<img src="bayeschange_logo.png" width="250" align="center"/>

<!-- badges: start -->
[![CRAN Status](http://www.r-pkg.org/badges/version/BayesChange)](https://cran.r-project.org/package=BayesChange)
<!-- badges: end -->

`BayesChange` provides C++ functions to perform Bayesian change points analysis.

## Installation

To install `BayesChange` the package `devtools` is needed. 

``` r
install.packages("devtools")
```

Now `BayesChange` can be installed through the [GitHub repository](https://github.com/lucadanese/BayesChange) of the package:

``` r
devtools::install_github("lucadanese/BayesChange")
```

## Package contents

The package contains two main functions: 

* `detect_cp` change points detection on time series and epidemic diffusions. 
* `clust_cp` clustering of time series or epidemic diffusions with common change points. 

Additional methods and functions are included: 

* `print()` and `summary()` return information about the algorithm.
* `posterior_estimate()` estimates the change points or the final partition of the data. 
* `plot()` provides a graphical representation of the results.
* `plot_psm()` provides the posterior similarity matrix for the output of `clust_cp`.
* `sim_epi_data()` generates an arbitrary number of simulated survival functions.
  

## Detect change points 


``` r
library(BayesChange)

## Univariate time series

data("stock_uni")

params_uni <- list(a = 1,
                   b = 1,
                   c = 1,
                   phi = 0.1)

out <- clust_cp(data = stock_uni[1:5,], n_iterations = 7500, n_burnin = 2500,
                L = 1, q = 0.5, B = 10000, params = params_uni, kernel = "ts")

print(out)
summary(out)
posterior_estimate(out)
plot(out)


``` 

``` r

## Multivariate time series

data("stock_multi")

params_multi <- list(m_0 = rep(0,2),
                     k_0 = 1,
                     nu_0 = 10,
                     S_0 = diag(1,2,2),
                     phi = 0.1)

out <- clust_cp(data = stock_multi[,,1:5], n_iterations = 7500, n_burnin = 2500,
                L = 1, B = 10000, params = params_multi, kernel = "ts")

print(out)
summary(out)
posterior_estimate(out)
plot(out)

```

``` r

## Epidemic diffusions

data("epi_synthetic_multi")

params_epi <- list(M = 250, xi = 1/8,
                   alpha_SM = 1,
                   a0 = 4,
                   b0 = 10,
                   I0_var = 0.1,
                   avg_blk = 2)

out <- clust_cp(epi_synthetic_multi, n_iterations = 5000, n_burnin = 2000,
                L = 1, B = 1000, params = params_epi, kernel = "epi")

print(out)
summary(out)
posterior_estimate(out)
plot(out)

``` 

## Cluster time dependent data with common change points 

``` r

## Univariate time series

data("stock_uni")

params_uni <- list(a = 1,
                   b = 1,
                   c = 1,
                   phi = 0.1)

out <- clust_cp(data = stock_uni[1:5,], n_iterations = 7500, n_burnin = 2500,
                L = 1, q = 0.5, B = 10000, params = params_uni, kernel = "ts")

print(out)
summary(out)
posterior_estimate(out)
plot(out)

``` 

``` r

## Multivariate time series

data("stock_multi")

params_multi <- list(m_0 = rep(0,2),
                     k_0 = 1,
                     nu_0 = 10,
                     S_0 = diag(1,2,2),
                     phi = 0.1)

out <- clust_cp(data = stock_multi[,,1:5], n_iterations = 7500, n_burnin = 2500,
                L = 1, B = 10000, params = params_multi, kernel = "ts")

print(out)
summary(out)
posterior_estimate(out)
plot(out)


``` 

``` r

## Epidemic diffusions

data("epi_synthetic_multi")

params_epi <- list(M = 250, xi = 1/8,
                   alpha_SM = 1,
                   a0 = 4,
                   b0 = 10,
                   I0_var = 0.1,
                   avg_blk = 2)

out <- clust_cp(epi_synthetic_multi, n_iterations = 5000, n_burnin = 2000,
                L = 1, B = 1000, params = params_epi, kernel = "epi")

print(out)
summary(out)
posterior_estimate(out)
plot(out)

```

