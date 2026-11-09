
<img src="bayeschange_logo.png" width="250" align="center"/>

<!-- badges: start -->
[![CRAN Status](http://www.r-pkg.org/badges/version/BayesChange)](https://cran.r-project.org/package=BayesChange)
<!-- badges: end -->

`BayesChange` provides C++ functions to perform Bayesian change points analysis.

## Installation

The released version is on CRAN. On Windows and macOS CRAN provides
pre-compiled binaries, so no additional software is needed:

``` r
install.packages("BayesChange")
```

### Installing from source

Installing from source (the GitHub version or a `.tar.gz` file) requires a
C++17 compiler and the [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/).
The R packages `Rcpp`, `RcppArmadillo` and `RcppGSL` are installed
automatically as dependencies.

1. **Install GSL and a compiler.**

   | System | Command |
   |---|---|
   | macOS | `xcode-select --install` (compiler), then `brew install gsl` ([Homebrew](https://brew.sh)) |
   | Debian / Ubuntu | `sudo apt-get install build-essential libgsl-dev` |
   | Fedora / RHEL | `sudo dnf install gcc-c++ gsl-devel` |
   | Windows | Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version; it already contains GSL. |

   On macOS and Linux, check that `gsl-config --version` works in a terminal:
   the installer uses it to locate GSL.

2. **Install the package.**

   ``` r
   install.packages("remotes")
   remotes::install_github("lucadanese/BayesChange")
   # or, from a source file:
   install.packages("BayesChange_x.y.z.tar.gz", repos = NULL, type = "source")
   ```

If GSL is installed in a non-standard location, point the installer to it:

``` sh
# macOS / Linux
R CMD INSTALL --configure-vars='GSL_CFLAGS=-I/path/to/gsl/include GSL_LIBS="-L/path/to/gsl/lib -lgsl -lgslcblas"' BayesChange_x.y.z.tar.gz
```

On Windows, set the environment variable `LIB_GSL` to the GSL installation
folder (containing `include/` and `lib/`) before installing.

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

