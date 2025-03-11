
# BayesChange

<!-- badges: start -->
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

* `detect_cp` detect change points on univariate and multivariate time series. 
* `clust_cp` cluster time series or survival functions with common change points. 

Additional methods and functions are included: 

* `print()` and `summary()` give informations about the algorithm.
* `posterior_estimate()` estimates the change points or the final partition of the data. 
* `plot()` provides a graphical representation of the results. 
* `sim_epi_data` generates an arbitrary number of simulated survival functions. 

## Detect change points 


``` r
library(BayesChange)

## Univariate time series

data_vec <- as.numeric(c(rnorm(50,0,0.1), rnorm(50,1,0.25)))


out <- detect_cp(data = data_vec, n_iterations = 2500, n_burnin = 500,
                 params = list(q = 0.25, phi = 0.1, a = 1, b = 1, c = 0.1))

print(out)
summary(out)
posterior_estimate(out)
plot(out)

``` 

``` r

## Multivariate time series

data_mat <- matrix(NA, nrow = 3, ncol = 100)

data_mat[1,] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
data_mat[2,] <- as.numeric(c(rnorm(50,0,0.125), rnorm(50,1,0.225)))
data_mat[3,] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))


out <- detect_cp(data = data_mat, n_iterations = 2500, n_burnin = 500,
                 params = list(q = 0.25, k_0 = 0.25, nu_0 = 4, phi_0 = diag(1,3,3), m_0 = rep(0,3),
                               par_theta_c = 2, par_theta_d = 0.2, prior_var_gamma = 0.1))

print(out)
summary(out)
posterior_estimate(out)
plot(out)

``` 

## Cluster time dependent data with common change points 

``` r

## Univariate time series

data_mat <- matrix(NA, nrow = 5, ncol = 100)

data_mat[1,] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
data_mat[2,] <- as.numeric(c(rnorm(50,0,0.125), rnorm(50,1,0.225)))
data_mat[3,] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
data_mat[4,] <- as.numeric(c(rnorm(25,0,0.135), rnorm(75,1,0.225)))
data_mat[5,] <- as.numeric(c(rnorm(25,0,0.155), rnorm(75,1,0.280)))

out <- clust_cp(data = data_mat, n_iterations = 5000, n_burnin = 1000,
                params = list(L = 1, B = 1000, gamma = 0.5), kernel = "ts")

print(out)
summary(out)
posterior_estimate(out)
plot(out)

``` 

``` r

## Multivariate time series


data_array <- array(data = NA, dim = c(3,100,5))

data_array[1,,1] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
data_array[2,,1] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
data_array[3,,1] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))

data_array[1,,2] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
data_array[2,,2] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
data_array[3,,2] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))

data_array[1,,3] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
data_array[2,,3] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
data_array[3,,3] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))

data_array[1,,4] <- as.numeric(c(rnorm(25,0,0.135), rnorm(75,1,0.225)))
data_array[2,,4] <- as.numeric(c(rnorm(25,0,0.135), rnorm(75,1,0.225)))
data_array[3,,4] <- as.numeric(c(rnorm(25,0,0.135), rnorm(75,1,0.225)))

data_array[1,,5] <- as.numeric(c(rnorm(25,0,0.155), rnorm(75,1,0.280)))
data_array[2,,5] <- as.numeric(c(rnorm(25,0,0.155), rnorm(75,1,0.280)))
data_array[3,,5] <- as.numeric(c(rnorm(25,0,0.155), rnorm(75,1,0.280)))

out <- clust_cp(data = data_array, n_iterations = 3000, n_burnin = 1000,
                params = list(B = 1000, L = 1, gamma = 0.5, k_0 = 0.25,
                              nu_0 = 5, phi_0 = diag(0.1,3,3),
                              m_0 = rep(0,3)), kernel = "ts")

print(out)
summary(out)
posterior_estimate(out)
plot(out)


``` 

``` r

## Epidemiological data

data_mat <- matrix(NA, nrow = 5, ncol = 50)

betas <- list(c(rep(0.45, 25),rep(0.14,25)),
              c(rep(0.55, 25),rep(0.11,25)),
              c(rep(0.50, 25),rep(0.12,25)),
              c(rep(0.52, 10),rep(0.15,40)),
              c(rep(0.53, 10),rep(0.13,40)))

inf_times <- list()

for(i in 1:5){

  inf_times[[i]] <- sim_epi_data(10000, 10, 50, betas[[i]], 1/8)

  vec <- rep(0,50)
  names(vec) <- as.character(1:50)

  for(j in 1:50){
    if(as.character(j) %in% names(table(floor(inf_times[[i]])))){
      vec[j] = table(floor(inf_times[[i]]))[which(names(table(floor(inf_times[[i]]))) == j)]
    }
  }
  data_mat[i,] <- vec
}

out <- clust_cp(data = data_mat, n_iterations = 3000, n_burnin = 1000,
                params = list(M = 250, L = 1, B = 1000), kernel = "epi")

print(out)print(out)
summary(out)
posterior_estimate(out)
plot(out)

```

