
# BayesChange

<!-- badges: start -->
<!-- badges: end -->

`BayesChange` provides C++ functions to perform Bayesian Change Points analysis. 

## Installation

To install `BayesChange` the package `devtools` is needed. 

``` r
install.packages("devtools")
```

Now `BayesChange` can be installed through the [GitHub repository](https://github.com/lucadanese/BayesChange) of the package:

``` r
devtools::install_github("lucadanese/BayesChange")
```

## Example

The package contains five main functions: 

* `detect_cp_uni` and `detect_cp_multi` perform change points detection on univariate and multivariate time series data.

* `clust_cp_uni` and `clust_cp_multi` cluster univariate and multivariate time series with common change points. 
* `clust_cp_epi` perform clustering of survival functions with common change points. 

Additional functions are included: 

* `sim_epi_data` generates an arbitrary number of simulated survival functions. 
* `get_clust_VI` returns a point estimate for the data. 


``` r
library(BayesChange)

# detect_cp_uni

data_vec <- as.numeric(c(rnorm(50,0,0.1), rnorm(50,1,0.25)))

out <- detect_cp_uni(data = data_vec,                             
                            n_iterations = 2500,
                            q = 0.25,
                            phi = 0.1, a = 1, b = 1, c = 0.1)

get_clust_VI(out$order)

``` 

``` r

# detect_cp_multi

data_mat <- matrix(NA, nrow = 3, ncol = 100)

data_mat[1,] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
data_mat[2,] <- as.numeric(c(rnorm(50,0,0.125), rnorm(50,1,0.225)))
data_mat[3,] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))

out <- detect_cp_multi(data = data_mat,
                              n_iterations = 2500,
                              q = 0.25,k_0 = 0.25, nu_0 = 4, phi_0 = diag(1,3,3), m_0 = rep(0,3),
                              par_theta_c = 2, par_theta_d = 0.2, prior_var_gamma = 0.1)

get_clust_VI(out$order)

``` 

``` r

# clust_cp_uni

data_mat <- matrix(NA, nrow = 5, ncol = 100)

data_mat[1,] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
data_mat[2,] <- as.numeric(c(rnorm(50,0,0.125), rnorm(50,1,0.225)))
data_mat[3,] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
data_mat[4,] <- as.numeric(c(rnorm(25,0,0.135), rnorm(75,1,0.225)))
data_mat[5,] <- as.numeric(c(rnorm(25,0,0.155), rnorm(75,1,0.280)))

out <- clust_cp_uni(data = data_mat, n_iterations = 5000, B = 1000, L = 1, gamma = 0.5)

get_clust_VI(out$clust[2500:5000,])

``` 

``` r

# clust_cp_multi

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

out <- clust_cp_multi(data = data_array, n_iterations = 5000, B = 1000, L = 1,
                        gamma = 0.1, k_0 = 0.25, nu_0 = 5, phi_0 = diag(0.1,3,3), m_0 = rep(0,3))

get_clust_VI(out$clust[2500:5000,])

``` 

``` r

# clust_cp_epi

data_mat <- matrix(NA, nrow = 5, ncol = 50)

betas <- list(c(rep(0.45, 25),rep(0.14,25)),
               c(rep(0.55, 25),rep(0.11,25)),
               c(rep(0.50, 25),rep(0.12,25)),
               c(rep(0.52, 10),rep(0.15,40)),
               c(rep(0.53, 10),rep(0.13,40)))

  inf_times <- list()

  for(i in 1:5){

    inf_times[[i]] <- sim_epi_data(S0 = 10000, I0 = 10, max_time = 50, beta_vec = betas[[i]], gamma_0 = 1/8)

    vec <- rep(0,50)
    names(vec) <- as.character(1:50)

    for(j in 1:50){
      if(as.character(j) %in% names(table(floor(inf_times[[i]])))){
      vec[j] = table(floor(inf_times[[i]]))[which(names(table(floor(inf_times[[i]]))) == j)]
      }
    }
    data_mat[i,] <- vec
  }

  out <- clust_cp_epi(data = data_mat, n_iterations = 5000, M = 500, B = 1000, L = 1)

get_clust_VI(out$clust[1000:5000,])



```

