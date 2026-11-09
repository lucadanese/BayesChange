# Regression tests for the issues raised in the JSS review (#6410).

params_uni <- list(prior_var_phi = 0.1, prior_delta_c = 1, prior_delta_d = 1)

make_data <- function(){
  set.seed(1)
  c(rnorm(100, 0, 1), rnorm(100, 2, 1))
}

test_that("detect_cp is reproducible given user_seed, whatever the R RNG state", {
  dat <- make_data()
  run <- function() detect_cp(dat, n_iterations = 200, n_burnin = 50,
                              params = params_uni, kernel = "ts",
                              print_progress = FALSE, user_seed = 1)
  set.seed(10); out1 <- run()
  set.seed(20); out2 <- run()
  out3 <- run()
  expect_identical(out1$orders, out2$orders)
  expect_identical(out1$orders, out3$orders)
  expect_identical(as.numeric(out1$phi_MCMC), as.numeric(out2$phi_MCMC))
})

test_that("detect_cp with user_seed leaves the caller's RNG stream untouched", {
  dat <- make_data()
  set.seed(5); x1 <- runif(3)
  set.seed(5)
  invisible(detect_cp(dat, n_iterations = 50, params = params_uni, kernel = "ts",
                      print_progress = FALSE, user_seed = 1))
  x2 <- runif(3)
  expect_identical(x1, x2)
})

test_that("detect_cp with user_seed = NULL follows set.seed()", {
  dat <- make_data()
  run <- function() detect_cp(dat, n_iterations = 100, params = params_uni, kernel = "ts",
                              print_progress = FALSE, user_seed = NULL)
  set.seed(3); out1 <- run()
  set.seed(3); out2 <- run()
  expect_identical(out1$orders, out2$orders)
})

test_that("clust_cp is reproducible given user_seed", {
  set.seed(1)
  dat <- t(replicate(3, c(rnorm(30, 0, 0.5), rnorm(30, 3, 0.5))))
  run <- function() clust_cp(dat, n_iterations = 100, n_burnin = 10, L = 1, B = 100,
                             params = list(a = 1, b = 1, c = 1, phi = 0.1),
                             kernel = "ts", print_progress = FALSE, user_seed = 2)
  set.seed(10); out1 <- run()
  set.seed(20); out2 <- run()
  expect_identical(out1$clust, out2$clust)
  expect_identical(out1$orders, out2$orders)
})

test_that("the Metropolis-Hastings update of phi rejects some proposals", {
  # Before the fix, the log full conditional of phi was -Inf (integer
  # truncation in the multivariate log-gamma), every proposal was accepted
  # and phi performed a pure random walk.
  data("eu_inflation", package = "BayesChange", envir = environment())
  out <- detect_cp(eu_inflation[1, ], n_iterations = 500, n_burnin = 100,
                   params = params_uni, kernel = "ts",
                   print_progress = FALSE, user_seed = 1)
  phi <- as.numeric(out$phi_MCMC)
  expect_true(all(is.finite(phi)))
  expect_gt(sum(diff(phi) == 0), 0)
})

test_that("posterior_estimate returns contiguous blocks and consistent change points", {
  dat <- make_data()
  out <- detect_cp(dat, n_iterations = 300, n_burnin = 100, params = params_uni,
                   kernel = "ts", print_progress = FALSE, user_seed = 1)
  est <- posterior_estimate(out, loss = "binder")
  cp <- posterior_estimate(out, loss = "binder", show_cp = TRUE)
  expect_length(est, length(dat))
  expect_true(all(diff(est) %in% c(0, 1)))
  expect_equal(est[1], 1)
  expect_identical(as.integer(cp), which(diff(est) == 1L) + 1L)
  expect_error(posterior_estimate(out, loss = "not_a_loss"))
})

test_that("computational time is recorded when print_progress = FALSE", {
  dat <- make_data()
  out <- detect_cp(dat, n_iterations = 3000, params = params_uni, kernel = "ts",
                   print_progress = FALSE, user_seed = 1)
  expect_gt(out$time, 0)
})

test_that("univariate integrated likelihood equals the closed-form marginal", {
  # Stationary OU/AR(1) block with Var = s2, Corr = phi^|i-j|,
  # s2 ~ IG(a, b), mu | s2 ~ N(0, s2 / c)  =>  x ~ multivariate t.
  exact <- function(x, phi, a, b, c){
    n <- length(x); R <- phi^abs(outer(1:n, 1:n, "-")) + matrix(1 / c, n, n)
    q <- drop(t(x) %*% solve(R, x))
    lgamma(a + n / 2) - lgamma(a) - (n / 2) * log(2 * pi * b) -
      0.5 * as.numeric(determinant(R)$modulus) - (a + n / 2) * log(1 + q / (2 * b))
  }
  set.seed(1)
  for(phi in c(0.1, 0.8)) for(n in c(1, 2, 3, 20)){
    x <- rnorm(n)
    expect_equal(BayesChange:::log_marginal_uni_ts(x, rep(0, n), phi, 2, 0.5, 0.7),
                 exact(x, phi, 2, 0.5, 0.7), tolerance = 1e-8)
  }
  x <- rnorm(30); o <- rep(0:2, c(2, 13, 15))
  expect_equal(BayesChange:::log_marginal_uni_ts(x, o, 0.3, 1, 1, 1),
               sum(sapply(split(x, o), exact, phi = 0.3, a = 1, b = 1, c = 1)), tolerance = 1e-8)
})
