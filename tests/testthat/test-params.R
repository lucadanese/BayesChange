test_that("constructors return the documented defaults", {
  p <- params_uni()
  expect_s3_class(p, "BayesChangeParams")
  expect_equal(unclass(p)[c("a", "b", "c", "phi", "prior_var_phi", "prior_delta_c", "prior_delta_d")],
               list(a = 1, b = 1, c = 0.1, phi = 0.1, prior_var_phi = 0.1, prior_delta_c = 1, prior_delta_d = 1))
  expect_equal(attr(params_uni(a = 2), "defaults"), setdiff(names(p), "a"))
  e <- params_epi()
  expect_equal(unclass(e)[c("M", "xi", "a0", "b0", "I0_var", "avg_blk")],
               list(M = 250, xi = 1/8, a0 = 4, b0 = 10, I0_var = 0.01, avg_blk = 2))
  expect_null(params_multi()$S_0)
})

test_that("constructors validate their arguments", {
  expect_error(params_uni(a = -1), "a must be a positive number")
  expect_error(params_uni(phi = 1), "phi must be in")
  expect_error(params_multi(S_0 = matrix(1, 2, 3)), "square")
  expect_error(params_multi(m_0 = c(0, 0), S_0 = diag(3)), "same dimension")
  expect_error(params_epi(xi = 2), "xi must be in")
  expect_error(params_uni(alpha = 1), "unused argument")
})

test_that("print shows the values and marks defaults", {
  out <- capture.output(print(params_uni(a = 2)))
  expect_match(out[1], "univariate time series")
  expect_true(any(grepl("^  a +2$", out)))
  expect_true(any(grepl("b +1  \\(default\\)", out)))
  expect_true(any(grepl("S_0 +diag\\(0.1, d\\)", capture.output(print(params_multi())))))
})

test_that("a constructor and a plain list give identical results", {
  set.seed(1); y <- c(rnorm(50), rnorm(50, 2))
  o1 <- detect_cp(y, n_iterations = 100, params = params_uni(a = 2, prior_var_phi = 0.2),
                  kernel = "ts", print_progress = FALSE, user_seed = 1)
  o2 <- detect_cp(y, n_iterations = 100, params = list(a = 2, prior_var_phi = 0.2),
                  kernel = "ts", print_progress = FALSE, user_seed = 1)
  expect_identical(o1$orders, o2$orders)
  Y <- t(replicate(3, c(rnorm(20), rnorm(20, 3))))
  c1 <- clust_cp(Y, n_iterations = 50, B = 50, params = params_uni(phi = 0.2), kernel = "ts",
                 print_progress = FALSE, user_seed = 1)
  c2 <- clust_cp(Y, n_iterations = 50, B = 50, params = list(phi = 0.2), kernel = "ts",
                 print_progress = FALSE, user_seed = 1)
  expect_identical(c1$clust, c2$clust)
})

test_that("unknown names, wrong types and sampler settings in params are reported", {
  set.seed(1); y <- c(rnorm(50), rnorm(50, 2))
  run <- function(p) detect_cp(y, n_iterations = 20, params = p, kernel = "ts",
                               print_progress = FALSE, user_seed = 1)
  expect_warning(run(list(prior_var_ph = 0.2)), "did you mean \"prior_var_phi\"")
  expect_error(run(params_epi()), "params_uni")
  expect_warning(run(list(q = 0.3)), "ignored")
  Y <- t(replicate(3, c(rnorm(20), rnorm(20, 3))))
  runc <- function(...) clust_cp(Y, n_iterations = 30, B = 30, kernel = "ts",
                                 print_progress = FALSE, user_seed = 1, ...)
  expect_warning(c1 <- runc(params = list(q = 0.3)), "deprecated")
  c2 <- runc(params = params_uni(), q = 0.3)
  expect_identical(c1$clust, c2$clust)
})

test_that("multivariate defaults follow the data dimension", {
  r <- BayesChange:::resolve_params(params_multi(), "multi", d = 3)$values
  expect_equal(r$m_0, rep(0, 3))
  expect_equal(r$nu_0, 4)
  expect_equal(r$S_0, diag(0.1, 3, 3))
  expect_error(BayesChange:::resolve_params(params_multi(m_0 = c(0, 0)), "multi", d = 3), "length 3")
})
