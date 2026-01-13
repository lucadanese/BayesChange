test_that("clust_cp works", {

  data("stock_multi")

  params_multi <- list(m_0 = rep(0,2),
                       k_0 = 1,
                       nu_0 = 10,
                       S_0 = diag(1,2,2),
                       phi = 0.1)

  out_test <- clust_cp(data = stock_multi[,,1:5], n_iterations = 7500, n_burnin = 2500,
                  L = 1, B = 10000, params = params_multi, kernel = "ts")


  est <- posterior_estimate(out_test, maxNClusters = 3)

  if(length(table(est)) <= 5 & length(table(est)) >= 1){
    check = TRUE
  } else {
    check = FALSE
  }

  expect_equal(check, TRUE)
})


