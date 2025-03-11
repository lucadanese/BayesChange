test_that("detect_cp works", {


  data_test <- as.numeric(c(rnorm(50,0,0.1), rnorm(50,1,0.25)))

  out_test <- detect_cp(data = data_test,
                        n_iterations = 100,
                        params = list(q = 0.25,
                                      phi = 0.1, a = 1, b = 1,
                                      c = 0.1),
                        print_progress = FALSE)

  est <- posterior_estimate(out_test, maxNClusters = 100)

  if((length(table(est))) <= 100 & (length(table(est))) >= 1){
    check = TRUE
  } else {
    check = FALSE
  }

  expect_equal(check, TRUE)
})
