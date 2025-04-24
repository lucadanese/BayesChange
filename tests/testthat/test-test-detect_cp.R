test_that("detect_cp works", {


  data_test <- as.numeric(c(rnorm(50,0,0.1), rnorm(50,1,0.25)))

  out_test <- detect_cp(data = data_test,
                        n_iterations = 100,
                        print_progress = FALSE, kernel = "ts")

  est <- posterior_estimate(out_test, maxNClusters = 100)

  if((length(table(est))) <= 100 & (length(table(est))) >= 1){
    check = TRUE
  } else {
    check = FALSE
  }

  expect_equal(check, TRUE)
})
