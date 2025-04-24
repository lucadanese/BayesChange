test_that("clust_cp works", {

  data_mat <- matrix(NA, nrow = 5, ncol = 50)

  data_mat[1,] <- as.numeric(c(rnorm(30,0,0.100), rnorm(20,1,0.250)))
  data_mat[2,] <- as.numeric(c(rnorm(30,0,0.125), rnorm(20,1,0.225)))
  data_mat[3,] <- as.numeric(c(rnorm(30,0,0.175), rnorm(20,1,0.280)))
  data_mat[4,] <- as.numeric(c(rnorm(10,0,0.135), rnorm(40,1,0.225)))
  data_mat[5,] <- as.numeric(c(rnorm(10,0,0.155), rnorm(40,1,0.280)))

  out_test <- clust_cp(data = data_mat,
                       n_iterations = 100, params =  list(B = 100, L = 1, phi = 0.1),
                       kernel = "ts", print_progress = FALSE)

  est <- posterior_estimate(out_test, maxNClusters = 3)

  if(length(table(est)) <= 5 & length(table(est)) >= 1){
    check = TRUE
  } else {
    check = FALSE
  }


  expect_equal(check, TRUE)
})


