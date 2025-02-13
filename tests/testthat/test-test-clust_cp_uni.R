test_that("clust_cp_uni works", {

  data_mat <- matrix(NA, nrow = 5, ncol = 50)

  data_mat[1,] <- as.numeric(c(rnorm(30,0,0.100), rnorm(20,1,0.250)))
  data_mat[2,] <- as.numeric(c(rnorm(30,0,0.125), rnorm(20,1,0.225)))
  data_mat[3,] <- as.numeric(c(rnorm(30,0,0.175), rnorm(20,1,0.280)))
  data_mat[4,] <- as.numeric(c(rnorm(10,0,0.135), rnorm(40,1,0.225)))
  data_mat[5,] <- as.numeric(c(rnorm(10,0,0.155), rnorm(40,1,0.280)))

  out_test <- clust_cp_uni(data = data_mat,
                               n_iterations = 100,
                               B = 100, L = 1, gamma = 0.1, print_progress = FALSE)

  if(length(table(get_clust_VI(out_test$clust))) <= 5 & length(table(get_clust_VI(out_test$clust))) >= 1){
    check = TRUE
  } else {
    check = FALSE
  }


  expect_equal(check, TRUE)
})
