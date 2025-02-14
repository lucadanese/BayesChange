test_that("clust_cp_multi works", {

  data_array <- array(data = NA, dim = c(2,50,3))

  data_array[1,,1] <- as.numeric(c(rnorm(20,0,0.100), rnorm(30,1,0.250)))
  data_array[2,,1] <- as.numeric(c(rnorm(20,0,0.100), rnorm(30,1,0.250)))

  data_array[1,,2] <- as.numeric(c(rnorm(20,0,0.100), rnorm(30,1,0.250)))
  data_array[2,,2] <- as.numeric(c(rnorm(20,0,0.100), rnorm(30,1,0.250)))

  data_array[1,,3] <- as.numeric(c(rnorm(10,0,0.155), rnorm(40,1,0.280)))
  data_array[2,,3] <- as.numeric(c(rnorm(10,0,0.155), rnorm(40,1,0.280)))

  out_test <- clust_cp_multi(data = data_array, n_iterations = 100, B = 100, L = 1,
                             gamma = 0.1, k_0 = 0.25, nu_0 = 5, phi_0 = as.matrix(diag(0.1,2,2)), m_0 = rep(0,2), print_progress = FALSE)

  if((length(table(get_clust_VI(out_test$clust))) <= 3) & (length(table(get_clust_VI(out_test$clust))) >= 1)){
    check = TRUE
  } else {
    check = FALSE
  }

  expect_equal(check, TRUE)
})
