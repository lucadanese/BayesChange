test_that("detect_cp_multi works", {


  data_mat <- matrix(NA, nrow = 3, ncol = 100)

  data_mat[1,] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
  data_mat[2,] <- as.numeric(c(rnorm(50,0,0.125), rnorm(50,1,0.225)))
  data_mat[3,] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))

  length_orders <- as.numeric()

  for(i in 1:5){

    out_test <- detect_cp_multi(data = data_mat,
                                       n_iterations = 100, q = 0.25,
                                       k_0 = 0.25, nu_0 = 4, phi_0 = diag(1,3,3), m_0 = rep(0,3),
                                       par_theta_c = 2, par_theta_d = 0.2, prior_var_gamma = 0.1, print_progress = FALSE)

    length_orders[i] <- length(table(get_clust_VI(out_test$order)))

  }


  if(median(length_orders) <= 100 & median(length_orders) >= 1){
    check = TRUE
  } else {
    check = FALSE
  }

  expect_equal(check, TRUE)
})
