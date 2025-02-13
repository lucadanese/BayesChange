test_that("clust_cp_epi works", {

  length_orders <- as.numeric()

  data_mat <- matrix(NA, nrow = 3, ncol = 20)

  betas <- list(c(rep(0.45, 10),rep(0.14,10)),
                c(rep(0.52, 5),rep(0.15,15)),
                c(rep(0.53, 5),rep(0.13,15)))

  inf_times <- list()

  for(i in 1:3){

    inf_times[[i]] <- sim_epi_data(S0 = 10000, I0 = 10, max_time = 20, beta_vec = betas[[i]], gamma_0 = 1/8)

    vec <- rep(0,20)
    names(vec) <- as.character(1:20)

    for(j in 1:20){
      if(as.character(j) %in% names(table(floor(inf_times[[i]])))){
        vec[j] = table(floor(inf_times[[i]]))[which(names(table(floor(inf_times[[i]]))) == j)]
      }
    }
    data_mat[i,] <- vec
  }

  out_test <-  clust_cp_epi(data = data_mat, n_iterations = 100, M = 50, B = 50, L = 1, print_progress = FALSE)

  if((length(table(get_clust_VI(out_test$clust)))) <= 3 & (length(table(get_clust_VI(out_test$clust)))) >= 1){
    check = TRUE
  } else {
    check = FALSE
  }

  expect_equal(check, TRUE)
})
