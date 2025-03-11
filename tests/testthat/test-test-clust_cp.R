test_that("clust_cp works", {

  data_mat <- matrix(NA, nrow = 5, ncol = 50)

  data_mat[1,] <- as.numeric(c(rnorm(30,0,0.100), rnorm(20,1,0.250)))
  data_mat[2,] <- as.numeric(c(rnorm(30,0,0.125), rnorm(20,1,0.225)))
  data_mat[3,] <- as.numeric(c(rnorm(30,0,0.175), rnorm(20,1,0.280)))
  data_mat[4,] <- as.numeric(c(rnorm(10,0,0.135), rnorm(40,1,0.225)))
  data_mat[5,] <- as.numeric(c(rnorm(10,0,0.155), rnorm(40,1,0.280)))

  out_test <- clust_cp(data = data_mat,
                       n_iterations = 100, params =  list(B = 100, L = 1, gamma = 0.1),
                       kernel = "ts", print_progress = FALSE)

  est <- posterior_estimate(out_test, maxNClusters = 3)

  if(length(table(est)) <= 5 & length(table(est)) >= 1){
    check = TRUE
  } else {
    check = FALSE
  }


  expect_equal(check, TRUE)
})

test_that("clust_cp epi works", {

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

  out_test <-  clust_cp(data = data_mat, n_iterations = 100, params = list(M = 50, B = 50, L = 1),
                        kernel = "epi", print_progress = FALSE)

  est <- posterior_estimate(out_test, maxNClusters = 3)

  if((length(table(est))) <= 3 & (length(table(est))) >= 1){
    check = TRUE
  } else {
    check = FALSE
  }

  expect_equal(check, TRUE)
})
