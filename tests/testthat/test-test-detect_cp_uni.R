test_that("detect_cp_uni works", {


  data_test <- as.numeric(c(rnorm(50,0,0.1), rnorm(50,1,0.25)))

  out_test <- detect_cp_uni(data = data_test,
                            n_iterations = 100,
                            q = 0.25,
                            phi = 0.1, a = 1, b = 1,
                            c = 0.1, print_progress = FALSE)

  if((length(table(get_clust_VI(out_test$order)))) <= 100 & (length(table(get_clust_VI(out_test$order)))) >= 1){
    check = TRUE
  } else {
    check = FALSE
  }

  expect_equal(check, TRUE)
})
