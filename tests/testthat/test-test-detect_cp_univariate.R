test_that("DetectCPsUnivariate works", {


  data_test <- as.numeric(c(rnorm(50,0,0.1), rnorm(50,1,0.25)))

  length_orders <- as.numeric()

  for(i in 1:10){

    out_test <- detect_cp_univariate(data = data_test,
                                     n_iterations = 2500,
                                     q = 0.25,
                                     phi = 0.1, a = 1, b = 1, c = 0.1)

    length_orders[i] <- length(table(get_order_VI(out_test$order[1000:1500,])))

  }


  if(mean(length_orders) < 3 & mean(length_orders) > 1.5){
    check = TRUE
  } else {
    check = FALSE
  }

  expect_equal(check, TRUE)
})
