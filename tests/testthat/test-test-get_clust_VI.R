test_that("get_clust_VI works", {

  clust <- matrix(c(0,0,0,0,0,
    0,0,0,0,1,
    1,1,1,1,1), 3,5)

  out_clust <- as.numeric(get_clust_VI(clust))

  expect_equal(out_clust, c(0,0,0,1,1))
})
