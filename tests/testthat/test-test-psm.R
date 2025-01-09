test_that("psm works", {

  clust <- matrix(c(0,0,0,0,0,
                    0,0,0,0,1,
                    1,1,1,1,1), 3,5)

  psm_clust <- psm(clust)

  expect_equal(psm_clust, matrix(c(1,1,1,0,0,
                                   1,1,1,0,0,
                                   1,1,1,0,0,
                                   0,0,0,1,1,
                                   0,0,0,1,1), 5,5))
})
