#' Test for the max matrix multiplication 

A <- matrix(c(1,1,0,1,0,1,0,1,1), 3, 3)
B <- matrix(c(0,1,1,1,0,1,1,1,0), 3, 3)
res <- matrix(1, 3, 3)

H <- matrix(c(1,0,0.5,0,1,0.5,0,0,0), 3, 3)

C <- matrix(c(1,2,3,4,5,6,7,8), 4, 2)
D <- matrix(c(6,2), 2, 1)
res2 <- matrix(c(10, 12, 18, 24), 4, 1)

MMM <- matrix( c(2,1,1,2), 2, 2)
v <- c(1.5, 1.5)

test_that("Testing max-matrix multiplication", {
  expect_equal(maxmatmul(A, B), res)
  expect_equal(maxmatmul(H, H), H)
  expect_equal(maxmatmul(C, D), res2)
  # expect_equal(maxmatmul(M, v)[1], res3[1])
  # expect_equal(maxmatmul(v, M)[2], res3[2])
  # expect_equal(maxmatmul(v, v), 2.25)
})


