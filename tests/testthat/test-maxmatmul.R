#' Test for the max matrix multiplication 

# first matrix
A <- matrix(c(1,1,0,1,0,1,0,1,1), 3, 3)
B <- matrix(c(0,1,1,1,0,1,1,1,0), 3, 3)

matprodAB <- maxmatmul(A, B)
res <- matrix(1, 3, 3)

# second matrix 
H <- matrix(c(1,0,0.5,0,1,0.5,0,0,0), 3, 3)
matprodHH <- maxmatmul(H,H)

# matrix vector
C <- matrix(c(1,2,3,4,5,6,7,8), 4, 2)
D <- matrix(c(6,2), 2, 1)

matprodCD <- maxmatmul(C,D)
res2 <- matrix(c(10, 12, 18, 24), 4, 1)

# vector matrix
v <- c(1.5, 1.5)
MMM <- matrix( c(2,1,1,2), 2, 2)

matprodvM <- maxmatmul(v, MMM)
res3 <- matrix(c(3,3), 1, 2)

# vector vector
vec1 <- 1:10
vec2 <- seq(1,0, -0.1)

matprodv1v2 <- maxmatmul(vec1, vec2)
res4 <- matrix(3, 1, 1)


matprodv1v2 <- maxmatmul(vec1, vec2)

test_that("Testing max-matrix multiplication", {

  # dimensions 
  expect_equal(dim(matprodAB), c(3,3))
  expect_equal(dim(matprodHH), c(3,3))
  expect_equal(dim(matprodCD), c(4,1))
  expect_equal(dim(matprodvM), c(1,2))
  expect_equal(dim(matprodv1v2), c(1,1))

  # check for correct results
  expect_equal(matprodAB, res)
  expect_equal(matprodHH, H)
  expect_equal(matprodCD, res2)
  expect_equal(matprodvM, res3)
  expect_equal(matprodv1v2, res4)
})


