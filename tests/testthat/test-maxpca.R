#
#
#
library(evd)
tol <- 1
set.seed <- 123456

set.seed <- 12345
A <- matrix(c(1,0,0.5, 0,1, 0.5), 3, 2)
z <- matrix(rfrechet(20000), 2, 10000)
sample <- t(maxmatmul(A, z))
maxstabPCA1 <- max_stable_prcomp(sample, 1)
maxstabPCA2 <- max_stable_prcomp(sample, 2)
maxstabPCA3 <- max_stable_prcomp(sample, 3)

maxstabPCA <- maxstabPCA2

zz <- matrix(rfrechet(200), 100, 2)
compr <- compress(maxstabPCA, zz)
reconstr <- reconstruct(maxstabPCA, compr)

zv <- rfrechet(2)
sampzv <- t(maxmatmul(A, zv))
print(sampzv)
compv <- compress(maxstabPCA, sampzv)
recv <- reconstruct(maxstabPCA, compv)



# TODO: make useful tests
test_that("Testing max-PCA and setup functions", {
  expect_equal(dim(reconstr)[2], 3)
  expect_equal(dim(reconstr)[1], 100)
  expect_equal(dim(compr)[2], 2)

  expect_equal(dim(recv)[2], 3)
  expect_equal(dim(recv)[1], 1)
  expect_equal(dim(compv)[2], 2)


})


