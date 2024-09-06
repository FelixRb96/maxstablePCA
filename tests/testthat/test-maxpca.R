#
#
#
library(evd)
tol <- 1
set.seed(123456)

A <- matrix(c(1,0,0.5, 0,1, 0.5), 3, 2)
z <- matrix(rfrechet(20000), 2, 10000)
sample <- t(maxmatmul(A, z))
maxstabPCA1 <- max_stable_prcomp(sample, 1)
maxstabPCA2 <- max_stable_prcomp(sample, 2)
maxstabPCA3 <- max_stable_prcomp(sample, 3)

maxstabPCA <- maxstabPCA2

zz <- matrix(rfrechet(300), 100, 3)
xx <- t(maxmatmul(A, t(zz))) 
compr <- compress(maxstabPCA, xx)
reconstr <- reconstruct(maxstabPCA, xx)

zv <- matrix(rfrechet(4), 2, 2)
sampzv <- t(maxmatmul(A, zv)) 
compv <- compress(maxstabPCA, sampzv)
recv <- reconstruct(maxstabPCA, sampzv)



test_that("Testing max-PCA and setup functions", {
  
  # dimension checks
  expect_equal(dim(reconstr)[2], 3)
  expect_equal(dim(reconstr)[1], 100)
  expect_equal(dim(compr)[2], 2)

  expect_equal(dim(recv)[2], 3)
  expect_equal(dim(recv)[1], 2)
  expect_equal(dim(compv)[2], 2)

  # check if optimizers converge 
  expect_true(maxstabPCA1$optim_conv_status > 0)
  expect_true(maxstabPCA2$optim_conv_status > 0)
  expect_true(maxstabPCA3$optim_conv_status > 0)

  # check if summary output is generated
  expect_output(summary(maxstabPCA1))
  expect_output(summary(maxstabPCA2))
  expect_output(summary(maxstabPCA3))

  # check if max_stable_prcomp fails if given negative values
  expect_error(max_stable_prcomp(matrix(-10:10, 10, 2), p = 1, s = 1))

  # check that all reconstructions and encodings are positive
  expect_true(all(compr > 0))
  expect_true(all(recv > 0))

  # check that all matrices of max_stable_prcomp are non-negative
  expect_true(all(maxstabPCA1$encoder_matrix >= 0))
  expect_true(all(maxstabPCA1$decoder_matrix >= 0))
  expect_true(all(maxstabPCA1$reconstr_matrix >= 0))
  
  expect_true(all(maxstabPCA2$encoder_matrix >= 0))
  expect_true(all(maxstabPCA2$decoder_matrix >= 0))
  expect_true(all(maxstabPCA2$reconstr_matrix >= 0))

  expect_true(all(maxstabPCA3$encoder_matrix >= 0))
  expect_true(all(maxstabPCA3$decoder_matrix >= 0))
  expect_true(all(maxstabPCA3$reconstr_matrix >= 0))

})


