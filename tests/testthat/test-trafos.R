#
#
#

sample <- matrix(rfrechet(2000), 1000, 2)
sample[,1] <- 2 * sample[,1] - 5
t1 <- transform_unitfrechet(sample)
t2 <- transform_unitfrechet(sample, empirical = F)


test_that(
          "Testing trafos basic properties",
          {
            expect_gte(t1[1,1], 0)
            expect_gte(t2[1,1], 0)
})
