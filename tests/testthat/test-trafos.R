#
#
#

set.seed(123)

sample <- matrix(evd::rfrechet(2000), 1000, 2)
sample[,1] <- 2 * sample[,1] - 5
t1 <- transform_unitfrechet(sample)

t1_retrafo <- transform_orig_margins(t1, sample)


test_that(
          "Testing trafos basic properties",
          {
            expect_gte(t1[1,1], 0)
})
