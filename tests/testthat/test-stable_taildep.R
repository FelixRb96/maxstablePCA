sample <- matrix(c(
  69.9399283, 3.960401, 69.939928,
  0.3952988, 4.416034,  4.416034,
  4.9373416,  4.937342,  1.698470,
  4.8028370, 4.802837,  2.991453,
  0.8265868,  8.180792,  8.180792,
  1.1155470, 31.144117, 31.144117
), 6, 3)

test_that(
          "Testing stable tail dependence basic  properties",
          {
            expect_lte(0, stable_tail_dependence(c(1,1,1), 3, sample))
            expect_lte(0, stable_tail_dependence(c(0,2,0.2), 2.5, sample))
            expect_lte(0, stable_tail_dependence(c(5,0,0), 5, sample))

})
