testthat::test_that("str_length is number of characters", {
  data('sim_PO_data')
  res = PO(Surv(X, delta) ~ Z[,1]+ Z[,2]+ Z[,3]+ Z[,4]+ Z[,5]+ Z[,6]+ Z[,7]
           + Z[,8]+ Z[,9]+ Z[,10]+ Z[,11]+ Z[,12]+ Z[,13],
           data = sim_PO_data,method = 'NPMLE')
  expect_output(str(res), "List of 7")

  beta = c(0.5,0,-0.5, rep(0,10))
  diff = abs(res$coefficients - beta)
  diff_check = diff < 10
  ture_vector = rep(TRUE,13)
  expect_equal(diff_check, ture_vector)

})
