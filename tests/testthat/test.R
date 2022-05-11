data('sim_PO_data')
beta = c(0.5,0,-0.5, rep(0,10))
testthat::test_that("difference between beta_hat and true beta", {
  res = PO(Surv(X, delta) ~ Z[,1]+ Z[,2]+ Z[,3]+ Z[,4]+ Z[,5]+ Z[,6]+ Z[,7]
           + Z[,8]+ Z[,9]+ Z[,10]+ Z[,11]+ Z[,12]+ Z[,13],
           data = sim_PO_data,method = 'NPMLE')
  diff_check = max(abs(res$coefficients - beta)) < 0.2
  names(diff_check) = c()
  expect_equal(diff_check, TRUE)

})

testthat::test_that("difference between beta_hat and true beta", {
  res = PO(Surv(X, delta) ~ Z[,1]+ Z[,2]+ Z[,3]+ Z[,4]+ Z[,5]+ Z[,6]+ Z[,7]
           + Z[,8]+ Z[,9]+ Z[,10]+ Z[,11]+ Z[,12]+ Z[,13],
           data = sim_PO_data,method = 'U-method')
  diff_check = max(abs(res$coefficients - beta)) < 0.2
  names(diff_check) = c()
  expect_equal(diff_check, TRUE)
})


