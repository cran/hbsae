context("contents of sae-object")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")
d <- generateFakeData()
x <- fSAE(y0 ~ x + x.area + area2, data=d$sam, area="area", popdata=d$Xpop, silent=TRUE)

test_that("output contains required components", {
  expect_true(is.vector(EST(x)))
  expect_true(is.vector(MSE(x)))
  expect_equal(names(EST(x)), names(MSE(x)))
})

test_that("residuals and fitted values are vectors of right length", {
  expect_equal(length(residuals(x)), nrow(d$sam))
  expect_equal(length(fitted(x)), nrow(d$sam))
})
