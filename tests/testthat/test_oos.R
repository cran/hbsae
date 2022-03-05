context("Out-of-sample areas")

dat <- data.frame(area=rep(c(0:4, 6:9, 11), each=2), y=1:20)
pop <- matrix(5, nrow=9, ncol=1, dimnames=list(c(1:3, 5:10), "(Intercept)"))

test_that("oos areas are handled correctly in case of mismatch between sampled and prediction areas", {
  res <- fSAE(y ~ 1, dat, area="area", popdata=pop, silent=TRUE)
  expect_identical(res$predAreaNames[res$inSampleAreas], intersect(res$sampledAreaNames, res$predAreaNames))
  expect_identical(res$sampledAreaNames[res$predSampledAreas], intersect(res$sampledAreaNames, res$predAreaNames))
  expect_identical(res$predAreaNames, rownames(pop))
  expect_identical(names(EST(res)), res$predAreaNames)
  expect_identical(names(raneff(res)), res$predAreaNames)
  expect_identical(names(raneff(res, pop=FALSE)), res$sampledAreaNames)
  expect_identical(names(wDirect(res)), res$predAreaNames)
  expect_identical(names(wDirect(res, pop=FALSE)), res$sampledAreaNames)
})

test_that("oos areas are handled correctly when sampled and prediction areas have no overlap", {
  pop5 <- pop["5", , drop=FALSE]
  res <- fSAE(y ~ 1, dat, area="area", popdata=pop5, silent=TRUE)
  expect_identical(length(res$predAreaNames[res$inSampleAreas]), 0L)
  expect_identical(length(res$sampledAreaNames[res$predSampledAreas]), 0L)
  expect_identical(res$predAreaNames, rownames(pop5))
  expect_identical(names(EST(res)), res$predAreaNames)
  expect_identical(names(raneff(res)), res$predAreaNames)
  expect_identical(names(raneff(res, pop=FALSE)), res$sampledAreaNames)
  expect_identical(names(wDirect(res)), res$predAreaNames)
  expect_identical(names(wDirect(res, pop=FALSE)), res$sampledAreaNames)
})

