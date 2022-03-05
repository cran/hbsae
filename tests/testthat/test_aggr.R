
context("Estimation of population aggregates")

test_that("methods synthetic and survreg yield same overall population estimates", {
  set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")
  dat <- generateFakeData()
  res1 <- fSAE(y0 ~ x + area2, dat$sam, popdata=colSums(dat$Xpop), method="synthetic", silent=TRUE)
  res2 <- fSAE(y0 ~ x + area2, dat$sam, popdata=colSums(dat$Xpop), method="survreg", silent=TRUE)
  excl <- which(names(res1) %in% c("gamma", "g1", "g2", "Vbeta", "Vraneff", "method", "call"))
  expect_equal(res1[-excl], res2[-excl])
})
