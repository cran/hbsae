
context("Posterior (im)propriety")

# conditions for posterior propriety for area-level model
# see Rao (2003) p. 238:
#   under p(beta) ~ 1 and p(sv2) ~ 1 (the default),
#     or more generally p(sv2) improper but bounded,
#     the posterior for sv2 is proper provided m > p + 2
#   note that the posterior mean for sv2 may only exist if m > p + 4
test_that("numerical integration fails when it should", {
  m <- 4
  expect_error(fSAE.Area(est.init=1:m, var.init=rep(1,m), silent=TRUE))
  m <- 5
  expect_error(fSAE.Area(est.init=1:m, var.init=rep(1,m), silent=TRUE))
  m <- 6
  expect_is(fSAE.Area(est.init=1:m, var.init=rep(1,m), silent=TRUE), "sae")
  m <- 8
  expect_error(fSAE.Area(est.init=1:m, var.init=rep(1,m), prior=function(x) 1/x, silent=TRUE))
  # for some reason this last test fails for m > 8, even though the posterior for sv2 (lambda) is improper? 
})
