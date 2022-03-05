
context("Hyper g-prior")
# prior for between-area variance that leads to closed-form expressions,
# see
# F. Liang, R. Paulo, G. Molina, M.A. Clyde and J.O. Berger (2008).
#   Mixtures of g priors for Bayesian variable selection.
#   Journal of the American Statistical Association, 103(481), 410-423.

library(hypergeo)

# for reproducibility, even across platforms:
set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

q <- 10
nj <- round(runif(q, 2, 20))
n <- sum(nj)
df <- data.frame(
  area=rep(1:q, nj),
  y=5 + 0.25*rep(rnorm(q), nj) + rnorm(n)
)
a <- 4
prior <- function(x) (1 + x)^(-0.5*a)  # # hyper-g prior

test_that("results under hyper g-prior correspond are correct", {

  # compute R-squared
  yb <- mean(df$y)
  ybj <- as.vector(tapply(df$y, df$area, mean))
  R2 <- sum(nj * (ybj - yb)^2)
  R2 <- R2 / (R2 + sum((df$y - rep(ybj, nj))^2))
  #R2 - summary(lm(y ~ as.factor(area), df))$r.squared
  
  res <- fSAE(y ~ 1, data=df, area="area", prior=prior, w=1/nj, silent=TRUE)
  
  # posterior mean of lambda
  f1 <- hypergeo((n-1)/2, 2, (q+a-1)/2, R2)
  f2 <- hypergeo((n-1)/2, 1, (q+a-1)/2, R2)
  f3 <- 2 / (q + a - 5)
  expect_equal(res$lambdahat, as.numeric(f3 * (f1 / f2)), tolerance=1e-3)

  # small area estimates
  Xpop <- matrix(round(runif(q, 1, 10) * nj), ncol=1, dimnames=list(1:q, "(Intercept)"))
  res <- fSAE(y ~ 1, data=df, area="area", prior=prior, popdata=Xpop, w=1/nj, wpop=1/nj, silent=TRUE)
  f4 <- hypergeo((n-1)/2, 2, 1 + (q+a-1)/2, R2)
  f5 <- 2 / (q + a - 1)
  Eaj <- f5 * (f4 / f2)
  fj <- nj / as.vector(Xpop)
  expect_equal(as.numeric(EST(res)), as.numeric(fj*ybj + (1-fj)*yb + (1-fj)*Eaj*(ybj-yb)), tolerance=1e-3)

  # MSEs
  A <- (q + a - 3) / (2 * hypergeo((n-1)/2, 1, (q+a-1)/2, R2))
  EV <- A * ((n-1)/(n-3)) * var(df$y) * (1 - fj)^2
  EV <- EV * ((1/(as.vector(Xpop) - nj)) * (2/(q+a-3)) * hypergeo((n-3)/2,1,(q+a-1)/2,R2)
              + (1/nj)*(2/(q+a-1))*(2/(q+a-3))*hypergeo((n-3)/2,2,(q+a+1)/2,R2)
              + (1/n)*(2/(q+a-1))*hypergeo((n-3)/2,1,(q+a+1)/2,R2))
  VE <- A * (ybj - yb)^2 * (1 - fj)^2
  VE <- VE * (2*(2/(q+a+1))*(2/(q+a-1))*(2/(q+a-3))*hypergeo((n-1)/2,3,(q+a+3)/2,R2)
              -2*Eaj*(2/(q+a-1))*(2/(q+a-3))*hypergeo((n-1)/2,2,(q+a+1)/2,R2)
              + (Eaj)^2 * (2/(q+a-3))*hypergeo((n-1)/2,1,(q+a-1)/2,R2) )
  expect_equal(as.numeric(MSE(res)), as.numeric(EV + VE), tolerance=1e-3)
})
