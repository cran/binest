\name{waic_hetop}
\alias{waic_hetop}

\title{WAIC for FH-HETOP model}

\description{
  \code{waic_hetop()} is a helper function used by
  \code{\link[=fh_hetop]{fh_hetop()}}. \code{waic_hetop()} computes the
  Watanabe-Akaike information criterion (WAIC) for the FH-HETOP model
  using the data and posterior samples of the group means, group
  standard deviations and cutpoints.
}

\usage{
waic_hetop(ngk, samps)
}

\arguments{
  \item{ngk}{
    Data giving bin counts. Column \code{k} of row \code{g} counts how
    many students from school or district \code{g} scored in bin
    \code{k}.
  }

  \item{samps}{
    A matrix of posterior samples that includes at least the
    group means, group standard deviations and the cutpoints.  Column
    names for these three collections of parameters must contain the
    strings 'mu', 'sigma' and 'cuts', respectively.
  }
}

\details{
  Although this function can be called directly by the user, it is
  primarily intended to be used to compute WAIC as part of the function
  \code{fh_hetop()}.  Details on the WAIC calculation are provided by
  Vehtari and Gelman (2017).
}

\value{
  A list with the following components:

  \item{lpd_hat}{Part 1 of the WAIC calculation: the estimated log
  pointwise predictive density, summed across groups.}

  \item{phat_waic}{Part 2 of the WAIC calculation: the effective number
  of parameters.}

  \item{waic}{The WAIC criterion: -2 times (lpd_hat - phat_waic).}
}

\references{
  Lockwood J.R., Castellano K.E. and Shear B.R. (2018).
  \dQuote{Flexible Bayesian models for inferences from coarsened,
    group-level achievement data,}
  \emph{Journal of Educational and Behavioral Statistics.} 43(6):663--692.

  Vehtari A., Gelman A. and Gabry J. (2017).
  \dQuote{Practical Bayesian model evaluation using leave-one-out
    cross-validation and WAIC,}
  \emph{Statistics and Computing.} 27(5):1413--1432.
}

\author{
  J.R. Lockwood \email{jrlockwood@ets.org}
}

\examples{
\donttest{
if (requireNamespace("R2jags", quietly = TRUE)) {
  set.seed(42)
  G <- 10
  ngk <- gendata_hetop(G = G, K = 4, ng = rep(50, G),
                       mug = rnorm(G), sigmag = exp(rnorm(G, 0, 0.2)),
                       cutpoints = c(-1, 0, 1))
  m <- fh_hetop(ngk,
                p = c(10, 10), m = c(100, 100),
                gridL = c(-5, log(0.10)), gridU = c(5, log(5.0)),
                n.iter = 200, n.burnin = 100, seed = 1)
  waic <- waic_hetop(ngk, m$BUGSoutput$sims.matrix)
  print(waic)
}
}
}

\keyword{utilities}
