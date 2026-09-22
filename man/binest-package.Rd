\name{binest-package}
\alias{binest-package}
\alias{binest}
\docType{package}

\title{Estimation of Group Means and SDs from Binned Count Data}

\description{
  Estimates group-level means and standard deviations from binned
  (coarsened) count data, where the within-bin scores are unobserved.
  All three functions fit the same heteroskedastic ordered probit
  (HETOP) model of Reardon, Shear, Castellano and Ho (2017), in which
  each group's values are normally distributed around a mean and SD of
  its own.  They differ only in how they fit it, and they share a
  common output structure:

  \itemize{
    \item \code{\link[=fast_hetop]{fast_hetop()}}: fits each group
      separately, solving closed-form truncated-normal score
      equations.  Runs in time linear in the number of groups times
      the number of bins.  This is the preferred function in the
      package.
    \item \code{\link[=mle_hetop]{mle_hetop()}}: maximizes the
      likelihood over all groups at once.  Returns the same estimates
      as \code{fast_hetop(estimator = "ML")}, far more slowly.
      Deprecated in favor of \code{\link[=fast_hetop]{fast_hetop()}}.
    \item \code{\link[=fh_hetop]{fh_hetop()}}: fits the model by MCMC,
      placing a Fay-Herriot hyperprior over the group parameters and
      reporting posterior means (Lockwood, Castellano and Shear 2018).
      Deprecated in favor of \code{\link[=fast_hetop]{fast_hetop()}}.
  }

  The \code{mle_hetop()} and \code{fh_hetop()} functions are forked from
  the HETOP package by J. R. Lockwood (CRAN, last released 2019).
  \code{mle_hetop()} has been modified to speed up its runtime via a
  vectorized inner loop and to remove two user-facing arguments
  (\code{fixedcuts} and \code{svals}) that some users found
  confusing; cutpoints and starting values are now derived internally
  from the data.  \code{fh_hetop()} has been modified the same way,
  removing its \code{fixedcuts} argument for the same reason: it was
  never a channel for supplying cutpoints known in advance on the
  test-score scale, only a coordinate-system normalization, and
  supplying native-scale cutpoints could send the MCMC sampler into a
  degenerate region rather than producing an error.

  \code{mle_hetop()} and \code{fh_hetop()} are deprecated in favor of
  \code{fast_hetop()} and remain in the package only for comparison
  purposes.  See \code{vignette("binest")} for an empirical
  comparison on Texas STAAR Grade 6 mathematics data.
}

\section{Bundled data}{
  The package ships with \code{\link{tx_g6_math_2018}}, a
  district-level dataset of bin counts and reported mean scores
  from the 2017-18 administration of the State of Texas Assessments
  of Academic Readiness (STAAR) Grade 6 mathematics test.  See the
  vignette for usage.
}

\references{
  Fisher, R. A. (1922). On the mathematical foundations of
  theoretical statistics. \emph{Philosophical Transactions of the
  Royal Society of London A}, 222, 309-368.

  Lockwood, J. R., Castellano, K. E., and Shear, B. R. (2018).
  Flexible Bayesian models for inferences from coarsened,
  group-level achievement data. \emph{Journal of Educational and
  Behavioral Statistics}, 43(6), 663-692.

  Reardon, S. F., Shear, B. R., Castellano, K. E., and Ho, A. D.
  (2017). Using heteroskedastic ordered probit models to recover
  moments of continuous test score distributions from coarsened
  data. \emph{Journal of Educational and Behavioral Statistics},
  42(1), 3-45.

  Sheppard, W. F. (1898). On the calculation of the most probable
  values of frequency-constants for data arranged according to
  equidistant divisions of a scale. \emph{Proceedings of the London
  Mathematical Society}, 29, 353-380.
}

\author{
  Paul T. von Hippel \email{ph3828@eid.utexas.edu},
  David J. Hunter, and J. R. Lockwood.
}

\keyword{package}
