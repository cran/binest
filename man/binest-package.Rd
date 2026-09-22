\name{binest-package}
\alias{binest-package}
\alias{binest}
\docType{package}

\title{Estimation of Group Means and SDs from Binned Count Data}

\description{
  The three main package functions fit the same HETOP model, but
  differ dramatically in speed:

  \itemize{
    \item \code{\link[=fast_hetop]{fast_hetop()}} is the fastest. It
      can fit 1,000 schools or districts in less than 1 second.
      \code{fast_hetop(estimator = "ML")} produces maximum likelihood
      estimates and \code{fast_hetop(estimator = "EB_shrunk")}
      produces shrunken empirical Bayes estimates.
    \item \code{\link[=mle_hetop]{mle_hetop()}} returns the same
      estimates as \code{fast_hetop(estimator = "ML")} but far more
      slowly. When given 1,000 districts, \code{mle_hetop()} fails to
      converge after an hour. Because of the speed difference,
      \code{mle_hetop()} is deprecated in favor of
      \code{fast_hetop(estimator = "ML")}.
    \item \code{\link[=fh_hetop]{fh_hetop()}} returns estimates
      similar to \code{fast_hetop(estimator = "EB_shrunk")}, but far
      more slowly, taking over 40 minutes to estimate 1,000
      districts. \code{fh_hetop()} fits the HETOP model using Markov
      Chain Monte Carlo, placing a Fay-Herriot prior over the group
      parameters and reporting Bayesian posterior means and SDs
      (Lockwood, Castellano and Shear 2018). Because
      \code{fast_hetop(estimator = "EB_shrunk")} can produce similar
      estimates more quickly, we have deprecated \code{fh_hetop()}.
  }

  The \code{mle_hetop()} and \code{fh_hetop()} functions are forked
  from the \code{HETOP} package by J. R. Lockwood, which was last
  updated in 2019 and archived on 2025-03-24 because email to the
  maintainer was undeliverable. When forking \code{mle_hetop()} and
  \code{fh_hetop()}, we removed two arguments (\code{fixedcuts} and
  \code{svals}) that some users found confusing. \code{mle_hetop()}
  and \code{fh_hetop()} are included in \pkg{binest} for comparison
  purposes (see \code{vignette("binest")}), but they are deprecated
  and will not be maintained.
}

\section{Bundled data}{
  The package ships with \code{\link{tx_g6_math_2018}}, which provides
  bin counts and mean scores for every Texas district that
  participated in the 2017-18 administration of the State of Texas
  Assessments of Academic Readiness (STAAR) Grade 6 mathematics test.
}

\references{
  Lockwood, J. R., Castellano, K. E., and Shear, B. R. (2018).
  Flexible Bayesian models for inferences from coarsened,
  group-level achievement data. \emph{Journal of Educational and
  Behavioral Statistics}, 43(6), 663-692.
  \doi{10.3102/1076998618795124}

  Reardon, S. F., Shear, B. R., Castellano, K. E., and Ho, A. D.
  (2017). Using heteroskedastic ordered probit models to recover
  moments of continuous test score distributions from coarsened
  data. \emph{Journal of Educational and Behavioral Statistics},
  42(1), 3-45.
  \doi{10.3102/1076998616666279}
}

\author{
  Paul T. von Hippel \email{ph3828@eid.utexas.edu},
  David J. Hunter, and J. R. Lockwood.
}

\keyword{package}
