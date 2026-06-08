\name{bin_means}
\alias{bin_means}

\title{Fast Estimation of Group Means and SDs from Binned Counts}

\description{
  Estimates \code{G} group means and standard deviations from count
  data in \code{K} ordinal categories, under the assumption that
  within each group the underlying scores are normally distributed.
  Cutpoints are either supplied by the caller or derived from the
  pooled bin proportions.

  \code{bin_means} is the preferred estimator in this package and
  supersedes \code{\link{mle_hetop}} and \code{\link{fh_hetop}}: it
  runs much faster than either HETOP variant, produces an estimate
  for every group with at least three populated bins, and on real
  data has been found to be at least as accurate as either HETOP
  variant.  The HETOP functions remain in the package for comparison
  purposes.
}

\usage{
bin_means(ngk, cutpoints = NULL, eb_shrink = FALSE,
          iterate = FALSE, tol = 1e-3, maxit = 100)
}

\arguments{
  \item{ngk}{
    Numeric matrix of dimension \code{G x K} in which column \code{k} of
    row \code{g} indicates the number of units from group \code{g}
    falling into category \code{k}.
  }

  \item{cutpoints}{
    Optional numeric vector of length \code{K-1}. If \code{NULL}
    (default), cutpoints are derived from the pooled bin proportions
    via \code{qnorm(cumsum(colSums(ngk)/sum(ngk))[1:(K-1)])} and the
    output is on a standardized scale (latent population mean 0,
    total variance 1). If supplied, the cutpoints are treated as known
    on the test-score scale, and the output \code{est_raw} is on the
    test-score scale.
  }

  \item{eb_shrink}{
    Logical. If \code{TRUE}, applies empirical-Bayes shrinkage to the
    per-group estimates. Each group's log-SD is shrunk toward the
    population-weighted mean of log-SDs with weight equal to
    \code{tau^2 / (tau^2 + s_g^2)}, where \code{tau^2} is the moment
    estimator of between-group log-SD variance and \code{s_g^2} is
    the per-group sampling variance of the log-SD MLE, derived from
    the Fisher information of the binned-normal likelihood. The
    per-group mean is shrunk analogously, with the additional
    constraint of a per-group MAP step under the actual multinomial
    bin likelihood. Default \code{FALSE}.
  }

  \item{iterate}{
    Logical. If \code{TRUE} and \code{cutpoints = NULL}, the function
    refines the data-derived cutpoints by alternating between (a)
    fitting per-group means and SDs given current cutpoints and (b)
    matching the pooled bin proportions to the implied mixture-of-
    normals CDF. Iterates until cutpoints change by less than
    \code{tol}. Not valid when cutpoints are supplied. Default
    \code{FALSE}; iteration did not meaningfully improve per-group
    estimates in our tests on real or simulated data.
  }

  \item{tol}{
    Convergence tolerance (on the standardized scale) for the
    per-group MLE iteration and, when \code{iterate = TRUE}, the
    cutpoint-refinement EM. Default \code{1e-3}.
  }

  \item{maxit}{
    Maximum number of iterations. Default \code{100}.
  }
}

\details{
  The function derives \code{K-1} cutpoints on a standardized scale
  (mean 0, SD 1) by applying \code{qnorm} to the cumulative pooled bin
  proportions \code{cumsum(colSums(ngk)/sum(ngk))[1:(K-1)]}.  For each
  bin \code{k} with cutpoints \code{c_\{k-1\}} and \code{c_k}, the function
  then computes the truncated-normal first and second moments:

    p_k       = Phi(c_k) - Phi(c_\{k-1\})
    E[Z|k]    = ( phi(c_\{k-1\}) - phi(c_k) ) / p_k
    E[Z^2|k]  = 1 + ( c_\{k-1\}*phi(c_\{k-1\}) - c_k*phi(c_k) ) / p_k

  For each group \code{g}, the estimated mean is the within-group
  bin-proportion-weighted average of \code{E[Z|k]}, and the estimated
  variance is the within-group bin-proportion-weighted average of
  \code{E[Z^2|k]} minus the squared estimated mean.

  Estimates are returned on two scales.
}

\value{
  A list with the following components:

  \item{est_raw}{Estimates on the raw (test-score) scale when
    \code{cutpoints} were supplied, or on the standardized scale
    otherwise.  A list with elements \code{group_mean_<suffix>},
    \code{group_sd_<suffix>}, \code{cutpoints}, and \code{icc}
    (intracluster correlation), where \code{<suffix>} is \code{mle}
    when \code{eb_shrink = FALSE} and \code{eb} when
    \code{eb_shrink = TRUE}.}

  \item{est_std}{Estimates on the standardized scale where the
    population-weighted state mean is 0 and the total (within +
    between) state SD is 1.  Same elements as \code{est_raw}.}

  \item{gof}{Per-group Pearson chi-square goodness-of-fit of the
    within-group normality assumption: a data frame with columns
    \code{chisq}, \code{df}, and \code{p}.  \code{NA} when the test
    is unidentified (fewer than four populated bins).}

  \item{iter_info}{Diagnostic flags and tuning parameters from the
    fit, including \code{within} and \code{eb_shrink}.}
}

\references{
  Sheppard W.F. (1898). On the calculation of the most probable values
  of frequency-constants for data arranged according to equidistant
  divisions of a scale.  \emph{Proceedings of the London Mathematical
  Society}, 29, 353-380.

  Fisher R.A. (1922). On the mathematical foundations of theoretical
  statistics.  \emph{Philosophical Transactions of the Royal Society
  A}, 222, 309-368.
}

\author{
  Paul T. von Hippel and David J. Hunter.
}

\examples{
set.seed(1001)
G   <- 10
mug <- seq(from = -2.0, to = 2.0, length = G)
sigmag <- seq(from =  2.0, to = 0.8, length = G)
cutpoints <- c(-1.0, 0.0, 0.8)
ng  <- rep(1000, G)
ngk <- gendata_hetop(G, K = 4, ng, mug, sigmag, cutpoints)

bm <- bin_means(ngk)
print(cbind(true = mug,       est = bm$est_raw$group_mean_mle))
print(cbind(true = sigmag,    est = bm$est_raw$group_sd_mle))
print(cbind(true = cutpoints, est = bm$est_raw$cutpoints))
}

\keyword{models}
