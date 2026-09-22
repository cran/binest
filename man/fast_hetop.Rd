\name{fast_hetop}
\alias{fast_hetop}

\title{fast_hetop(): Fast Estimation of the HETOP model for Binned Test Scores}

\description{
  \code{fast_hetop()} is the preferred function in this package. It can
  produce estimates for hundreds or thousands of schools or districts
  in less than a second. It can provide either maximum likelihood
  estimates (\code{estimator = "ML"}) or shrunken empirical Bayes
  estimates (\code{estimator = "EB_shrunk"}). It also has several
  features that the other package functions lack, including

  \itemize{
    \item a test for whether scores really follow a normal
      distribution, as the HETOP model assumes;
    \item the ability to accept published cut scores;
    \item the ability to adjust standard errors downward when the data
      represent a population rather than a sample.
  }

  \code{fast_hetop(estimator = "ML")} provides identical estimates to
  \code{\link[=mle_hetop]{mle_hetop()}}, but runs much more quickly,
  especially when the number of schools or districts is large.

  \code{fast_hetop(estimator = "EB_shrunk")} provides similar but not
  identical estimates to \code{\link[=fh_hetop]{fh_hetop()}}, and again runs much
  more quickly.

}

\usage{
fast_hetop(ngk, cutpoints_known = FALSE, cutpoints = NULL,
          pooled_mean = 0, pooled_sd = 1,
          scope,
          estimator = "ML",
          tol = 1e-4, maxit = 100,
          conf.level = 0.95,
          estimate_unidentified_districts = TRUE)
}

\arguments{
  \item{ngk}{
    Data giving bin counts. Column \code{k} of row \code{g} counts how
    many students from school or district \code{g} scored in bin
    \code{k}.
  }

  \item{cutpoints_known}{
    Whether the cut scores are known.

    \code{TRUE} means they are known and supplied through the
    \code{cutpoints} argument.

    \code{FALSE} means they are unknown and must be estimated, assuming
    that the distribution pooled across all schools or districts has
    the mean and SD given by \code{pooled_mean} and \code{pooled_sd}.
  }

  \item{cutpoints}{
    A vector of cut scores. Required when
    \code{cutpoints_known = TRUE}. Cut scores must be listed in
    ascending order, and there must be one fewer cut score than the
    number of bins.
  }

  \item{pooled_mean, pooled_sd}{
    Used only when \code{cutpoints_known = FALSE}: the mean and SD to
    give the score distribution pooled across schools or districts.
    Defaults \code{0} and \code{1} report results on a standardized
    scale.
  }

  \item{scope}{
    Whether the data represent a sample of students or the whole
    population. Required, with no default, because it determines what
    the reported standard errors mean: \code{"sample"} counts both
    sampling and binning error, while \code{"population"} counts
    binning error alone.
  }

  \item{estimator}{
    Which estimator to compute. \code{"ML"} (default) returns maximum
    likelihood estimates for each school or district.
    \code{"EB_shrunk"} returns empirical Bayes estimates, which shrink
    each school or district toward the average, by an amount that
    depends on how precisely it is estimated. Shrinkage requires
    \code{scope = "sample"}.
  }

  \item{tol}{
    Stop iterating when the mean and SD estimates change by less than
    this fraction of an SD. Default \code{1e-4}.
  }

  \item{maxit}{
    Maximum number of iterations. Default \code{100}.
  }

  \item{conf.level}{
    Confidence level for reported confidence intervals. Default
    \code{0.95}.
  }

  \item{estimate_unidentified_districts}{
    If \code{TRUE} (default), districts with fewer than three
    populated bins will receive estimates. The HETOP model cannot
    identify both the mean and the SD of a district with fewer than
    three populated bins, but the mean can be estimated by assuming a
    typical value for the SD. Specifically, the SD is assumed to be the
    geometric mean of the SDs of the other districts.
  }
}

\value{
  A list with the following components:

  \item{est_raw}{Returned only when \code{cutpoints_known = TRUE}.
    Lists each school or district's \code{mean},
    \code{sd}, \code{mean_se}, \code{sd_se},
    \code{mean_ci_lower}, \code{mean_ci_upper},
    \code{sd_ci_lower}, \code{sd_ci_upper}, along with
    \code{scope}, \code{estimator}, \code{conf.level},
    \code{cutpoints}, and \code{icc} (the intracluster correlation).}

  \item{est_std}{Estimates on the standardized scale, where the
    population-weighted state mean is 0 and the total (within plus
    between) state SD is 1. Same elements as \code{est_raw}.

    That total SD is computed from the observed between-group variance
    of the fitted means. Some implementations, Stata's
    \code{hetop.ado} among them, first subtract the estimation error in
    those means, and also apply a small-sample adjustment to the
    within-group component, so their standardized scale differs from
    this one by a small constant factor. Neither convention is more
    correct than the other, and quantities that do not depend on the
    scale (correlations, ratios of SDs, cut spacing ratios) are
    unaffected.}

  \item{gof}{For each school or district, a Pearson chi-square testing
    the goodness of fit (gof) of the normal distribution that the HETOP
    model assumes. A data frame with columns \code{chisq}, \code{df},
    \code{p}, and \code{min_exp} (the smallest expected count, for
    callers who wish to screen on it). Degrees of freedom are
    \code{K - 3}, where \code{K} is the number of bins, so at least
    four bins are needed for the test to have anything to test. The
    power of the test is low unless the school or district is large, so
    a non-significant p value does not establish that scores are
    normally distributed.}

  \item{iter_info}{Diagnostics from the fit, useful for checking that
    it behaved. Not needed for ordinary use. It records

    \itemize{
      \item how the per-group iteration went: \code{within_iterations}
        (total passes), \code{within_iterations_per_district},
        \code{within_converged}, \code{within_exact3} (groups solved
        in closed form because exactly three bins were populated), and
        \code{within_unidentified};
      \item when the cut scores were estimated, how that outer loop
        went: \code{iterated}, \code{iterations}, \code{converged},
        and the two tolerances used, \code{tol} and
        \code{tol_cutpoints_std};
      \item how any underidentified groups were salvaged:
        \code{salvage_counts}, the individual counts
        \code{salvaged_two_bin}, \code{salvaged_one_interior},
        \code{salvaged_one_bottom} and \code{salvaged_one_top}, and
        the borrowed values \code{mu_pool} and \code{sigma_pool};
      \item under \code{estimator = "EB_shrunk"}, the shrinkage
        quantities: \code{eb_shrink}, \code{eb_weights},
        \code{eb_tau2} and \code{eb_tau2_mu} (between-group variance
        of the log-SDs and of the means), and the shrinkage targets
        \code{eb_log_mean} and \code{eb_mean_mu}.
    }}
}

\references{
  Reardon S., Shear B.R., Castellano K.E. and Ho A.D. (2017).
  \dQuote{Using heteroskedastic ordered probit models to recover moments
    of continuous test score distributions from coarsened data,}
  \emph{Journal of Educational and Behavioral Statistics} 42(1):3--45.

  Lockwood J.R., Castellano K.E. and Shear B.R. (2018).
  \dQuote{Flexible Bayesian models for inferences from coarsened,
    group-level achievement data,}
  \emph{Journal of Educational and Behavioral Statistics.} 43(6):663--692.
}

\author{
  Paul T. von Hippel and David J. Hunter.
}

\examples{
set.seed(1001)
G   <- 10

## Let means and SDs vary across the groups.
mug       <- seq(from = -2.0, to = 2.0, length = G)
sigmag    <- seq(from =  2.0, to = 0.8, length = G)
cutpoints <- c(-1.0, 0.0, 0.8)
ng        <- rep(1000, G)
ngk       <- gendata_hetop(G, K = 4, ng, mug, sigmag, cutpoints)

## Here the counts were simulated by sampling ng scores per group, so
## scope = "sample" and the reported SEs will include sampling as well
## as binning error.
##
## Cutpoints known: both est_raw (test-score scale) and est_std
## (standardized scale) are returned.
bm <- fast_hetop(ngk, cutpoints_known = TRUE, cutpoints = cutpoints,
                 scope = "sample")
print(cbind(true = mug,       est = bm$est_raw$mean))
print(cbind(true = sigmag,    est = bm$est_raw$sd))
print(cbind(true = cutpoints, est = bm$est_raw$cutpoints))
print(cbind(est = bm$est_raw$mean,
            se  = bm$est_raw$mean_se,
            lo  = bm$est_raw$mean_ci_lower,
            hi  = bm$est_raw$mean_ci_upper))

## If the data represented the whole population, then scope =
## "population", and the reported standard errors are smaller because
## they reflect only binning error.
bm_pop <- fast_hetop(ngk, cutpoints_known = TRUE, cutpoints = cutpoints,
                     scope = "population")
print(cbind(sample_se     = bm$est_raw$mean_se,
            population_se = bm_pop$est_raw$mean_se))

## If cutpoints_known = FALSE, then cutpoints are unknown and are
## estimated on a scale with mean and SD given by pooled_mean
## (default 0) and pooled_sd (default 1).
bm2 <- fast_hetop(ngk, scope = "sample")
print(bm2$est_std$cutpoints)
print(bm2$est_std$mean)

## If the cut scores are unknown but the pooled mean and SD are known,
## you can pass the pooled mean and SD and the estimates will come back
## on that scale.
bm3 <- fast_hetop(ngk, pooled_mean = 1640.2, pooled_sd = 138.6,
                  scope = "sample")
print(bm3$est_std$cutpoints)
}

\keyword{models}
