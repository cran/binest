\name{mle_hetop}
\alias{mle_hetop}

\title{mle_hetop(): Joint Maximum Likelihood Estimation of the HETOP Model from Bin Counts}

\section{Deprecated}{
  \code{mle_hetop()} is deprecated in favor of
  \code{fast_hetop(estimator = "ML")}, which produces ML estimates
  much faster. \code{mle_hetop()} is preserved for comparison but will
  not be maintained.
}

\description{
  \code{mle_hetop()} returns maximum likelihood estimates from the
  HETOP model. Estimates are identical to those from
  \code{fast_hetop(estimator = "ML")}, but \code{mle_hetop()} runs much
  more slowly because it tries to estimate all districts at once and
  uses numerical approximations for the first and second derivatives.

  This implementation is forked, with small changes, from the
  \code{HETOP} package by J. R. Lockwood.
}

\usage{
mle_hetop(ngk, iterlim = 1500, ...)
}

\arguments{
  \item{ngk}{
    Data giving bin counts. Column \code{k} of row \code{g} counts how
    many students from school or district \code{g} scored in bin
    \code{k}.
  }

  \item{iterlim}{
    Maximum number of iterations used in optimization (passed to
    \code{\link{nlm}}).
  }

  \item{...}{
    Any other arguments for \code{\link{nlm}}.
  }
}

\details{
  This function requires at least 3 populated bins to identify the
  model. For districts with fewer than 3 populated bins, we identify
  the mean by assuming that the SD is equal to the geometric mean of
  the SDs for other districts.

  This function can report estimates on four different scales:

  \enumerate{
    \item the original estimation scale with two fixed cutpoints;
    \item a scale defined by forcing the group means and log group
      standard deviations each to have weighted mean of zero, where
      weights are proportional to the total count for each group;
    \item a scale where the population mean of the latent variable is
      zero and the population standard deviation is one; and
    \item a scale similar to (3) but where a bias correction is
      applied. See Reardon et al. (2017) for details on this bias
      correction.
  }

  The function also returns an estimated intracluster correlation
  (ICC) of the latent variable, defined as the ratio of the
  between-group variance of the latent variable to its marginal
  variance. Scales (1)-(3) above lead to the same estimated ICC;
  scale (4) uses a bias-corrected estimate of the ICC which will not
  in general equal the estimate from scales (1)-(3).
}

\value{
  A list with the following components:

  \item{est_fc}{A list of estimated group means, group standard
    deviations, cutpoints and ICC on scale (1).}

  \item{est_zero}{A list of estimated group means, group standard
    deviations, cutpoints and ICC on scale (2).}

  \item{est_star}{A list of estimated group means, group standard
    deviations, cutpoints and ICC on scale (3).}

  \item{est_starbc}{A list of estimated group means, group standard
    deviations, cutpoints and ICC on scale (4).}
  
  \item{nlmdetails}{The object returned by \code{\link{nlm}} that
    summarizes detailed of the optimization.}

  \item{pstatus}{A dataframe, with one row for each group, summarizing
    the estimation status of the mean and standard deviation for each
    group.  A value of \code{est} means that the parameter was estimated
    without constraints.  A value of \code{mean}, used for the group
    standard deviations, indicates that the parameter was constrained.
    Values of \code{min} or \code{max}, used for the group means,
    indicate that the parameter was constrained.}
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
  J. R. Lockwood (original implementation); David J. Hunter and
  Paul T. von Hippel \email{ph3828@eid.utexas.edu} (modifications for
  \pkg{binest}).
}

\examples{
set.seed(1001)

## define true parameters
G         <- 10
mug       <- seq(from= -2.0, to= 2.0, length=G)
sigmag    <- seq(from=  2.0, to= 0.8, length=G)
cutpoints <- c(-1.0, 0.0, 0.8)

## generate data with large counts
ng   <- rep(100000,G)
ngk  <- gendata_hetop(G, K = 4, ng, mug, sigmag, cutpoints)
print(ngk)

## compute MLE and check parameter recovery (cutpoints derived from data):
m    <- mle_hetop(ngk)
print(cbind(true = mug,       est = m$est_fc$mug))
print(cbind(true = sigmag,    est = m$est_fc$sigmag))
print(cbind(true = cutpoints, est = m$est_fc$cutpoints))

## estimates on other scales:
p    <- ng/sum(ng)
print(sum(p * m$est_zero$mug))
print(sum(p * log(m$est_zero$sigmag)))

print(sum(p * m$est_star$mug))
print(sum(p * (m$est_star$mug^2 + m$est_star$sigmag^2)))

## dealing with sparse counts
ngk_sparse <- matrix(rpois(G*4, lambda=5), ncol=4)
ngk_sparse[1,] <- c(5,8,0,0)
ngk_sparse[2,] <- c(0,10,10,0)
ngk_sparse[3,] <- c(12,0,0,0)
ngk_sparse[4,] <- c(0,0,0,10)
print(ngk_sparse)

m    <- mle_hetop(ngk_sparse)
print(m$pstatus)
print(unique(m$est_fc$sigmag[1:4]))
print(exp(mean(log(m$est_fc$sigmag[5:10]))))
print(m$est_fc$mug[3])
print(min(m$est_fc$mug[-3]))
print(m$est_fc$mug[4])
print(max(m$est_fc$mug[-4]))
}

\keyword{models}
