# binest

Estimation of group-level means and standard deviations from binned
(coarsened) count data. All three functions fit the same
heteroskedastic ordered probit (HETOP) model, in which each group's
values are normally distributed around a mean and SD of its own. They
differ only in how they fit it, and they share a common output
structure:

* `fast_hetop()` — fits each group separately, solving closed-form
  truncated-normal score equations. Linear in the number of groups
  times the number of bins. This is the preferred function in the
  package.
* `mle_hetop()` — maximizes the likelihood over all groups at once.
  Returns the same estimates as `fast_hetop(estimator = "ML")`, far
  more slowly. Deprecated in favor of `fast_hetop()`.
* `fh_hetop()` — fits the model by MCMC, placing a prior over the
  group parameters and reporting posterior means. Deprecated in favor
  of `fast_hetop()`.

This package was previously called HETOP and was maintained by J. R.
Lockwood; it is renamed and extended to reflect the broader
functionality now included.

## Installation

```r
# install.packages("remotes")
remotes::install_github("paulvonhippel/binest")
```

(Or `install.packages("binest")` for the released version on CRAN.)

## A quick example

```r
library(binest)
data(tx_g6_math_2018)

ngk <- with(tx_g6_math_2018,
            cbind(unsatisfactory, approaches, meets, masters))
cuts <- c(1536, 1653, 1772)

## `scope` is required: Texas reports counts for every tested student,
## so each district's students are its whole population and the only
## uncertainty is the binning. Use scope = "sample" when the units are
## a sample from a larger population, which adds sampling error to the
## reported SEs.
fit <- fast_hetop(ngk, cutpoints_known = TRUE, cutpoints = cuts,
                  scope = "population")
cor(fit$est_raw$mean, tx_g6_math_2018$reported_mean)
```

See `vignette("binest")` for a full comparison of the three functions
on the Texas STAAR Grade 6 mathematics data.

## References

* Fisher, R. A. (1922). On the mathematical foundations of theoretical
  statistics. *Philosophical Transactions of the Royal Society of
  London A*, 222, 309-368.
* Lockwood, J. R., Castellano, K. E., & Shear, B. R. (2018). Flexible
  Bayesian models for inferences from coarsened, group-level
  achievement data. *JEBS*, 43(6), 663-692.
* Reardon, S. F., Shear, B. R., Castellano, K. E., & Ho, A. D. (2017).
  Using heteroskedastic ordered probit models to recover moments of
  continuous test score distributions from coarsened data. *JEBS*,
  42(1), 3-45.
* Sheppard, W. F. (1898). On the calculation of the most probable
  values of frequency-constants for data arranged according to
  equidistant divisions of a scale. *PLMS*, 29, 353-380.

## License

GPL (>= 2).
