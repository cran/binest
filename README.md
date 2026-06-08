# binest

Estimation of group-level means and standard deviations from binned
(coarsened) count data. The package implements three methods with a
common output structure:

* `bin_means()` — fast per-group estimator under within-group
  normality. Linear in the number of groups times the number of bins.
* `mle_hetop()` — maximum-likelihood fit of the heteroskedastic
  ordered probit (HETOP) model.
* `fh_hetop()` — Bayesian variant of HETOP via MCMC.

This package was previously called HETOP and was maintained by J. R.
Lockwood; it is renamed and extended to reflect the broader set of
estimators now included.

## Installation

```r
# install.packages("remotes")
remotes::install_github("paulvonhippel/binest")
```

(Or `install.packages("binest")` once on CRAN.)

## A quick example

```r
library(binest)
data(tx_g6_math_2018)

ngk <- with(tx_g6_math_2018,
            cbind(unsatisfactory, approaches, meets, masters))
cuts <- c(1536, 1653, 1772)

fit <- bin_means(ngk, cutpoints = cuts)
cor(fit$est_raw$group_mean_mle, tx_g6_math_2018$reported_mean)
```

See `vignette("binest")` for a full comparison of the three estimators
on the Texas STAAR Grade-6 mathematics data.

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
