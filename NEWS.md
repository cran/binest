# binest 0.2-1 (resubmission)

Corrections requested by CRAN reviewer Benjamin Altmann:

* DESCRIPTION: function names now written with parentheses (e.g.
  `bin_means()`); package and software names ('HETOP', 'CRAN') now
  in single quotes per CRAN style.
* `man/waic_hetop.Rd`: replaced commented-out example with a
  runnable example wrapped in `\donttest{}`.
* `R/fh_hetop.R`: suppressed JAGS' default console output. The
  progress bar now defaults to `"none"` rather than `"text"`, and
  the `R2jags::jags()` call is wrapped in `capture.output()` so
  JAGS' C++ initialization messages do not write to stdout.

# binest 0.2-0 (2026-05-29)

## API changes to `bin_means()`

The `bin_means()` function signature has been simplified.

* The `within` argument has been removed; the function now always uses
  per-district maximum-likelihood estimation (the former
  `within = TRUE` behavior).
* The `sampling_var` argument has been removed; empirical-Bayes
  shrinkage always uses the Fisher information from the binned-normal
  likelihood (the former `sampling_var = "fisher"` behavior).
* The default convergence tolerance `tol` is now `1e-3` on the
  standardized scale, in place of the earlier `1e-6`.

## Output structure

* `bin_means()` now returns components `est_raw` and `est_std` in
  place of the earlier `est_fc`, `est_zero`, and `est_star`.
* The per-group estimates inside `est_raw` and `est_std` are named
  `group_mean_mle` and `group_sd_mle` when `eb_shrink = FALSE`, and
  `group_mean_eb` and `group_sd_eb` when `eb_shrink = TRUE`, in place
  of the earlier `mug` and `sigmag`.
* `bin_means()` now also returns a per-group `gof` data frame with the
  chi-square statistic, degrees of freedom, and p-value of a
  goodness-of-fit test for the within-group normality assumption.

## Cleaner unidentified-district handling

* Districts with fewer than three populated bins are not jointly
  identified under within-group normality.  `bin_means()` now returns
  \code{NA} for both the mean and the SD of such districts, matching
  the methodological recommendation of Reardon et al. (2017).

# binest 0.1-0 (2026-05-27)

## Renamed from HETOP

This package was previously called `HETOP`, originally authored and
maintained by J.R. Lockwood (last CRAN release: HETOP 0.2-6, June 2019;
archived from CRAN in March 2025).  It is renamed `binest` to reflect
the broader scope of methods now included for estimating distributional
moments from binned (coarsened) count data.

## New estimator

* `bin_means()` — a fast per-group estimator of means and SDs under
  within-group normality.  Uses the pooled bin proportions to derive
  standardized cutpoints, then weights truncated-normal moments
  within each bin by within-group bin proportions.  Returns estimates
  on the same four scales as `mle_hetop()`.  Runs in `O(GK)` time.

## Bug fix in `mle_hetop()` carried over from the HETOP 0.3-0 patch

`mle_hetop()` previously took two user-facing arguments, `fixedcuts`
and `svals`.  Internal starting values for group means and log SDs are
on a standardized scale; when the supplied cutpoints were on a native
test-score scale (e.g. STAAR scale scores like 1536 and 1653), the log
likelihood was numerically flat near the starting values, `nlm()`
exited with zero iterations, and the function silently returned the
starting values disguised as MLE estimates.

`mle_hetop()` no longer accepts `fixedcuts` or `svals`.  Cutpoints are
derived internally from the pooled bin proportions via

    qnorm(cumsum(colSums(ngk)/sum(ngk))[1:2])

so the cutpoints, starting means, and starting log SDs all live on the
same standardized scale.  The function warns if `nlm()` reports zero
iterations.  The cell-probability computation is also vectorized for
roughly an 18x speedup per likelihood evaluation.

The bug was first identified by David J. Hunter in January 2022 and
independently reproduced by Benjamin R. Shear in July 2024.

## History of HETOP

The `HETOP` package was originally authored by J.R. Lockwood, with
substantial methodological contributions from Sean F. Reardon,
Benjamin R. Shear, Katherine E. Castellano, and Andrew D. Ho.  The
implementations of `fh_hetop`, `gendata_hetop`, `triple_goal`, and
`waic_hetop` in this package are unchanged from HETOP 0.2-6
(J.R. Lockwood, June 2019).
