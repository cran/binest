# binest 0.3-0 (2026-09-21)

* **Breaking change.** The `iterate` argument is withdrawn from
  `fast_hetop()`. Whether the outer cutpoint loop runs was never a
  modelling choice: with cutpoints supplied there is nothing to refine,
  and with cutpoints estimated the refinement is always wanted. It is
  now derived internally as `!cutpoints_known`, which also removes the
  invalid `iterate = TRUE` with `cutpoints_known = TRUE` combination
  that previously had to be guarded against.

* **Breaking change.** The `robust` argument is withdrawn from
  `fast_hetop()`. Standard errors now always come from the
  model-implied Fisher information. A Huber-White sandwich remains
  implemented internally but is no longer reachable from the API: a
  sandwich corrects the variance, not the estimand, and where the
  within-group normal is wrong the damage is bias in the fitted mean,
  which a wider interval absorbs rather than measures. The `robust`
  field is no longer echoed back in `est_raw` and `est_std`, and the
  per-call message no longer reports it.

* `estimator` now has the literal default `"MLE"` rather than the
  `c("MLE", "EB_shrunk")` idiom, so the usage line states the default
  outright. Behavior is unchanged; `match.arg()` still validates
  against both values.


* `fast_hetop()`'s `gof` component now reports a chi-square test for
  **every** identified group with four or more bins, and gains a
  `min_exp` column.  Previously the Pearson sum ran only over bins
  with a positive observed count and degrees of freedom were taken
  from that count, so a four-bin group with one empty bin got `df = 0`
  and no test.  That screened on the wrong quantity: an observed zero
  is data (the fitted normal gives the bin positive probability, so
  seeing nobody there is evidence against the fit), while it is a
  small *expected* count that threatens the chi-square approximation.
  Degrees of freedom are now `K - 3`, counting bins rather than
  occupied bins, and `min_exp` lets the caller apply Cochran's rule or
  any other screen. Users comparing rejection rates with earlier
  versions should expect more groups to be tested.

## Explicit `cutpoints_known`; optional pooled moments for the reported scale

* **Breaking change.** `fast_hetop()` takes a new `cutpoints_known`
  argument (logical, default `FALSE`) stating which mode the call is in,
  rather than inferring it from whether `cutpoints` happens to be
  `NULL`. With `cutpoints_known = TRUE`, `cutpoints` is required and its
  length is checked against the number of bins in `ngk`
  (`length(cutpoints)` must be `K - 1`). With `cutpoints_known = FALSE`,
  supplying `cutpoints` is an error. Contradictory combinations are
  rejected rather than silently resolved.

* New `pooled_mean` and `pooled_sd` arguments, used only when
  `cutpoints_known = FALSE`, giving the mean and SD to place the pooled
  distribution at on the reported scale. Defaults `0` and `1` reproduce
  the previous standardized output. Other values are useful when the cut
  scores are unknown but the pooled moments are published: the group
  estimates and the estimated cutpoints then come back on that scale
  directly, instead of the caller rescaling afterwards. Both are echoed
  back in `est_std`.

* The cutpoint-refinement EM's tolerance is no longer a user argument
  (see the convergence-criteria section below). It is worth knowing that
  it is in *standardized* units regardless of how results are reported:
  the EM always iterates on the internal latent scale, and
  `pooled_mean`/`pooled_sd` are applied only afterwards as an affine
  relabeling of the converged solution. The likelihood is invariant to
  that reparameterization, so the fit is identical either way and only
  the reported numbers change.

* `iter_info` gains `within_iterations_per_district`, the per-group
  iteration counts (accumulated across EM passes when the cutpoints are
  estimated). Previously only the maximum survived, which made it
  impossible to tell whether a long run reflected many slow groups or
  one pathological one.

## New required `scope` argument; one SE/CI instead of four

* **Breaking change.** `fast_hetop()` now takes a required `scope`
  argument, either `"sample"` or `"population"`, with no default. Calls
  that omit it stop with an error explaining the choice. Which standard
  errors are appropriate depends on whether the units behind each
  group's counts are a sample from a larger population or that group's
  entire population, and that is a fact about the data that only the
  caller knows — defaulting it would silently return the wrong
  uncertainty.

    * `scope = "sample"`: the estimand is the parameter of the
      population the units were drawn from, so the error has two
      components. Returns the total SE, `Var_total = Var_sampling +
      Var_binning`.
    * `scope = "population"`: the units *are* the group's population, so
      its true mean is the actual mean of their scores and there is no
      sampling error. Returns the binning-only SE, `Var_total −
      Var_sampling`.

* **Breaking change.** The output now carries exactly one standard error
  and one confidence interval per parameter — `mean_se`, `sd_se`,
  `mean_ci_lower`/`_upper`, and `sd_ci_lower`/`_upper` — replacing the
  `_total`, `_sampling`, `_binning` and `_sandwich` variants. All are
  still computed internally; only the one matching `scope` is
  returned, so results cannot be read with the wrong notion of
  uncertainty. The choice is echoed back as a `scope` field, and
  reported via `message()` on each call. Code needing more than one flavor should
  fit more than once; the point estimates are identical either way.

* New `robust` argument (logical, default `FALSE`), orthogonal to
  `scope`: `scope` selects which error sources the SE covers, `robust`
  selects how the variance is estimated — model-implied Fisher
  information, or a Huber-White sandwich whose meat uses the observed
  bin counts. All four combinations are supported. Under `scope =
  "population"` the same `Var_total − Var_sampling` subtraction is
  applied to the sandwich; because the sandwich equals the Fisher
  variance only in expectation, that can go negative, and affected
  groups fall back to the Fisher binning-only variance with the count
  reported in `iter_info$robust_binning_fisher_fallback_mu` and `_sd`.
  `robust = TRUE` warns under `estimator = "EB_shrunk"`, where the
  sandwich has no analog.

* Documentation now notes that `est_std` standardizes by the *observed*
  between-group variance, whereas some implementations (Stata's
  `hetop.ado`) first remove the estimation error in the group means and
  adjust the within-group component for small samples. The two
  standardized scales therefore differ by a small constant factor;
  neither is more correct, and scale-free quantities are unaffected.

## Estimates are named `mean` and `sd`

* **Breaking change.** The per-group estimates in `est_raw` and
  `est_std` are now named `mean` and `sd`, with `mean_se`, `sd_se`,
  `mean_ci_lower`/`_upper` and `sd_ci_lower`/`_upper` alongside them.
  Through 0.2-1 they carried an estimator-dependent suffix
  (`group_mean_mle` and `group_sd_mle`, or `group_mean_eb` and
  `group_sd_eb`), which forced callers to build the variable name from
  the estimator they had asked for. Which estimator produced them is
  now reported in the `estimator` element instead, so the same code
  reads the result whichever estimator was used. The `group_` prefix
  is dropped as redundant inside a per-group result.

* Note that `mean` and `sd` shadow the base functions of those names
  when the result list is used inside `with()`: `with(fit$est_raw,
  sd(mean))` fails, because `sd` there is the vector, not the
  function. Extract with `$` if you need both.

## Parameter-change convergence criteria; a single `tol`

* **Behavior change.** The per-group iteration now stops on a
  *parameter* change rather than a likelihood change: a group is done
  once both `|Δμ_g| / σ_g < tol` and `|Δ log σ_g| < tol`, with `tol`
  defaulting to `1e-4`. Both terms are scale-free, so `tol` means the
  same thing whether the internal scale is test-score points or the
  standardized latent scale, and a whole fit can be described as
  "settled to `tol` SD". Note the mean term is in units of each group's
  *own* σ, not the pooled SD.

* **`tol` is now the only tolerance argument.** The outer
  cutpoint-refinement EM has its own tolerance, but it is derived
  internally as `10 * tol` rather than exposed. It has to be looser than
  the per-group one: the cutpoint solve depends on group fits converged
  only to `tol`, so their residual wobble reaches the cutpoints at
  roughly the same order, and a threshold at or below that floor can
  never be met — the EM then runs to `maxit` however large `maxit` is.
  Since that relationship is a property of the algorithm rather than a
  modelling choice, there is nothing for a caller to decide and no way
  to set an inconsistent pair.

* The previous criterion — change in the per-observation log likelihood
  — had two failure modes. It stopped early on flat ridges, where the
  objective settles while μ and σ are still drifting, which is common in
  small groups whose mass sits in an open-ended extreme bin and whose σ
  is weakly identified. And because it was per-observation, its implied
  tolerance on a group's total log likelihood was `n_g * tol`, looser
  for larger groups.

* Note that a parameter-change rule bounds the last *step*, not the
  distance remaining. For an iteration contracting at rate ρ, what
  remains is on the order of step/(1−ρ), so on flat ridges the achieved
  accuracy can be an order of magnitude looser than `tol`. For EM, ρ is
  the fraction of missing information, so the slowest groups are exactly
  those whose bin counts say least about σ.

* `iter_info` records both the `tol` supplied and the derived cutpoint
  tolerance, so a fitted object can report what it actually used.

* The per-group iteration no longer evaluates the log likelihood at all.
  The update is an EM step — E-step conditional moments, closed-form
  M-step — so the objective never had to be computed except for the old
  stopping test, and dropping it saves a density and CDF evaluation per
  group per iteration.

## `tol` now controls the per-district MLE; new `tol_cutpoints` for the EM

* **Behavior change.** `fast_hetop()`'s `tol` argument previously
  governed only the outer cutpoint-refinement EM, while the per-district
  MLE in `.refine_within()` used a convergence criterion fixed
  internally at `1e-6`. `tol` now controls the per-district criterion
  and defaults to `1e-10`; the EM tolerance moves to a new `tol_cutpoints`
  argument, defaulting to `5e-3` so that fits with estimated cutpoints
  are unchanged.

* The per-district criterion is a change in the **per-observation** log
  likelihood, so the implied tolerance on a district's total log
  likelihood is `n_g * tol`, and the residual shortfall below the true
  maximum grows with district size. At the old `1e-6` that was around
  `7e-5` for a median-sized district. Comparing against Stata's
  `hetop.ado` (Newton-Raphson, `ltolerance` on the total) across 1,134
  Texas districts, `hetop.ado` attained a marginally higher likelihood
  in about 90% of them, by a mean of `1.7e-05` and a maximum of
  `4.5e-03` -- margins far too small to move any estimate (the two
  agreed to 0.004 scale-score points on the mean), but systematic.
  Estimates from the new default should be at or above `hetop.ado`'s
  likelihood in nearly all districts.

* `tol` is also now passed to `.refine_mu_only()`, used to salvage
  two-populated-bin districts, which likewise had `1e-6` hard-coded.

## Huber-White sandwich SE (fourth SE flavor)

* `fast_hetop()` now reports a fourth per-district standard-error
  flavor alongside `_total`, `_sampling`, and `_binning`: a Huber-White
  sandwich (robust / model-misspecification-robust) SE, exposed as
  `group_mean_se_sandwich` and `group_sd_se_sandwich` with matching
  Wald / log-scale CIs `group_mean_ci_lower_sandwich`/`_upper_sandwich`
  and `group_sd_ci_lower_sandwich`/`_upper_sandwich`. The sandwich
  variance uses observed bin counts (not model probabilities) in the
  "meat" of the sandwich, `V_sandwich = I_matrix^{-1} J_matrix
  I_matrix^{-1}` where `J = sum_k n_{g,k} * s_k s_k^T` and
  `I_matrix = n_g * I` is the per-district Fisher information. Under
  a correctly specified within-group normal, `E[J] = n_g * I` and the
  sandwich collapses to the Fisher `_total` variance; under
  misspecification, the sandwich SE is calibrated even when
  within-group scores are not exactly normal, so it is the
  appropriate SE to report when the goodness-of-fit test rejects
  normality. All four sandwich fields (two SEs and two CI pairs)
  appear in both `est_raw` and `est_std`. Under
  `estimator = "EB_shrunk"` the sandwich framework does not cleanly
  apply to a shrunk posterior SD, so the sandwich fields are `NA`.
  For salvaged districts, the sigma-known sandwich SE
  `J11 / (n_g * I11)^2` is used for `mu` on two-bin and
  one-interior-bin salvages, and the sandwich SE is `NA` for
  one-extreme-bin salvages and (for `sigma`) for every salvaged
  district.

## Opt-out for the underidentified-district salvage

* `fast_hetop()` gains a new logical argument
  `estimate_unidentified_districts`, defaulting to `TRUE`. When `TRUE`
  (the default, and the behavior of every 0.2-3 development release
  until now), districts with fewer than three populated bins are
  salvaged by borrowing `log(sigma)` from the pool of identified
  districts. When `FALSE`, those districts are returned as `NA`,
  reproducing the pre-0.2-3 behavior of the package for callers who
  prefer to handle underidentified districts themselves.

## Two flavors of confidence interval per parameter

* `fast_hetop()` now reports two confidence intervals per per-district
  parameter instead of one, corresponding to the two sampling-model
  interpretations of the input data. The former
  `group_mean_ci_lower`/`_upper` is renamed
  `group_mean_ci_lower_total`/`_upper_total` and joined by
  `group_mean_ci_lower_binning`/`_upper_binning`; the SD CIs are
  renamed and joined analogously. The `_total` CI is built on the
  total (sampling+binning) SE and is the appropriate interval for
  sample data where the district's students are drawn from a
  within-district superpopulation. The `_binning` CI is built on the
  binning-only SE (using `mug_se_binning` for the mean and
  `log_sigmag_se_binning` on the log scale for the SD, then
  exponentiated back) and is the appropriate interval for population /
  administrative data (e.g. state education-agency reports) where the
  district's students are the full population and only the coarsening
  into bins is a real unknown. Both intervals appear in both `est_raw`
  and `est_std` at the requested `conf.level`. Under `estimator =
  "EB_shrunk"` the `_binning` SE is NA, so the `_binning` CI is NA as
  well; likewise for salvaged districts whose underlying binning SE
  is NA.

## Sampling-vs-binning decomposition of standard errors

* `fast_hetop()` now reports each per-district standard error as three
  fields instead of one, so a reader can see how much of the SE comes
  from finite sample size vs. from having binned rather than raw
  scores. The former `group_mean_se` is renamed `group_mean_se_total`
  and joined by `group_mean_se_sampling` (what the SE would be from
  `n_g` continuous unbinned normal observations,
  `sigma_g / sqrt(n_g)`) and `group_mean_se_binning` (the residual,
  the extra SE attributable to coarsening those scores into `K` bins).
  `group_sd_se` is renamed and decomposed analogously via the
  log-scale delta method, using `1 / (2 n_g)` for the log-SD sampling
  variance. Under the Pythagorean identity
  `Var_total = Var_sampling + Var_binning` (binning cannot add
  information), the three add in quadrature. All six fields appear in
  both `est_raw` and `est_std`. Confidence intervals are still built
  on the `_total` SE; only the field name inside the CI construction
  changed. Under `estimator = "EB_shrunk"` the `_total` SE remains the
  approximate empirical-Bayes posterior SD, and the sampling/binning
  split does not apply cleanly, so `_sampling` and `_binning` are
  returned as NA. For salvaged districts (`< 3` populated bins), sigma
  is imputed rather than estimated, so all three sigma SE components
  are NA; mu SE components are populated normally for two-bin and
  one-interior-bin salvages and NA for one-extreme-bin salvages.

## Salvage of districts with fewer than three populated bins

* `fast_hetop()` now recovers estimates for districts with only one
  or two populated bins, which previously returned `NA`. Sigma is
  borrowed from the identified pool as
  `sigma_pool = exp(mean(log(sigma_g)))` over districts with at
  least three populated bins (unweighted per district, matching
  `mle_hetop()`'s pooling), and mu is then estimated with sigma
  fixed. Four sub-rules cover the possible populated-bin patterns:
  a fixed-point iteration on mu alone for two-bin districts; a
  closed-form `(cutpoints[j-1] + cutpoints[j]) / 2` for one interior
  populated bin; and the running max (resp. min) of mu over all
  identified and already-salvaged districts for a single populated
  bin at the top (resp. bottom) bin. Salvaged districts flow through
  the empirical-Bayes shrinkage, Fisher-information SEs, and
  confidence-interval machinery unchanged; their goodness-of-fit
  chi-square remains `NA` (the test requires >= 4 populated bins).
  The salvage counts under each rule and the pooled sigma are
  reported in `iter_info$salvage_counts` and
  `iter_info$sigma_pool`. If no district in the input is identified,
  salvage is skipped with a warning and the unidentified rows remain
  `NA` as before.

# binest 0.2-2 (unreleased)

## Cutpoint-refinement algorithm rewritten

* When `cutpoints = NULL`, cutpoints are still initialized to
  `qnorm(cumsum(pooled_props))`, but the EM refinement now alternates
  between (a) updating each cutpoint `C_i` to the quantile of the
  current mixture-of-normals CDF at the pooled cumulative proportion
  `P_i`, and (b) refitting each district's `(mu_g, sigma_g)` to its
  full truncated-normal MLE at the new cutpoints. Previously, step
  (b) used a one-step moment shortcut that computed each district's
  `(mu_g, sigma_g)` as a weighted average of *standard-normal*
  truncated-bin moments -- a biased approximation that pulled
  per-district means and SDs toward the pooled mean and SD. The
  shortcut function `.fast_hetop_oneshot()` has been removed.

* `tol` default relaxed from `1e-3` to `5e-3`. The new EM does a
  full per-district MLE at every iteration, so each outer step is
  more expensive; a slightly looser cutpoint tolerance keeps overall
  runtime comparable while still resolving cutpoints to a couple of
  decimal places on the standardized scale.

* `est_raw` is now omitted from the return value when
  `cutpoints = NULL`. The EM lands on some arbitrary location and
  scale for the cutpoints, and reporting per-district estimates on
  that arbitrary scale alongside `est_std` was duplicate output on a
  scale with no external interpretation. When cutpoints are supplied,
  both `est_raw` (test-score scale) and `est_std` (standardized
  scale) are returned as before.

## `bin_means()` renamed to `fast_hetop()`

* The function formerly exported as `bin_means()` is now exported as
  `fast_hetop()`, better reflecting that it recovers essentially the
  same estimates as `mle_hetop()` (and `fh_hetop()`) in closed form.
  `mle_hetop()` and `fh_hetop()` are now deprecated in favor of
  `fast_hetop()` and remain in the package only for comparison.

## New output from `fast_hetop()`

* `fast_hetop()` now returns `group_mean_se` and `group_sd_se` in both
  `est_raw` and `est_std`: the Fisher-information sampling SE of the
  unshrunk per-district MLE (`estimator = "MLE"`), or an approximate
  empirical-Bayes posterior SE `sqrt(w_g * s2_g)` (`estimator = "EB_shrunk"`).
  The underlying Fisher-information calculation already existed
  internally to set the empirical-Bayes shrinkage weights; it was
  simply never returned to the caller. `NA` for unidentified
  districts. Calibration of these SEs (via standardized residuals in
  a simulation) is reported in `code/normal_simulation.R`.

* `fast_hetop()`'s `iterate` argument now defaults to
  `is.null(cutpoints)` rather than always `FALSE`: cutpoint iteration
  is on by default whenever cutpoints are being estimated from the
  data (since it is cheap), and off by default when cutpoints are
  supplied (where it was already invalid and is now simply not
  attempted rather than requiring the caller to pass
  `iterate = FALSE` to avoid an error). This is a behavior change for
  existing callers who relied on the implicit default and did not pass
  `iterate` explicitly: estimated-cutpoints calls will now iterate
  unless `iterate = FALSE` is passed. Iteration did not meaningfully
  change per-group estimates in testing on the Texas data or in
  simulation, so this is not expected to change results in typical
  use, but it is a default-value change worth noting in case it does
  in some edge case.

* `fast_hetop()` gains a `conf.level` argument (default `0.95`) and now
  also returns `group_mean_ci_lower`/`_upper` and
  `group_sd_ci_lower`/`_upper` in both `est_raw` and `est_std`. The
  mean's interval is the ordinary symmetric Wald interval around
  `group_mean_se`. The SD's interval is instead constructed on the log
  scale and exponentiated back, which is asymmetric around `group_sd`
  but better calibrated: the sampling distribution of an estimated SD
  is right-skewed in finite samples (most noticeably for small
  groups), and simulation showed the naive symmetric interval
  under-covers as a result, while the log-scale interval does not.

## API changes to `fh_hetop()`

* The `fixedcuts` argument has been removed, matching the fix already
  applied to `mle_hetop()` in 0.2-0. Cutpoints are now derived
  internally from the pooled bin proportions via
  `qnorm(cumsum(colSums(ngk)/sum(ngk))[1:2])`. `fixedcuts` was never a
  channel for supplying cutpoints that are known in advance on the
  test-score scale (e.g. published cut scores): the model is
  unidentified without pinning down the location and scale of the
  latent normal somehow, and fixing two cutpoints was only ever that
  normalization, not a constraint reflecting external knowledge about
  the data. Supplying cutpoints on the native test-score scale rather
  than the standardized scale the model uses internally reliably sent
  the MCMC sampler into a degenerate region it could not escape
  (cell probabilities driven to 0 or 1 in floating point), the same
  failure mode documented for `mle_hetop()` below. Deriving the
  cutpoints internally removes both the failure mode and the mistaken
  impression that known cutpoints could be supplied this way.

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
  `group_mean_mle` and `group_sd_mle` when `estimator = "MLE"`, and
  `group_mean_eb` and `group_sd_eb` when `estimator = "EB_shrunk"`, in place
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
