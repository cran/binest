###############################################################################
## fast_hetop: fast estimator of group means and SDs from binned count
## data, under the assumption that within each group the underlying
## scores are normally distributed. Recovers the same estimates as the
## heteroskedastic ordered probit (HETOP) model in closed form, without
## the iterative/MCMC machinery mle_hetop()/fh_hetop() use.
##
## Arguments:
##   ngk        GxK numeric matrix of bin counts.
##   cutpoints_known
##              Logical, default FALSE.  States which of the two modes
##              the call is in, and the other cutpoint arguments are
##              checked against it rather than inferred:
##
##                TRUE   `cutpoints` is required and must contain the
##                       K-1 cut scores on the test-score scale.
##                       `pooled_mean` and `pooled_sd` are rejected,
##                       since the cut scores already fix the scale.
##                FALSE  `cutpoints` must be absent; the cutpoints are
##                       estimated from the pooled bin proportions.
##                       `pooled_mean` and `pooled_sd` set the scale the
##                       results are reported on.
##
##              Contradictory combinations are errors rather than being
##              silently resolved.
##   cutpoints  Numeric vector of length K-1: the known cut scores on
##              the test-score scale.  Required when
##              cutpoints_known = TRUE and rejected otherwise.  Must be
##              strictly increasing with no missing values; its length
##              is checked against the number of bins in `ngk`.
##   pooled_mean, pooled_sd
##              Used only when cutpoints_known = FALSE: the mean and SD
##              to give the pooled (population-weighted, within +
##              between) distribution on the reported scale.  Default
##              (0, 1), the usual standardized scale.
##
##              Supplying other values is useful when the cut scores are
##              unknown but the pooled moments are not -- an agency may
##              publish an overall mean and SD without publishing the
##              cut scores.  Passing them puts the group estimates and
##              the estimated cutpoints directly on that scale instead
##              of leaving the caller to rescale afterwards.
##   scope      REQUIRED, no default: "sample" or "population".  States
##              whether the units behind each group's counts are a sample
##              from a larger population or that group's entire
##              population, and thereby which standard errors and
##              confidence intervals are returned.
##
##                "sample"     the estimand is the parameter of the
##                             population the units were drawn from, so
##                             the error has two parts -- drawing these
##                             units rather than others, and observing
##                             only their bin counts.  Returns the TOTAL
##                             SE, satisfying
##                             Var_total = Var_sampling + Var_binning.
##                "population" the group's units ARE the population, so
##                             its true mean is the actual mean of those
##                             units' scores and there is no sampling
##                             error.  The only uncertainty is the
##                             coarsening.  Returns the BINNING-ONLY SE,
##                             Var_total - Var_sampling, with
##                             Var_sampling = sigma^2/n_g for the mean
##                             and 1/(2 n_g) for log sigma.
##
##              There is no sensible default: which one applies is a
##              fact about the data that only the caller knows, and
##              guessing would silently return the wrong uncertainty.
##              Population and administrative extracts are "population";
##              survey samples and simulated draws are "sample".
##   tol        Convergence tolerance for the per-district iteration in
##              `.refine_within` (default 1e-4).  A district stops once
##              BOTH
##
##                 |Delta mu_g| / sigma_g < tol
##                 |Delta log sigma_g|    < tol
##
##              between successive iterations -- that is, once its
##              estimates have settled to within `tol` of its own SD.
##
##              This is the ONLY tolerance a caller sets.  When the
##              cutpoints are estimated, the outer EM's tolerance is
##              derived from it internally (10 * tol, see the note where
##              `tol_cutpoints_std` is computed): the outer loop has to
##              be looser than the inner one, and that relationship is a
##              property of the algorithm rather than something to
##              decide, so there is no second knob and no way to set an
##              invalid pair.  The whole fit can then be described as
##              "settled to `tol` SD".
##
##              Both terms are scale-free, so `tol` means the same thing
##              whether the internal scale is test-score points (when
##              cutpoints are supplied) or the standardized latent scale
##              (when they are estimated).  Note that the mean term is
##              in units of each district's OWN sigma, not the pooled
##              SD.
##
##              Caveat: this bounds the size of the last STEP, not the
##              distance still to travel.  For an iteration contracting
##              at rate rho per pass, what remains is on the order of
##              step/(1 - rho), so on a flat ridge -- a small district
##              with most of its mass in an open-ended extreme bin,
##              whose sigma is weakly identified -- the accuracy reached
##              can be an order of magnitude looser than `tol`.  For EM,
##              rho is the fraction of missing information, so the
##              districts that converge slowest are exactly those whose
##              bin counts say least about sigma.
##
##              Earlier versions stopped on the change in the
##              per-observation log likelihood instead.  That stopped
##              early on those same flat ridges, where the objective
##              settles while mu and sigma are still drifting, and its
##              implied tolerance on a district's TOTAL log likelihood
##              was n_g * tol, hence looser for larger districts.
##   maxit      Iteration cap (default 100).  Applies both to the outer
##              EM and to each district's inner iteration.
##   conf.level Confidence level for mean_ci_*/sd_ci_*
##              (default 0.95).
##   estimate_unidentified_districts
##              Logical.  If TRUE (default), districts with fewer than
##              three populated bins are salvaged by borrowing log(sigma)
##              from the pool of identified districts and estimating
##              mu (or setting it to a heuristic value for one-extreme-
##              populated-bin districts).  If FALSE, those
##              districts are returned as NA, reproducing the pre-0.2-3
##              behavior of the package.
##
## Output: a list with components est_raw, est_std, gof, iter_info.
##   est_raw   estimates on the raw (test-score) scale, when cutpoints
##             were supplied; otherwise absent (the EM-converged scale
##             has no external interpretation, so only est_std is
##             returned). Carries exactly ONE standard error per
##             parameter -- mean_se and sd_se -- and ONE
##             confidence interval per parameter at the requested
##             conf.level -- mean_ci_lower/_upper and
##             sd_ci_lower/_upper. Which quantity those hold is
##             determined by `scope`, echoed back as
##             est_raw$scope so a saved object is
##             self-describing. Earlier versions returned several
##             flavors side by side (_total, _sampling, _binning,
##             _sandwich); they are still computed internally but only
##             the requested one is returned, so the output cannot be
##             read with the wrong notion of uncertainty.
##
##             Under estimator = "ML" the underlying total variance is
##             from the Fisher information for the binned-normal model,
##             and the sampling component it is compared against is the
##             variance that would apply to n_g continuous unbinned
##             observations (sigma_g/sqrt(n_g) for mu, and
##             sigma_g/sqrt(2 n_g) for sigma via the log-scale delta
##             method). Under estimator = "EB_shrunk" the total is the
##             approximate empirical-Bayes posterior SE sqrt(w_g * s2_g);
##             the sampling/binning split does not apply cleanly there,
##             so scope = "population" returns NA. All fields are NA for
##             unidentified districts.
##
##             Mean CIs are symmetric Wald intervals mug +/- z * SE; SD
##             CIs are built on the log scale around log(sd) and
##             exponentiated back, so they are asymmetric around
##             sd.
##   est_std   the same quantities on the standardized scale, where the
##             population-weighted state mean is 0 and the total
##             (within + between) state SD is 1. Same field names,
##             rescaled. Note this total SD is computed from the
##             observed between-group variance; some implementations
##             (Stata's hetop.ado among them) subtract the estimation
##             error in the group means first, which makes their
##             standardized scale differ from this one by a small
##             constant factor. Neither convention is more correct, and
##             scale-free quantities are unaffected.
##   gof       per-district Pearson chi-square goodness-of-fit of the
##             within-district normality assumption (chisq, df, p,
##             min_exp).  df = K - 3, counting bins rather than
##             occupied bins, so every identified district with K >= 4
##             is tested.  min_exp is the smallest expected count, for
##             callers who want Cochran's screen.
##   iter_info diagnostic flags including `eb_shrink` and (when
##             applicable) EM convergence counts and shrinkage tuning
##             parameters.
###############################################################################

.district_3bin_start <- function(ngk, cutpoints) {
    ## Starting values for each district's (mu_g, sigma_g) ML search.
    ## Used for every district on every call (both known-cutpoints and
    ## data-derived-cutpoints modes), and inside every EM iteration.
    ## No pooled quantity of any kind is used: each district is seeded
    ## entirely from its own bin counts.
    ##
    ## For a district with >= 3 populated bins (the package's existing
    ## identifiability threshold), collapse its own bins to 3
    ## categories using its own leftmost and rightmost POPULATED bins
    ## (not necessarily bin 1 and bin K) as the two surviving
    ## boundaries, pooling everything in between into a middle
    ## category. Because the two boundary bins are populated by
    ## construction, and the >= 3-populated-bin threshold guarantees at
    ## least one more populated bin strictly between them, all three
    ## collapsed-category proportions are strictly between 0 and 1 --
    ## so matching them exactly to a normal CDF via qnorm always gives
    ## a finite, closed-form solution to the otherwise-transcendental
    ## two-parameter equations (no continuity correction needed, and
    ## none of the failure modes that a raw bin-1/bin-K version would
    ## have when an extreme bin happens to be empty).
    ##
    ## Districts with < 3 populated bins are unidentified under the
    ## package's rule and will be set to NA later regardless, so they
    ## (and the rare case where the matched quantiles are numerically
    ## indistinguishable) get a generic cutpoint-range fallback: this
    ## only needs to hand `.refine_within` *some* finite value, never a
    ## good one.
    G <- nrow(ngk)
    K <- ncol(ngk)
    mu0    <- numeric(G)
    sigma0 <- numeric(G)
    fallback_mu    <- (cutpoints[1] + cutpoints[K - 1]) / 2
    fallback_sigma <- max((cutpoints[K - 1] - cutpoints[1]) / 2, 1e-3)
    for (g in seq_len(G)) {
        counts    <- ngk[g, ]
        n_g       <- sum(counts)
        populated <- which(counts > 0)
        if (length(populated) < 3) {
            mu0[g] <- fallback_mu; sigma0[g] <- fallback_sigma
            next
        }
        j <- populated[1]
        k <- populated[length(populated)]
        c_lo <- cutpoints[j]        # right edge of bin j
        c_hi <- cutpoints[k - 1]    # left edge of bin k
        P1   <- sum(counts[1:j]) / n_g
        Ptop <- sum(counts[k:K]) / n_g
        a <- stats::qnorm(P1)
        b <- stats::qnorm(1 - Ptop)
        if (!is.finite(a) || !is.finite(b) || (b - a) < 1e-6) {
            mu0[g] <- fallback_mu; sigma0[g] <- fallback_sigma
            next
        }
        sigma_g <- (c_hi - c_lo) / (b - a)
        mu_g    <- c_hi - b * sigma_g
        if (!is.finite(sigma_g) || sigma_g <= 0 || !is.finite(mu_g)) {
            mu0[g] <- fallback_mu; sigma0[g] <- fallback_sigma
            next
        }
        mu0[g]    <- mu_g
        sigma0[g] <- sigma_g
    }
    list(mu = mu0, sigma = sigma0)
}

.refine_within <- function(qgk, cuts_for_moments, mug_z, sigmag_z,
                           tol = 1e-4, maxit = 100) {
    ## For each group g, iteratively refine (mug_g, sigmag_g) via
    ## truncated-normal moment matching.
    ##
    ## Convergence is checked per district on the PARAMETERS, not on the
    ## likelihood.  A district is flagged converged once both
    ##
    ##   |Delta mug_g| / sigmag_g < tol   and   |Delta log sigmag_g| < tol
    ##
    ## between successive iterations, and is then skipped in all future
    ## iterations.  The outer loop exits when every district has
    ## converged (or maxit is hit).
    ##
    ## Why parameters rather than the likelihood.  A likelihood-change
    ## rule stops when the objective stops moving, which is exactly what
    ## happens on a flat ridge -- and some districts have very flat
    ## ridges, typically small ones with most of their mass in an
    ## open-ended extreme bin, whose sigma is only weakly identified.
    ## There the likelihood settles while mu and sigma are still
    ## drifting by amounts large enough to matter, so the fit stops
    ## early and silently.  A parameter-change rule cannot stop early
    ## that way; on those districts it runs to maxit instead, which at
    ## least shows up in `converged` rather than passing unnoticed.
    ## It also removes an n-dependence: the old criterion used the
    ## PER-OBSERVATION log-likelihood, so its implied tolerance on a
    ## district's total was n_g * tol, looser for large districts.
    ##
    ## Both terms are scale-free by construction, which matters because
    ## the internal scale is the test-score scale when the caller
    ## supplied cutpoints and an approximately standardized latent scale
    ## when they did not.  Dividing the mean change by that district's
    ## own sigma, and using the change in log sigma, makes `tol` mean
    ## "settled to this fraction of the district's own SD" either way --
    ## so 1e-3 does not silently become a thousandth of a scale-score
    ## point when cutpoints are supplied in score units.
    ##
    ## Fast-converging districts stop being updated once they settle, so
    ## the total work is roughly (sum over g of iterations to converge
    ## district g), not (G * iterations for the slowest district).
    G <- length(mug_z)
    K <- length(cuts_for_moments) + 1
    c_left  <- c(-Inf, cuts_for_moments)
    c_right <- c(cuts_for_moments, Inf)

    converged  <- logical(G)
    iter_count <- integer(G)

    for (iter in 1:maxit) {
        active <- which(!converged)
        if (length(active) == 0) break

        for (g in active) {
            sg <- max(sigmag_z[g], 1e-8)
            z_left  <- (c_left  - mug_z[g]) / sg
            z_right <- (c_right - mug_z[g]) / sg
            phi_l <- dnorm(z_left)
            phi_r <- dnorm(z_right)
            Phi_l <- pnorm(z_left)
            Phi_r <- pnorm(z_right)
            p_k   <- Phi_r - Phi_l
            ## E[Y | bin k]
            bm <- mug_z[g] + sg * (phi_l - phi_r) / p_k
            ## E[Y^2 | bin k]
            zphi_l <- ifelse(is.infinite(z_left),  0, z_left  * phi_l)
            zphi_r <- ifelse(is.infinite(z_right), 0, z_right * phi_r)
            b2 <- mug_z[g]^2 +
                  2 * mug_z[g] * sg * (phi_l - phi_r) / p_k +
                  sg^2 * (1 + (zphi_l - zphi_r) / p_k)
            new_mu <- sum(qgk[g, ] * bm)
            new_v  <- sum(qgk[g, ] * b2) - new_mu^2
            new_sg <- sqrt(max(new_v, 0))

            iter_count[g] <- iter_count[g] + 1L

            ## Scale-free parameter-change test.  The mean change is
            ## measured in units of the district's own SD, and the SD
            ## change on the log scale, so `tol` means the same thing
            ## whether the internal scale is test-score points or the
            ## standardized latent scale.  Guard the denominators: a
            ## district whose sigma has collapsed toward zero would
            ## otherwise divide by (almost) nothing and never converge.
            sg_ref  <- max(sg, new_sg, 1e-8)
            d_mu    <- abs(new_mu - mug_z[g]) / sg_ref
            d_logsg <- if (new_sg > 1e-8 && sg > 1e-8) {
                abs(log(new_sg) - log(sg))
            } else {
                abs(new_sg - sg) / sg_ref
            }
            if (d_mu < tol && d_logsg < tol) {
                converged[g] <- TRUE
            }

            mug_z[g]    <- new_mu
            sigmag_z[g] <- new_sg
        }
    }

    list(mug = mug_z, sigmag = sigmag_z,
         iterations = max(iter_count),
         converged  = all(converged),
         per_district_iterations = iter_count)
}

.refine_mu_only <- function(qgk, cuts_for_moments, mug_z, sigma_fixed,
                            tol = 1e-4, maxit = 100) {
    ## Refine only mu_g for each group, holding sigma fixed at
    ## `sigma_fixed` (a scalar).  Mirrors `.refine_within` but drops
    ## the sigma update; used to salvage districts with fewer than
    ## three populated bins by borrowing sigma from the identified
    ## pool.  Convergence uses the same scale-free parameter-change test
    ## as `.refine_within`: a group is done once |Delta mu| / sigma
    ## falls below `tol`.  There is no sigma term because sigma is held
    ## fixed here rather than estimated.
    G <- length(mug_z)
    K <- length(cuts_for_moments) + 1
    c_left  <- c(-Inf, cuts_for_moments)
    c_right <- c(cuts_for_moments, Inf)
    sg <- max(sigma_fixed, 1e-8)

    converged  <- logical(G)
    iter_count <- integer(G)

    for (iter in 1:maxit) {
        active <- which(!converged)
        if (length(active) == 0) break
        for (g in active) {
            z_left  <- (c_left  - mug_z[g]) / sg
            z_right <- (c_right - mug_z[g]) / sg
            phi_l <- dnorm(z_left)
            phi_r <- dnorm(z_right)
            Phi_l <- pnorm(z_left)
            Phi_r <- pnorm(z_right)
            p_k   <- pmax(Phi_r - Phi_l, 1e-300)
            bm <- mug_z[g] + sg * (phi_l - phi_r) / p_k
            ## Empty bins have qgk = 0 and drop out; use safe products.
            contrib <- qgk[g, ] * bm
            contrib[qgk[g, ] == 0] <- 0
            new_mu <- sum(contrib)

            iter_count[g] <- iter_count[g] + 1L
            if (abs(new_mu - mug_z[g]) / sg < tol) converged[g] <- TRUE
            mug_z[g]  <- new_mu
        }
    }

    list(mug = mug_z,
         iterations = max(iter_count, 0L),
         converged  = all(converged),
         per_district_iterations = iter_count)
}

.refine_cutpoints <- function(mug_z, sigmag_z, pg, pooled_cum, init_cuts,
                              tol = 1e-3, maxit = 100) {
    ## Update cutpoints by matching pooled bin proportions to the implied
    ## mixture-of-normals CDF.  For each target cumulative proportion
    ## P_k in pooled_cum, find c such that
    ##    sum_g pg * pnorm((c - mug_z[g]) / sigmag_z[g]) = P_k
    K <- length(init_cuts) + 1
    cuts <- init_cuts
    for (iter in 1:maxit) {
        new_cuts <- numeric(K - 1)
        for (k in 1:(K - 1)) {
            target <- pooled_cum[k]
            f <- function(c) sum(pg * pnorm((c - mug_z) / pmax(sigmag_z, 1e-8))) - target
            ## Bracket the root.  Start broad: any cutpoint should be
            ## within a few standard deviations of the standardized mean.
            lo <- min(mug_z) - 10
            hi <- max(mug_z) + 10
            new_cuts[k] <- stats::uniroot(f, c(lo, hi), tol = tol)$root
        }
        if (max(abs(new_cuts - cuts)) < tol) {
            cuts <- new_cuts
            attr(cuts, "iterations") <- iter
            attr(cuts, "converged")  <- TRUE
            return(cuts)
        }
        cuts <- new_cuts
    }
    attr(cuts, "iterations") <- maxit
    attr(cuts, "converged")  <- FALSE
    cuts
}

.fisher_s2_log_sigma <- function(mug, sigmag, ng, cuts_std) {
    ## Per-district sampling variance of log(sigma_g) under the
    ## cutpoint-known normal model.  Uses the (2,2) entry of the
    ## inverse Fisher information per observation, divided by n_g
    ## and by sigma_g^2 (delta-method on the log scale).
    G <- length(mug)
    s2 <- numeric(G)
    for (g in 1:G) {
        sg <- max(sigmag[g], 1e-6)
        z  <- c(-Inf, (cuts_std - mug[g]) / sg, Inf)
        Phi <- pnorm(z)
        phi <- dnorm(z)
        p_k <- diff(Phi)
        p_k[p_k < 1e-12] <- 1e-12
        zphi <- ifelse(is.infinite(z), 0, z * phi)
        a_k <- -diff(phi)  / sg
        b_k <-  diff(zphi) / sg  # note sign: derivative of p_k wrt sigma_g
        b_k <- -b_k
        I11 <- sum(a_k^2 / p_k)
        I12 <- sum(a_k * b_k / p_k)
        I22 <- sum(b_k^2 / p_k)
        det_I <- I11 * I22 - I12^2
        if (det_I <= 0 || I22 <= 0) {
            s2[g] <- 1 / (2 * ng[g])  # fall back to continuous-data
        } else {
            ## var(hat sigma_g) per observation = I11 / det_I
            var_sig_per_obs <- I11 / det_I
            s2[g] <- var_sig_per_obs / (ng[g] * sg^2)  # delta to log
        }
    }
    s2
}

.eb_shrink <- function(log_sig_hat, n_g, mug, sigmag, cuts_std) {
    ## Fay-Herriot empirical-Bayes shrinkage of per-district log SDs.
    ## Sampling variance s_g^2 comes from the Fisher information for the
    ## binned-normal model with known cutpoints.  NA entries in
    ## log_sig_hat are treated as unidentified and passed through to the
    ## output; they are excluded from tau^2 / overall-mean estimation.
    G <- length(log_sig_hat)
    ok <- !is.na(log_sig_hat)
    s2 <- rep(NA_real_, G)
    s2[ok] <- .fisher_s2_log_sigma(mug[ok], sigmag[ok], n_g[ok], cuts_std)
    wts_ok <- n_g[ok] / sum(n_g[ok])
    overall_mean <- sum(wts_ok * log_sig_hat[ok])
    var_obs <- sum(wts_ok * (log_sig_hat[ok] - overall_mean)^2)
    tau2 <- max(var_obs - sum(wts_ok * s2[ok]), 0)
    w_g <- rep(NA_real_, G)
    w_g[ok] <- tau2 / (tau2 + s2[ok])
    log_sig_shrunk <- rep(NA_real_, G)
    log_sig_shrunk[ok] <- w_g[ok] * log_sig_hat[ok] + (1 - w_g[ok]) * overall_mean
    list(log_sig = log_sig_shrunk, w = w_g, tau2 = tau2,
         overall_log_mean = overall_mean, s2 = s2)
}

.fisher_s2_mu <- function(mug, sigmag, ng, cuts_std) {
    ## Per-district sampling variance of mug under the cutpoint-known
    ## normal model.  Uses the (1,1) entry of the inverse Fisher
    ## information per observation, divided by n_g.  No delta method
    ## needed because mu is already on the natural scale.
    G <- length(mug)
    s2 <- numeric(G)
    for (g in 1:G) {
        sg <- max(sigmag[g], 1e-6)
        z  <- c(-Inf, (cuts_std - mug[g]) / sg, Inf)
        Phi <- pnorm(z)
        phi <- dnorm(z)
        p_k <- diff(Phi)
        p_k[p_k < 1e-12] <- 1e-12
        zphi <- ifelse(is.infinite(z), 0, z * phi)
        a_k <- -diff(phi)  / sg          # d p_k / d mu
        b_k <-  diff(zphi) / sg
        b_k <- -b_k                       # d p_k / d sigma
        I11 <- sum(a_k^2 / p_k)
        I12 <- sum(a_k * b_k / p_k)
        I22 <- sum(b_k^2 / p_k)
        det_I <- I11 * I22 - I12^2
        if (det_I <= 0 || I11 <= 0) {
            s2[g] <- sg^2 / ng[g]         # fall back to continuous-data
        } else {
            ## var(hat mu_g) per observation = I22 / det_I
            var_mu_per_obs <- I22 / det_I
            s2[g] <- var_mu_per_obs / ng[g]
        }
    }
    s2
}

.fisher_s2_sandwich <- function(mug, sigmag, ng, cuts_std, ngk_sub) {
    ## Per-district Huber-White SANDWICH variance of (mug_g, sigma_g)
    ## under the cutpoint-known normal model, for the sigma-unknown
    ## districts (>= 3 populated bins).  Uses observed bin counts in
    ## the "meat" J, model probabilities in the "bread" I:
    ##
    ##   J11 = sum_k n_{g,k} * (a_k / p_k)^2
    ##   J12 = sum_k n_{g,k} * (a_k * b_k) / p_k^2
    ##   J22 = sum_k n_{g,k} * (b_k / p_k)^2
    ##
    ##   I_matrix = n_g * [I11 I12; I12 I22]   (Fisher info for district)
    ##   V_sandwich = I_matrix^{-1} %*% J_matrix %*% I_matrix^{-1}
    ##
    ##   s2_mu_sandwich     = V_sandwich[1,1]
    ##   s2_logsig_sandwich = V_sandwich[2,2] / sigma_g^2   (delta method)
    ##
    ## Sanity check: under a correctly specified within-group normal,
    ## E[n_{g,k}] = n_g * p_k, so E[J11] = n_g * I11 (and analogously
    ## for J12, J22), i.e. E[J_matrix] = n_g * I_matrix = I_matrix / 1.
    ## Then V_sandwich = I_matrix^{-1} * (n_g * I) * I_matrix^{-1}
    ## = (1/n_g) * I^{-1} = I_matrix^{-1}, the Fisher-info variance.
    ## Under misspecification, J diverges from n_g * I and the sandwich
    ## picks up the discrepancy -- yielding an SE calibrated even when
    ## the within-group scores are not exactly normal.
    ##
    ## Returns NA for a district when the Fisher information is
    ## non-invertible or the resulting diagonal is non-positive.
    G <- length(mug)
    s2_mu_sw     <- rep(NA_real_, G)
    s2_logsig_sw <- rep(NA_real_, G)
    for (g in seq_len(G)) {
        sg <- max(sigmag[g], 1e-6)
        z  <- c(-Inf, (cuts_std - mug[g]) / sg, Inf)
        Phi <- pnorm(z)
        phi <- dnorm(z)
        p_k <- diff(Phi)
        p_k[p_k < 1e-12] <- 1e-12
        zphi <- ifelse(is.infinite(z), 0, z * phi)
        a_k <- -diff(phi)  / sg
        b_k <-  diff(zphi) / sg
        b_k <- -b_k
        I11 <- sum(a_k^2 / p_k)
        I12 <- sum(a_k * b_k / p_k)
        I22 <- sum(b_k^2 / p_k)
        counts <- ngk_sub[g, ]
        J11 <- sum(counts * (a_k / p_k)^2)
        J12 <- sum(counts * (a_k * b_k) / p_k^2)
        J22 <- sum(counts * (b_k / p_k)^2)
        n_g <- ng[g]
        Imat <- n_g * matrix(c(I11, I12, I12, I22), nrow = 2)
        Jmat <-       matrix(c(J11, J12, J12, J22), nrow = 2)
        Iinv <- tryCatch(solve(Imat), error = function(e) NULL)
        if (is.null(Iinv)) next
        V <- Iinv %*% Jmat %*% Iinv
        v_mu  <- V[1, 1]
        v_sig <- V[2, 2]
        if (is.finite(v_mu)  && v_mu  > 0) s2_mu_sw[g]     <- v_mu
        if (is.finite(v_sig) && v_sig > 0) s2_logsig_sw[g] <- v_sig / sg^2
    }
    list(s2_mu_sandwich = s2_mu_sw, s2_logsig_sandwich = s2_logsig_sw)
}

.fisher_s2_mu_known_sigma_sandwich <- function(mug, sigmag, ng, cuts_std,
                                               ngk_sub) {
    ## Sandwich variance of mu for salvaged districts with sigma held
    ## fixed at the pool-borrowed value (two-bin and one-interior-bin
    ## salvages).  V_sandwich_mu = J11 / (n_g * I11)^2.  Returns NA when
    ## I11 is non-positive / non-finite.
    G <- length(mug)
    s2_mu_sw <- rep(NA_real_, G)
    for (g in seq_len(G)) {
        sg <- max(sigmag[g], 1e-6)
        z  <- c(-Inf, (cuts_std - mug[g]) / sg, Inf)
        Phi <- pnorm(z)
        phi <- dnorm(z)
        p_k <- diff(Phi)
        p_k[p_k < 1e-12] <- 1e-12
        a_k <- -diff(phi) / sg
        I11 <- sum(a_k^2 / p_k)
        counts <- ngk_sub[g, ]
        J11 <- sum(counts * (a_k / p_k)^2)
        n_g <- ng[g]
        if (!is.finite(I11) || I11 <= 0) next
        val <- J11 / (n_g * I11)^2
        if (is.finite(val) && val > 0) s2_mu_sw[g] <- val
    }
    s2_mu_sw
}

.fisher_s2_mu_known_sigma <- function(mug, sigmag, ng, cuts_std) {
    ## Per-district sampling variance of mug under the cutpoint-known
    ## normal model with SIGMA KNOWN (not estimated).  Uses 1/I11 --
    ## the reciprocal of the (1,1) entry of the Fisher information per
    ## observation -- since with sigma held fixed we do NOT marginalize
    ## the Fisher information over sigma estimation uncertainty.  This
    ## is the appropriate SE for salvaged districts (< 3 populated
    ## bins) that use the pool-borrowed sigma, since for them mu is the
    ## only quantity actually estimated.
    G <- length(mug)
    s2 <- numeric(G)
    for (g in 1:G) {
        sg <- max(sigmag[g], 1e-6)
        z  <- c(-Inf, (cuts_std - mug[g]) / sg, Inf)
        Phi <- pnorm(z)
        phi <- dnorm(z)
        p_k <- diff(Phi)
        p_k[p_k < 1e-12] <- 1e-12
        a_k <- -diff(phi) / sg
        I11 <- sum(a_k^2 / p_k)
        if (!is.finite(I11) || I11 <= 0) {
            s2[g] <- NA_real_
        } else {
            s2[g] <- 1 / (ng[g] * I11)
        }
    }
    s2
}

.eb_shrink_mu <- function(mu_hat, n_g, mug, sigmag, cuts_std) {
    ## Fay-Herriot empirical-Bayes shrinkage of per-district means.
    ## Sampling variance s_g^2 comes from the Fisher information for the
    ## binned-normal model with known cutpoints.
    s2 <- .fisher_s2_mu(mug, sigmag, n_g, cuts_std)
    wts <- n_g / sum(n_g)
    overall_mean <- sum(wts * mu_hat)
    var_obs <- sum(wts * (mu_hat - overall_mean)^2)
    tau2 <- max(var_obs - sum(wts * s2), 0)
    w_g <- tau2 / (tau2 + s2)
    mu_shrunk <- w_g * mu_hat + (1 - w_g) * overall_mean
    list(mu = mu_shrunk, w = w_g, tau2 = tau2,
         overall_mean = overall_mean, s2 = s2)
}

.map_shrink_mu <- function(ngk, mu_hat, sigmag, cuts_std) {
    ## Per-district MAP estimate of mu under the actual multinomial bin
    ## likelihood with fixed sigma and a Gaussian prior on mu with mean
    ## bar_mu and variance tau^2 estimated from the unshrunk mu_hat.
    ##
    ## tau^2 is estimated by method of moments: tau^2 = max(0, Var(mu_hat) -
    ## mean(V_g)) where V_g is the Fisher information variance of mu_hat.
    ## NA entries in mu_hat / sigmag indicate unidentified districts; they
    ## are excluded from tau^2 estimation and passed through as NA.
    G <- length(mu_hat)
    n_g <- rowSums(ngk)
    ok <- !is.na(mu_hat) & !is.na(sigmag)

    s2 <- rep(NA_real_, G)
    s2[ok] <- .fisher_s2_mu(mu_hat[ok], sigmag[ok], n_g[ok], cuts_std)
    wts_ok <- n_g[ok] / sum(n_g[ok])
    overall_mean <- sum(wts_ok * mu_hat[ok])
    var_obs <- sum(wts_ok * (mu_hat[ok] - overall_mean)^2)
    tau2 <- max(var_obs - sum(wts_ok * s2[ok]), 0)

    mu_map <- rep(NA_real_, G)

    if (tau2 <= 0) {
        mu_map[ok] <- overall_mean
        return(list(mu = mu_map, tau2 = tau2,
                    overall_mean = overall_mean, s2 = s2))
    }

    ## Negative log posterior: bin-multinomial NLL + Gaussian prior penalty
    for (g in which(ok)) {
        sg <- max(sigmag[g], 1e-6)
        counts <- ngk[g, ]
        neg_lp <- function(mu) {
            z <- c(-Inf, (cuts_std - mu) / sg, Inf)
            p <- diff(pnorm(z))
            p[p < 1e-12] <- 1e-12
            -sum(counts * log(p)) + (mu - overall_mean)^2 / (2 * tau2)
        }
        lo <- min(overall_mean, mu_hat[g]) - 6 * sqrt(s2[g])
        hi <- max(overall_mean, mu_hat[g]) + 6 * sqrt(s2[g])
        opt <- stats::optimize(neg_lp, lower = lo, upper = hi, tol = 1e-3)
        mu_map[g] <- opt$minimum
    }

    list(mu = mu_map, tau2 = tau2, overall_mean = overall_mean, s2 = s2)
}

fast_hetop <- function(ngk, cutpoints_known = FALSE, cutpoints = NULL,
                     pooled_mean = 0, pooled_sd = 1,
                     scope,
                     estimator = "ML",
                     tol = 1e-4, maxit = 100,
                     conf.level = 0.95,
                     estimate_unidentified_districts = TRUE) {
    ## `scope` is required and deliberately has no default: which
    ## standard errors are appropriate depends on whether the units
    ## behind each group's counts are a sample or the whole population,
    ## and that is a fact about the data that only the caller knows.
    ## Guessing it would silently return the wrong uncertainty, so we
    ## refuse to proceed without it.
    if (missing(scope)) {
        stop("`scope` is required and has no default. Set scope = \"sample\" ",
             "if the units behind each group's counts are a sample from a ",
             "larger population, so the standard errors should reflect both ",
             "sampling and binning error; set scope = \"population\" if the ",
             "counts cover every unit in the group (population or ",
             "administrative data), so the only uncertainty comes from ",
             "seeing bin counts rather than individual scores.",
             call. = FALSE)
    }
    scope <- match.arg(scope, c("sample", "population"))

    ## Standard errors come from the model-implied Fisher information.
    ##
    ## A Huber-White sandwich alternative is implemented below and was
    ## once exposed as a `robust` argument, but it is no longer offered.
    ## A sandwich corrects the variance, not the estimand: where the
    ## within-group normal is wrong the damage is bias in the fitted
    ## mean, and widening an interval to absorb bias is not what a
    ## confidence interval is for.  The code is kept, pinned off here,
    ## so the alternative stays available to anyone studying it.
    robust <- FALSE

    ## Match the user-facing `estimator` argument to a boolean the rest
    ## of the function uses; the internal name `eb_shrink` is kept for
    ## minimal disruption to the plumbing below.
    estimator <- match.arg(estimator, c("ML", "EB_shrunk"))
    eb_shrink <- (estimator == "EB_shrunk")

    ## Empirical Bayes is not offered for population data.
    ##
    ## Fay-Herriot shrinkage assumes mu_hat_g is drawn from a sampling
    ## distribution centred on mu_g with a known variance V_g, and
    ## weights each district by tau^2 / (tau^2 + V_g).  When the counts
    ## cover every student, there is no sampling distribution: the
    ## counts are fixed, mu_hat_g is a deterministic function of them,
    ## and the only uncertainty left is that the bins do not pin mu_g
    ## down exactly.  That uncertainty is real, but it is not a
    ## sampling variance, and the variance we could put in its place --
    ## the total Fisher variance minus sigma^2/n -- is obtained by
    ## subtraction rather than derived as the posterior variance of a
    ## finite-population mean given its bin counts.  Shrinking by a
    ## weight we cannot derive would give the analyst a number with no
    ## clear justification, so we decline instead.
    if (eb_shrink && scope == "population") {
        stop("estimator = \"EB_shrunk\" is not available with ",
             "scope = \"population\".\n",
             "  Empirical Bayes shrinkage is derived for estimates that ",
             "carry sampling error.\n",
             "  With population counts there is no sampling ",
             "distribution to shrink against, and\n",
             "  the shrinkage weights would rest on an approximation ",
             "we cannot justify.\n",
             "  Use estimator = \"ML\", or scope = \"sample\" if the ",
             "counts really are a sample.",
             call. = FALSE)
    }

    ## iterate defaults to TRUE when cutpoints are estimated from the
    ## data (cutpoints = NULL) and FALSE when cutpoints are supplied,
    ## since iteration only refines data-derived cutpoints and is not
    ## meaningful (and is rejected below) when cutpoints are already
    ## known. Iterating is cheap, so it is on by default whenever it
    ## applies, even though it did not meaningfully change the
    ## per-district estimates in either the Texas data or the
    ## simulation; users who want the faster non-iterative estimate can
    ## still pass iterate = FALSE explicitly.

    if (!is.numeric(ngk))  stop("ngk must be a GxK numeric matrix of category counts")
    if (!is.matrix(ngk))   stop("ngk must be a GxK numeric matrix of category counts")
    if (any(is.na(ngk)))   stop("ngk cannot contain missing values")
    if (any(ngk < 0))      stop("ngk must contain only non-negative values")
    if (any(rowSums(ngk) <= 0)) stop("ngk contains at least one row with insufficient data")

    G <- nrow(ngk)
    K <- ncol(ngk)
    if (K <= 2) stop("Function requires K >= 3 categories")
    if (!is.numeric(conf.level) || length(conf.level) != 1 ||
        conf.level <= 0 || conf.level >= 1) {
        stop("conf.level must be a single number strictly between 0 and 1")
    }
    z_crit <- qnorm(1 - (1 - conf.level) / 2)

    ## `cutpoints_known` states which of the two modes the caller is in, and
    ## the other cut-related arguments are validated against it rather
    ## than being inferred from whether `cutpoints` happens to be NULL.
    ## Contradictory combinations are errors, not silently resolved.
    if (!isTRUE(cutpoints_known) && !isFALSE(cutpoints_known)) {
        stop("`cutpoints_known` must be TRUE or FALSE.", call. = FALSE)
    }

    if (cutpoints_known) {
        if (is.null(cutpoints)) {
            stop("cutpoints_known = TRUE requires `cutpoints`: supply the K-1 = ",
                 K - 1, " cut scores on the test-score scale.", call. = FALSE)
        }
        if (!identical(pooled_mean, 0) || !identical(pooled_sd, 1)) {
            stop("`pooled_mean` and `pooled_sd` apply only when ",
                 "cutpoints_known = FALSE. With known cutpoints the scale is ",
                 "already fixed by the cut scores themselves.",
                 call. = FALSE)
        }
    } else {
        if (!is.null(cutpoints)) {
            stop("`cutpoints` were supplied but cutpoints_known = FALSE. Set ",
                 "cutpoints_known = TRUE to use them, or drop them to have the ",
                 "cutpoints estimated from the pooled bin proportions.",
                 call. = FALSE)
        }
        if (!is.numeric(pooled_mean) || length(pooled_mean) != 1 ||
            !is.finite(pooled_mean)) {
            stop("`pooled_mean` must be a single finite number.",
                 call. = FALSE)
        }
        if (!is.numeric(pooled_sd) || length(pooled_sd) != 1 ||
            !is.finite(pooled_sd) || pooled_sd <= 0) {
            stop("`pooled_sd` must be a single positive number.",
                 call. = FALSE)
        }
    }

    ## Whether the outer cutpoint loop runs is not a user choice: it is
    ## determined by whether there are cutpoints to refine.  When the
    ## caller supplies them there is nothing to iterate toward, and when
    ## they are derived from the pooled proportions the refinement is
    ## always wanted.  Exposing this as an argument only created an
    ## invalid combination to guard against.
    iterate <- !cutpoints_known

    ## Tolerance for the OUTER cutpoint-refinement EM, derived from `tol`
    ## rather than exposed as an argument.
    ##
    ## It has to be looser than the per-group tolerance.  The cutpoint
    ## solve is a function of the current group fits, which are converged
    ## only to `tol`, so their residual wobble reaches the cutpoints at
    ## roughly the same order.  Ask the EM for precision at or below that
    ## floor and the cutpoints jitter rather than settle: the criterion
    ## can never be met and the loop simply exhausts `maxit`.  Setting
    ## the two equal is already marginal, so the ratio below leaves an
    ## order of magnitude of headroom.
    ##
    ## That relationship is a property of the algorithm, not a modelling
    ## choice, so there is nothing for a caller to decide and no way for
    ## them to get it wrong.  `tol` is the single knob; this follows.
    tol_cutpoints_std <- 10 * tol

    pooled_props <- colSums(ngk) / sum(ngk)
    pooled_cum   <- cumsum(pooled_props)

    if (!cutpoints_known) {
        ## Cutpoints are derived from the pooled cumulative proportions
        ## via qnorm(), which returns +/-Inf for a 0 or 1 argument.  The
        ## check therefore only matters in this branch.  When the caller
        ## supplies known cutpoints on the native scale, individual
        ## districts may legitimately have empty bins (e.g. tail bins in
        ## small districts) and are still fittable from their populated
        ## bins alone.
        if (any(pooled_props == 0)) {
            stop("At least one pooled bin has zero count; fast_hetop cannot derive cutpoints from pooled proportions when a bin is empty.  Set cutpoints_known = TRUE and supply them via the `cutpoints` argument.")
        }
        scale_known <- FALSE
        ## Starting cuts for the outer alternation: treat the pooled
        ## distribution as normal and read the cuts off its quantiles.
        ## This is close to the solution whenever the pooled mixture is
        ## near normal, which is usual, since a mixture of normals whose
        ## means are themselves near-normal is very nearly normal.  The
        ## fixed point does not depend on this choice: restarting the
        ## alternation from cuts shifted, stretched, compressed, or set
        ## to a crude equally-spaced guess reaches the same cuts to a
        ## tenth of a score point on the Texas data, differing only in
        ## how many passes it takes to get there.
        cutpoints   <- qnorm(pooled_cum[1:(K-1)])
    } else {
        scale_known <- TRUE
        if (!is.numeric(cutpoints)) {
            stop("`cutpoints` must be a numeric vector.", call. = FALSE)
        }
        if (length(cutpoints) != (K - 1)) {
            stop("`cutpoints` has length ", length(cutpoints), " but `ngk` has ",
                 K, " bins, so K-1 = ", K - 1, " cut scores are required.",
                 call. = FALSE)
        }
        if (any(is.na(cutpoints))) {
            stop("`cutpoints` cannot contain missing values.", call. = FALSE)
        }
        if (any(diff(cutpoints) <= 0)) {
            stop("`cutpoints` must be strictly increasing.", call. = FALSE)
        }
    }

    ## Starting values for each district's ML search.
    ## `.district_3bin_start` seeds each district's iteration from a
    ## closed-form solution computed entirely from that district's own
    ## bin counts, using the current cutpoints (whether they were
    ## supplied by the caller or inferred from the pooled bin
    ## proportions).  No pooled reference is used.
    start <- .district_3bin_start(ngk, cutpoints)
    fit <- list(mug_z = start$mu, sigmag_z = start$sigma,
                pg = rowSums(ngk) / sum(ngk),
                pooled_mean = 0, pooled_sd = 1)

    iter_info <- list(iterated = FALSE, eb_shrink = FALSE)

    ## Districts with fewer than three populated bins are not
    ## jointly identified under a 2-parameter (mu_g, sigma_g) normal
    ## model, so they are excluded from any iterative update and set
    ## to NA at the end.
    ng_v         <- rowSums(ngk)
    qgk_all      <- ngk / ng_v
    bins_used    <- rowSums(ngk > 0)
    unidentified <- bins_used < 3
    identified   <- !unidentified

    ## EM-style refinement when cutpoints are unknown.  Alternates
    ## between (a) refitting cutpoints to match the pooled bin
    ## proportions under the implied mixture-of-normals CDF, using the
    ## current per-district fits, and (b) refitting per-district
    ## (mu_g, sigma_g) via the truncated-normal MLE at the new
    ## cutpoints.  Iterates until cutpoints stop changing by more than
    ## `tol_cutpoints_std`.
    if (iterate) {
        em_converged  <- FALSE
        em_iterations <- 0
        within_iter_total <- 0
        within_conv <- TRUE
        ## Per-district iteration counts, accumulated across EM passes:
        ## each pass refits every identified district, so a district's
        ## total work is the sum over passes, not the count from the
        ## last one.  Kept so callers can see the spread rather than
        ## just the maximum.
        per_district_iters <- rep(NA_integer_, G)
        per_district_iters[identified] <- 0L
        for (em_iter in seq_len(maxit)) {
            old_cuts <- cutpoints
            ## Step (a): update cutpoints
            new_cuts <- .refine_cutpoints(fit$mug_z, fit$sigmag_z, fit$pg,
                                          pooled_cum, init_cuts = cutpoints,
                                          tol = tol_cutpoints_std, maxit = maxit)
            attributes(new_cuts) <- NULL
            cutpoints <- new_cuts
            ## Step (b): update per-district (mu, sigma) via full MLE,
            ## on identified districts only.
            if (any(identified)) {
                ref <- .refine_within(qgk_all[identified, , drop = FALSE], cutpoints,
                                      fit$mug_z[identified], fit$sigmag_z[identified],
                                      tol = tol, maxit = maxit)
                fit$mug_z[identified]    <- ref$mug
                fit$sigmag_z[identified] <- ref$sigmag
                within_iter_total <- within_iter_total + ref$iterations
                within_conv <- within_conv && ref$converged
                per_district_iters[identified] <-
                    per_district_iters[identified] + ref$per_district_iterations
            }
            em_iterations <- em_iter
            ## Outer EM stopping rule: a change in CUTPOINT LOCATION, so
            ## it takes tol_cutpoints_std, not the per-district likelihood tol.
            if (max(abs(cutpoints - old_cuts), na.rm = TRUE) < tol_cutpoints_std) {
                em_converged <- TRUE
                break
            }
        }
        iter_info <- list(
            iterated   = TRUE,
            iterations = em_iterations,
            converged  = em_converged,
            within_iterations = within_iter_total,
            within_converged  = within_conv,
            within_iterations_per_district = per_district_iters
        )
    }

    ## Per-district (mu_g, sigma_g) MLE.  When iterate = TRUE the EM
    ## loop above already polished each district's (mu, sigma) at the
    ## converged cutpoints; when iterate = FALSE we do the polish here.
    if (!iterate) {
        if (scale_known) {
            ## The closed-form 3-bin solution is the exact MLE only
            ## when the data set itself has exactly K = 3 bins -- not
            ## merely when a given district happens to have exactly 3
            ## *populated* bins out of a larger K.
            exact3     <- (K == 3) & (bins_used == K) & identified
            needs_iter <- identified & !exact3

            if (any(needs_iter)) {
                ref <- .refine_within(qgk_all[needs_iter, , drop = FALSE], cutpoints,
                                      fit$mug_z[needs_iter], fit$sigmag_z[needs_iter],
                                      tol = tol, maxit = maxit)
                fit$mug_z[needs_iter]    <- ref$mug
                fit$sigmag_z[needs_iter] <- ref$sigmag
                iter_info$within_iterations <- ref$iterations
                iter_info$within_converged  <- ref$converged
                pdi <- rep(NA_integer_, G)
                pdi[needs_iter] <- ref$per_district_iterations
                ## Districts solved in closed form did no iterating.
                pdi[exact3] <- 0L
                iter_info$within_iterations_per_district <- pdi
            } else {
                iter_info$within_iterations <- 0
                iter_info$within_converged  <- TRUE
                pdi <- rep(NA_integer_, G)
                pdi[identified] <- 0L
                iter_info$within_iterations_per_district <- pdi
            }
            iter_info$within_exact3 <- which(exact3)
        } else {
            if (any(identified)) {
                ref <- .refine_within(qgk_all[identified, , drop = FALSE], cutpoints,
                                      fit$mug_z[identified], fit$sigmag_z[identified],
                                      tol = tol, maxit = maxit)
                fit$mug_z[identified]    <- ref$mug
                fit$sigmag_z[identified] <- ref$sigmag
                iter_info$within_iterations <- ref$iterations
                iter_info$within_converged  <- ref$converged
                pdi <- rep(NA_integer_, G)
                pdi[identified] <- ref$per_district_iterations
                iter_info$within_iterations_per_district <- pdi
            } else {
                iter_info$within_iterations <- 0
                iter_info$within_converged  <- TRUE
                pdi <- rep(NA_integer_, G)
                pdi[identified] <- 0L
                iter_info$within_iterations_per_district <- pdi
            }
        }
    }

    ## Salvage districts with fewer than three populated bins by
    ## borrowing sigma from the identified pool.  Log-sigma is pooled
    ## unweighted across identified districts (matching HETOP's
    ## `mle_hetop`), then mu is estimated with sigma held fixed:
    ##   - two populated bins: fixed-point iteration on mu alone;
    ##   - one interior populated bin (1 < j < K):
    ##       mu = (cutpoints[j-1] + cutpoints[j]) / 2, the MLE with
    ##       sigma fixed (does not actually depend on sigma);
    ##   - one populated bin at the top (bin K):
    ##       mu = max(mu) over districts already fit (identified plus
    ##       already-salvaged interior/two-bin), since the MLE diverges;
    ##   - one populated bin at the bottom (bin 1): symmetric min.
    ## Salvaged districts flow through the shrinkage / SE / GoF paths
    ## unchanged.  If no district is identified (nothing to pool),
    ## salvage is skipped with a warning and the unidentified rows
    ## remain NA.
    still_unidentified <- unidentified
    salvage_counts <- c(two_bin = 0L, one_interior = 0L,
                        one_top = 0L, one_bottom = 0L)
    sigma_pool     <- NA_real_
    ## Skip the salvage block entirely if the caller doesn't want
    ## unidentified districts salvaged.  In that case the NA fill below
    ## does exactly what the pre-salvage 0.2-2 release did: return NA
    ## for every district with fewer than 3 populated bins.
    if (any(unidentified) && estimate_unidentified_districts) {
        if (!any(identified)) {
            warning("No districts have >= 3 populated bins; ",
                    "cannot pool sigma to salvage unidentified districts. ",
                    "Returning NA for those districts.")
        } else {
            log_sig_pool <- mean(log(pmax(fit$sigmag_z[identified],
                                          .Machine$double.eps)))
            sigma_pool   <- exp(log_sig_pool)
            fit$sigmag_z[unidentified] <- sigma_pool

            uid <- which(unidentified)
            ## Classify by populated-bin pattern.
            pop_list <- lapply(uid, function(g) which(ngk[g, ] > 0))
            n_pop    <- lengths(pop_list)

            two_bin       <- uid[n_pop == 2]
            one_bin_group <- uid[n_pop == 1]
            one_top       <- one_bin_group[
                vapply(pop_list[n_pop == 1],
                       function(w) length(w) == 1 && w == K,
                       logical(1))]
            one_bottom    <- one_bin_group[
                vapply(pop_list[n_pop == 1],
                       function(w) length(w) == 1 && w == 1L,
                       logical(1))]
            one_interior  <- setdiff(one_bin_group, c(one_top, one_bottom))

            ## Rule 2: two populated bins.  Initialize mu at the pool
            ## mean of identified districts' mus (analogous to how we
            ## pool sigma).  The fixed-point iteration converges to the
            ## MLE from any reasonable start; the pool mean is a
            ## defensible prior-informed choice and consistent with the
            ## sigma pooling.
            if (length(two_bin) > 0) {
                mu_pool <- mean(fit$mug_z[identified])
                mu_init <- rep(mu_pool, length(two_bin))
                ref_mu <- .refine_mu_only(
                    qgk_all[two_bin, , drop = FALSE],
                    cutpoints, mu_init, sigma_pool,
                    tol = tol, maxit = maxit)
                fit$mug_z[two_bin] <- ref_mu$mug
                salvage_counts["two_bin"] <- length(two_bin)
                iter_info$mu_pool <- mu_pool
            }

            ## Rule 3: one populated interior bin.
            if (length(one_interior) > 0) {
                mu_int <- vapply(one_interior, function(g) {
                    j <- which(ngk[g, ] > 0)
                    (cutpoints[j - 1] + cutpoints[j]) / 2
                }, numeric(1))
                fit$mug_z[one_interior] <- mu_int
                salvage_counts["one_interior"] <- length(one_interior)
            }

            ## Rules 4 & 5: top / bottom singletons -- use max/min over
            ## districts with a finite mu at this point (identified +
            ## just-salvaged).
            fit_ok <- is.finite(fit$mug_z)
            fit_ok[one_top]    <- FALSE
            fit_ok[one_bottom] <- FALSE
            if (length(one_top) > 0) {
                if (any(fit_ok)) {
                    fit$mug_z[one_top] <- max(fit$mug_z[fit_ok])
                    salvage_counts["one_top"] <- length(one_top)
                } else {
                    fit$mug_z[one_top]    <- NA_real_
                    fit$sigmag_z[one_top] <- NA_real_
                }
            }
            if (length(one_bottom) > 0) {
                if (any(fit_ok)) {
                    fit$mug_z[one_bottom] <- min(fit$mug_z[fit_ok])
                    salvage_counts["one_bottom"] <- length(one_bottom)
                } else {
                    fit$mug_z[one_bottom]    <- NA_real_
                    fit$sigmag_z[one_bottom] <- NA_real_
                }
            }

            ## Any still-NA rows (shouldn't happen unless the fit_ok
            ## fallback above tripped) will be re-NA'd by the safety
            ## net below.
            still_unidentified <- is.na(fit$mug_z) | is.na(fit$sigmag_z)
        }
    }

    iter_info$salvage_counts <- salvage_counts
    iter_info$sigma_pool     <- sigma_pool
    ## Per-district salvage indices, so downstream SE code can apply
    ## sigma-known Fisher information (for two-bin and one-interior-bin
    ## salvages) or NA (for one-extreme-bin salvages) instead of the
    ## degenerate sigma-unknown Fisher information at the imputed sigma.
    if (!exists("two_bin",      inherits = FALSE)) two_bin      <- integer(0)
    if (!exists("one_interior", inherits = FALSE)) one_interior <- integer(0)
    if (!exists("one_top",      inherits = FALSE)) one_top      <- integer(0)
    if (!exists("one_bottom",   inherits = FALSE)) one_bottom   <- integer(0)
    iter_info$salvaged_two_bin      <- two_bin
    iter_info$salvaged_one_interior <- one_interior
    iter_info$salvaged_one_top      <- one_top
    iter_info$salvaged_one_bottom   <- one_bottom

    fit$mug_z[still_unidentified]    <- NA_real_
    fit$sigmag_z[still_unidentified] <- NA_real_
    fit$mug    <- fit$mug_z
    fit$sigmag <- fit$sigmag_z
    iter_info$within_unidentified <- which(still_unidentified)

    ## Optional empirical-Bayes shrinkage of the per-district SDs and means.
    ##
    ## Regardless of eb_shrink, this block also computes a standard error
    ## (on the standardized "z" search scale) for whichever point estimate
    ## is actually returned:
    ##   - eb_shrink = FALSE: the Fisher-information sampling SE of the
    ##     unshrunk per-district MLE.
    ##   - eb_shrink = TRUE: the approximate empirical-Bayes posterior SE
    ##     sqrt(w_g * s2_g), where w_g is the shrinkage weight and s2_g is
    ##     the raw (unshrunk) Fisher sampling variance. This equals
    ##     tau2*s2/(tau2+s2), the normal-normal conjugate posterior
    ##     variance, and ignores uncertainty in the estimated tau2 itself
    ##     -- the standard EB simplification (the same one implicit in
    ##     not propagating hyperparameter uncertainty in mle_hetop() or
    ##     the Efron-prior hyperparameters in fh_hetop()).
    ng_v <- rowSums(ngk)
    ## Fisher-SE and EB shrinkage compute on whatever scale (mu, sigma)
    ## and cutpoints share -- the test-score scale when cutpoints were
    ## supplied, or the EM-converged scale otherwise.  Standardization
    ## for est_std happens once at the end.
    ok_se <- !is.na(fit$mug_z) & !is.na(fit$sigmag_z)

    if (eb_shrink) {
        log_sig <- ifelse(is.na(fit$sigmag), NA_real_,
                          log(pmax(fit$sigmag, .Machine$double.eps)))
        ## SD shrinkage
        sh <- .eb_shrink(log_sig, ng_v,
                         mug = fit$mug_z, sigmag = fit$sigmag_z,
                         cuts_std = cutpoints)
        fit$sigmag    <- exp(sh$log_sig)
        fit$sigmag_z  <- fit$sigmag
        ## Mean shrinkage via per-district MAP under the actual bin
        ## likelihood with shrunk sigma_g and a Gaussian prior on mu.
        sh_mu <- .map_shrink_mu(ngk, fit$mug_z,
                                sigmag = fit$sigmag_z,
                                cuts_std = cutpoints)
        fit$mug_z <- sh_mu$mu
        fit$mug   <- fit$mug_z
        iter_info$eb_shrink     <- TRUE
        iter_info$eb_tau2       <- sh$tau2
        iter_info$eb_log_mean   <- sh$overall_log_mean
        iter_info$eb_weights    <- sh$w
        iter_info$eb_tau2_mu    <- sh_mu$tau2
        iter_info$eb_mean_mu    <- sh_mu$overall_mean

        w_mu <- rep(NA_real_, G)
        w_mu[ok_se] <- sh_mu$tau2 / (sh_mu$tau2 + sh_mu$s2[ok_se])
        mug_se_total        <- sqrt(w_mu * sh_mu$s2)
        log_sigmag_se_total <- sqrt(sh$w * sh$s2)
        sigmag_se_total     <- fit$sigmag * log_sigmag_se_total

        ## The sampling-vs-binning decomposition is a decomposition of a
        ## Fisher-information sampling SE, which is not what the EB
        ## posterior SE is.  For consistency with the MLE branch's
        ## output schema, expose _sampling and _binning fields but set
        ## them to NA in EB mode.
        mug_se_sampling        <- rep(NA_real_, G)
        mug_se_binning         <- rep(NA_real_, G)
        log_sigmag_se_sampling <- rep(NA_real_, G)
        log_sigmag_se_binning  <- rep(NA_real_, G)
        sigmag_se_sampling     <- rep(NA_real_, G)
        sigmag_se_binning      <- rep(NA_real_, G)

        ## The Huber-White sandwich SE is a model-misspecification-robust
        ## SE for the per-district MLE.  It has no analog for a shrunk
        ## posterior SD, so it is NA in EB mode.
        mug_se_sandwich        <- rep(NA_real_, G)
        log_sigmag_se_sandwich <- rep(NA_real_, G)
        sigmag_se_sandwich     <- rep(NA_real_, G)
        ## Same for the binning-only robust variant: it is derived from
        ## the sandwich, so it is unavailable in EB mode as well.
        mug_se_sw_binning        <- rep(NA_real_, G)
        log_sigmag_se_sw_binning <- rep(NA_real_, G)
        sigmag_se_sw_binning     <- rep(NA_real_, G)
    } else {
        s2_mu     <- rep(NA_real_, G)
        s2_logsig <- rep(NA_real_, G)
        s2_mu[ok_se]     <- .fisher_s2_mu(fit$mug_z[ok_se], fit$sigmag_z[ok_se],
                                           ng_v[ok_se], cutpoints)
        s2_logsig[ok_se] <- .fisher_s2_log_sigma(fit$mug_z[ok_se], fit$sigmag_z[ok_se],
                                                  ng_v[ok_se], cutpoints)

        ## Salvaged-district SE adjustments.
        ##
        ## Two-bin and one-interior-bin salvages have real MLEs for mu
        ## under the assumption that sigma equals the pool-borrowed
        ## value.  Their appropriate mu SE therefore comes from the
        ## sigma-known Fisher information (1/I11), NOT from the
        ## sigma-unknown formula I22/det_I above (which is degenerate
        ## for these districts because sigma isn't identified from
        ## their data).
        salvaged_mle_mu <- c(iter_info$salvaged_two_bin,
                             iter_info$salvaged_one_interior)
        if (length(salvaged_mle_mu) > 0) {
            s2_mu[salvaged_mle_mu] <- .fisher_s2_mu_known_sigma(
                fit$mug_z[salvaged_mle_mu],
                fit$sigmag_z[salvaged_mle_mu],
                ng_v[salvaged_mle_mu], cutpoints)
        }
        ## One-extreme-bin salvages assign a heuristic max/min mu that
        ## is not an MLE, so no valid Fisher SE exists.
        salvaged_extreme <- c(iter_info$salvaged_one_top,
                              iter_info$salvaged_one_bottom)
        if (length(salvaged_extreme) > 0) {
            s2_mu[salvaged_extreme] <- NA_real_
        }
        ## All salvaged districts: sigma was imputed, not estimated,
        ## so no valid sigma SE.  This propagates NA to sigmag_se_total
        ## (and its sampling/binning components) and to the sigma CI.
        all_salvaged <- c(salvaged_mle_mu, salvaged_extreme)
        if (length(all_salvaged) > 0) {
            s2_logsig[all_salvaged] <- NA_real_
        }

        ## Sampling-vs-binning decomposition of the Fisher-information
        ## SE.  Var_total = Var_sampling + Var_binning by the
        ## Pythagorean identity for information: Var_sampling is what
        ## you would report from n_g continuous unbinned normal
        ## observations, and Var_binning is the additional variance
        ## from having only bin counts.
        ##
        ## For mu:    Var_sampling = sigma_g^2 / n_g   (MLE of the mean of a normal)
        ## For log-sigma: Var_sampling = 1 / (2 n_g)   (delta method on chi-square(s^2))
        ##
        ## The max(., 0) guards are for numerical safety.  Var_total
        ## must be >= Var_sampling in exact arithmetic because binning
        ## cannot add information, but small numerical residuals can
        ## flip the sign.
        s2_mu_sampling     <- fit$sigmag_z^2 / ng_v
        s2_mu_binning      <- pmax(s2_mu - s2_mu_sampling, 0, na.rm = FALSE)
        s2_logsig_sampling <- 1 / (2 * ng_v)
        s2_logsig_binning  <- pmax(s2_logsig - s2_logsig_sampling, 0, na.rm = FALSE)

        ## Cascade the total-NA pattern (unidentified districts, and
        ## salvaged-one-extreme districts for mu; all salvaged districts
        ## for sigma) into the sampling/binning components too.
        mu_na  <- is.na(s2_mu)
        sig_na <- is.na(s2_logsig)
        s2_mu_sampling[mu_na]      <- NA_real_
        s2_mu_binning[mu_na]       <- NA_real_
        s2_logsig_sampling[sig_na] <- NA_real_
        s2_logsig_binning[sig_na]  <- NA_real_

        mug_se_total        <- sqrt(s2_mu)
        mug_se_sampling     <- sqrt(s2_mu_sampling)
        mug_se_binning      <- sqrt(s2_mu_binning)

        log_sigmag_se_total    <- sqrt(s2_logsig)
        log_sigmag_se_sampling <- sqrt(s2_logsig_sampling)
        log_sigmag_se_binning  <- sqrt(s2_logsig_binning)

        ## delta method: Var(sigma) = sigma^2 * Var(log sigma)
        sigmag_se_total     <- fit$sigmag * log_sigmag_se_total
        sigmag_se_sampling  <- fit$sigmag * log_sigmag_se_sampling
        sigmag_se_binning   <- fit$sigmag * log_sigmag_se_binning

        ## Gated off: `robust` is pinned FALSE, so none of this is used.
        if (robust) {
            ## Huber-White sandwich SE.  A fourth SE flavor, computed on the
            ## same MLE point estimate as the _total Fisher SE but using
            ## observed bin counts (not model probabilities) in the "meat".
            ## Robust to misspecification of the within-group normal shape:
            ## under a correctly specified model, E[n_{g,k}] = n_g * p_k, so
            ## the sandwich variance collapses to the Fisher-info variance
            ## (V_sandwich = I_matrix^{-1} * (n_g * I) * I_matrix^{-1} =
            ## I_matrix^{-1}).  Under misspecification, J diverges from
            ## n_g * I and the sandwich picks up the discrepancy.  This is
            ## the appropriate SE when the goodness-of-fit test rejects
            ## within-district normality.
            s2_mu_sandwich     <- rep(NA_real_, G)
            s2_logsig_sandwich <- rep(NA_real_, G)
            if (any(ok_se)) {
                sw <- .fisher_s2_sandwich(fit$mug_z[ok_se], fit$sigmag_z[ok_se],
                                          ng_v[ok_se], cutpoints,
                                          ngk[ok_se, , drop = FALSE])
                s2_mu_sandwich[ok_se]     <- sw$s2_mu_sandwich
                s2_logsig_sandwich[ok_se] <- sw$s2_logsig_sandwich
            }
            ## Salvaged 2-bin / 1-interior districts: use sigma-known
            ## sandwich (J11 / (n_g * I11)^2); log-sigma sandwich is NA
            ## because sigma is imputed, not estimated.
            if (length(salvaged_mle_mu) > 0) {
                s2_mu_sandwich[salvaged_mle_mu] <-
                    .fisher_s2_mu_known_sigma_sandwich(
                        fit$mug_z[salvaged_mle_mu],
                        fit$sigmag_z[salvaged_mle_mu],
                        ng_v[salvaged_mle_mu], cutpoints,
                        ngk[salvaged_mle_mu, , drop = FALSE])
                s2_logsig_sandwich[salvaged_mle_mu] <- NA_real_
            }
            ## Salvaged one-extreme-bin districts: no valid sandwich either
            ## for mu (heuristic max/min, not an MLE) or for sigma.
            if (length(salvaged_extreme) > 0) {
                s2_mu_sandwich[salvaged_extreme]     <- NA_real_
                s2_logsig_sandwich[salvaged_extreme] <- NA_real_
            }
            mug_se_sandwich        <- sqrt(s2_mu_sandwich)
            log_sigmag_se_sandwich <- sqrt(s2_logsig_sandwich)
            sigmag_se_sandwich     <- fit$sigmag * log_sigmag_se_sandwich

            ## Binning-only counterpart of the sandwich, by the same
            ## subtraction used for the Fisher version above: remove the
            ## variance the estimator would have with continuous scores,
            ## leaving the part attributable to seeing bin counts only.
            ## Var_sampling = sigma^2/n_g for mu and 1/(2 n_g) for log sigma
            ## are distribution-free given the true within-group variance,
            ## so subtracting them from a robust total is coherent.
            ##
            ## The subtraction can go negative here, unlike the Fisher
            ## version: coarsening cannot add information, so the Fisher
            ## total always exceeds the sampling variance, but the sandwich
            ## only equals the Fisher variance IN EXPECTATION (under correct
            ## specification E[J] = n_g * I), and with only K cell counts
            ## behind J it fluctuates around it.  A district goes negative
            ## when its observed counts imply less variability than the
            ## sampling variance alone -- which requires giving up the whole
            ## binning share, so it is rare, and likeliest in well-fitting
            ## districts where the sandwich sits near the Fisher variance
            ## rather than being inflated above it.
            ##
            ## Where that happens we fall back to the Fisher binning-only
            ## variance for that district.  Returning 0 would assert
            ## certainty and NA would discard a district needlessly; the
            ## Fisher value is a legitimate non-robust variance and is
            ## guaranteed non-negative.  The count of substitutions is
            ## exposed in iter_info so the fallback is never silent.
            neg_mu  <- is.finite(s2_mu_sandwich) & is.finite(s2_mu_sampling) &
                       (s2_mu_sandwich <= s2_mu_sampling)
            neg_sig <- is.finite(s2_logsig_sandwich) &
                       is.finite(s2_logsig_sampling) &
                       (s2_logsig_sandwich <= s2_logsig_sampling)

            s2_mu_sw_binning            <- s2_mu_sandwich - s2_mu_sampling
            s2_mu_sw_binning[neg_mu]    <- s2_mu_binning[neg_mu]
            s2_logsig_sw_binning        <- s2_logsig_sandwich - s2_logsig_sampling
            s2_logsig_sw_binning[neg_sig] <- s2_logsig_binning[neg_sig]

            s2_mu_sw_binning[is.na(s2_mu_sandwich)]         <- NA_real_
            s2_logsig_sw_binning[is.na(s2_logsig_sandwich)] <- NA_real_

            mug_se_sw_binning        <- sqrt(s2_mu_sw_binning)
            log_sigmag_se_sw_binning <- sqrt(s2_logsig_sw_binning)
            sigmag_se_sw_binning     <- fit$sigmag * log_sigmag_se_sw_binning
        }
    }

    pg    <- fit$pg
    mug   <- fit$mug
    sigmag <- fit$sigmag

    ## Per-district Pearson chi-square goodness-of-fit test of the
    ## within-district normality assumption.  For each district g, the
    ## statistic compares observed bin counts n_{gk} to the expected
    ## counts under the fitted normal N(mug, sigmag).  Under the null,
    ##     df = K - 1 - 2 = K - 3
    ## where K is the number of BINS, not the number of bins that
    ## happen to be occupied.
    ##
    ## An earlier version summed only over bins with a positive
    ## observed count and set df from that count, so a four-bin
    ## district with one empty bin got df = 0 and no test.  That was
    ## wrong on two counts.  An observed zero is data: the fitted
    ## normal has infinite support, so it assigns that bin a positive
    ## probability, and seeing nobody there is evidence against the
    ## fit.  Dropping the cell also breaks the correspondence with the
    ## model that was fitted, since the remaining probabilities are
    ## then no longer normalized to what the likelihood used.  What
    ## small cells actually threaten is the chi-square APPROXIMATION,
    ## and that depends on the EXPECTED count, not the observed one --
    ## so `min_exp` is reported alongside the test and the caller
    ## applies whatever screen they prefer (Cochran's rule of thumb is
    ## min_exp >= 5).
    ##
    ## Power is the more important limitation in practice.  With df = 1
    ## the detectable misfit scales as 1/n_g, so in a district of a few
    ## dozen students only a grossly non-normal distribution will
    ## register.  A district that does not reject is not thereby shown
    ## to be normal.
    G_ <- nrow(ngk)
    K_      <- ncol(ngk)
    chisq_g <- rep(NA_real_, G_)
    df_g    <- rep(NA_real_, G_)
    pval_g  <- rep(NA_real_, G_)
    minexp_g <- rep(NA_real_, G_)
    ng_v <- rowSums(ngk)
    d    <- K_ - 3
    for (g in seq_len(G_)) {
        if (is.na(mug[g]) || is.na(sigmag[g]) || sigmag[g] <= 0) next
        if (d < 1) next
        z       <- c(-Inf, (cutpoints - mug[g]) / sigmag[g], Inf)
        p_k     <- diff(stats::pnorm(z))
        exp_k   <- ng_v[g] * p_k
        obs_k   <- ngk[g, ]
        ## Guard only against a strictly zero expected count, which
        ## would divide by zero; that happens only in degenerate fits.
        keep    <- exp_k > 0
        if (!all(keep)) next
        chisq_g[g]  <- sum((obs_k - exp_k)^2 / exp_k)
        df_g[g]     <- d
        pval_g[g]   <- stats::pchisq(chisq_g[g], df = d, lower.tail = FALSE)
        minexp_g[g] <- min(exp_k)
    }
    gof <- data.frame(chisq = chisq_g, df = df_g, p = pval_g,
                      min_exp = minexp_g)

    ## Summary statistics across districts use a renormalized weight
    ## vector that excludes unidentified (NA) districts.
    ok <- !is.na(mug) & !is.na(sigmag)
    pg_ok <- pg
    pg_ok[!ok] <- 0
    if (sum(pg_ok) > 0) pg_ok <- pg_ok / sum(pg_ok)

    ## ICC on the raw scale: between-group variance / total variance
    raw_mean    <- sum(pg_ok * mug, na.rm = TRUE)
    between_var <- sum(pg_ok * (mug - raw_mean)^2, na.rm = TRUE)
    within_var  <- sum(pg_ok * sigmag^2, na.rm = TRUE)
    icc_raw     <- between_var / (between_var + within_var)

    ## Reporting scale for est_std: the pooled (population-weighted,
    ## within + between) distribution is placed at mean `pooled_mean`
    ## and SD `pooled_sd`.  The defaults (0, 1) give the usual
    ## standardized scale.
    ##
    ## Supplying other values is useful when the cutpoints are unknown
    ## but the pooled moments are not -- a state may publish its overall
    ## mean and SD without publishing the cut scores.  Passing them here
    ## puts the group estimates and the estimated cutpoints directly on
    ## the native score scale, rather than leaving the caller to rescale
    ## afterwards.  When cutpoints_known = TRUE the scale is already pinned by
    ## the cut scores, so these arguments are rejected above and the
    ## defaults leave this a pure standardization.
    raw_total_sd <- sqrt(between_var + within_var)
    scale_mult   <- pooled_sd / raw_total_sd
    mug_std      <- (mug - raw_mean) * scale_mult + pooled_mean
    sigmag_std   <- sigmag * scale_mult
    cuts_std     <- (cutpoints - raw_mean) * scale_mult + pooled_mean

    ## Confidence intervals, at the requested conf.level. The mean's CI
    ## is the ordinary symmetric Wald interval (mu's sampling
    ## distribution is already close to normal). The SD's CI is built
    ## on the log scale -- where the sampling distribution is much
    ## closer to normal than sigma itself -- and exponentiated back, so
    ## it is asymmetric around sigma_g rather than assuming symmetry
    ## where none exists.
    ##
    ## Two flavors of CI are reported per parameter:
    ##   `_total`   built on the *total* SE (sampling+binning combined).
    ##              Appropriate when the students are treated as a
    ##              sample from a within-district superpopulation, so
    ##              finite-sample sampling variability is a genuine
    ##              source of uncertainty about the parameter.
    ##   `_binning` built on the *binning-only* SE, and thus reflects
    ##              only the uncertainty attributable to coarsening
    ##              scores into K bins. Appropriate for population /
    ##              administrative data (e.g. state education-agency
    ##              reports) where the district's students are the full
    ##              population and there is no sampling variability to
    ##              account for -- only the binning is a real unknown.
    ##
    ## log_sigmag_se_total and log_sigmag_se_binning are both scale-free
    ## (see note above), so only the log-scale center shifts between
    ## raw and standardized.
    mug_ci_lower_raw_total    <- mug - z_crit * mug_se_total
    mug_ci_upper_raw_total    <- mug + z_crit * mug_se_total
    sigmag_ci_lower_raw_total <- exp(log(sigmag) - z_crit * log_sigmag_se_total)
    sigmag_ci_upper_raw_total <- exp(log(sigmag) + z_crit * log_sigmag_se_total)

    mug_ci_lower_raw_binning    <- mug - z_crit * mug_se_binning
    mug_ci_upper_raw_binning    <- mug + z_crit * mug_se_binning
    sigmag_ci_lower_raw_binning <- exp(log(sigmag) - z_crit * log_sigmag_se_binning)
    sigmag_ci_upper_raw_binning <- exp(log(sigmag) + z_crit * log_sigmag_se_binning)

    ## Test-score-scale sandwich intervals; see the `robust` note above.
    if (robust) {
        mug_ci_lower_raw_sandwich    <- mug - z_crit * mug_se_sandwich
        mug_ci_upper_raw_sandwich    <- mug + z_crit * mug_se_sandwich
        sigmag_ci_lower_raw_sandwich <- exp(log(sigmag) - z_crit * log_sigmag_se_sandwich)
        sigmag_ci_upper_raw_sandwich <- exp(log(sigmag) + z_crit * log_sigmag_se_sandwich)
    }

    mug_se_total_std           <- mug_se_total   * scale_mult
    mug_se_binning_std         <- mug_se_binning * scale_mult
    mug_ci_lower_std_total     <- mug_std - z_crit * mug_se_total_std
    mug_ci_upper_std_total     <- mug_std + z_crit * mug_se_total_std
    sigmag_ci_lower_std_total  <- exp(log(sigmag_std) - z_crit * log_sigmag_se_total)
    sigmag_ci_upper_std_total  <- exp(log(sigmag_std) + z_crit * log_sigmag_se_total)

    mug_ci_lower_std_binning     <- mug_std - z_crit * mug_se_binning_std
    mug_ci_upper_std_binning     <- mug_std + z_crit * mug_se_binning_std
    sigmag_ci_lower_std_binning  <- exp(log(sigmag_std) - z_crit * log_sigmag_se_binning)
    sigmag_ci_upper_std_binning  <- exp(log(sigmag_std) + z_crit * log_sigmag_se_binning)

    ## Standardized-scale sandwich intervals; see the `robust` note above.
    if (robust) {
        mug_se_sandwich_std          <- mug_se_sandwich * scale_mult
        mug_ci_lower_std_sandwich    <- mug_std - z_crit * mug_se_sandwich_std
        mug_ci_upper_std_sandwich    <- mug_std + z_crit * mug_se_sandwich_std
        sigmag_ci_lower_std_sandwich <- exp(log(sigmag_std) - z_crit * log_sigmag_se_sandwich)
        sigmag_ci_upper_std_sandwich <- exp(log(sigmag_std) + z_crit * log_sigmag_se_sandwich)

        ## Robust (sandwich) binning-only intervals, for scope =
        ## "population" with robust = TRUE.
        mug_ci_lower_raw_sw_binning    <- mug - z_crit * mug_se_sw_binning
        mug_ci_upper_raw_sw_binning    <- mug + z_crit * mug_se_sw_binning
        sigmag_ci_lower_raw_sw_binning <- exp(log(sigmag) - z_crit * log_sigmag_se_sw_binning)
        sigmag_ci_upper_raw_sw_binning <- exp(log(sigmag) + z_crit * log_sigmag_se_sw_binning)

        mug_se_sw_binning_std          <- mug_se_sw_binning * scale_mult
        mug_ci_lower_std_sw_binning    <- mug_std - z_crit * mug_se_sw_binning_std
        mug_ci_upper_std_sw_binning    <- mug_std + z_crit * mug_se_sw_binning_std
        sigmag_ci_lower_std_sw_binning <- exp(log(sigmag_std) - z_crit * log_sigmag_se_sw_binning)
        sigmag_ci_upper_std_sw_binning <- exp(log(sigmag_std) + z_crit * log_sigmag_se_sw_binning)
    }

    ## The group-level estimates are always called mean and
    ## sd, whichever estimator produced them.  Earlier versions
    ## suffixed them _mle or _eb, which made every downstream script
    ## hard-code the estimator into its field names.  The estimator is
    ## echoed back as a field instead, as `scope` is, so the object is
    ## still self-describing without the name churn.
    mean_nm <- "mean"
    sd_nm   <- "sd"

    ## Standard errors accompanying whichever point estimate was actually
    ## returned (unshrunk MLE Fisher SE, or approximate EB posterior SE;
    ## see the computation above). NA for unidentified districts.
    ## Confidence-interval field names carry the same conf.level for
    ## both the mean and the SD, and come in `_total` and `_binning`
    ## flavors (see the CI block above for what each flavor means).
    ## Select the one SE/CI flavor that matches `scope`.  All flavors
    ## were computed above; only the appropriate one is returned, so the
    ## output cannot be read with the wrong notion of uncertainty.
    ##
    ##   scope = "population": the group's units ARE the population, so
    ##     its true mean is the actual mean of those units' scores.
    ##     There is no sampling error -- the only uncertainty is that we
    ##     observe bin counts rather than individual scores.  Report the
    ##     binning-only SE.
    ##   scope = "sample": the units are a sample from a larger
    ##     population, so the target is that population's parameter and
    ##     the error has both components.  Report the total SE, which
    ##     satisfies Var_total = Var_sampling + Var_binning.
    if (scope == "population" && !robust) {
        se_mu_raw      <- mug_se_binning
        se_sd_raw      <- sigmag_se_binning
        ci_mu_lo_raw   <- mug_ci_lower_raw_binning
        ci_mu_hi_raw   <- mug_ci_upper_raw_binning
        ci_sd_lo_raw   <- sigmag_ci_lower_raw_binning
        ci_sd_hi_raw   <- sigmag_ci_upper_raw_binning
        ci_mu_lo_std   <- mug_ci_lower_std_binning
        ci_mu_hi_std   <- mug_ci_upper_std_binning
        ci_sd_lo_std   <- sigmag_ci_lower_std_binning
        ci_sd_hi_std   <- sigmag_ci_upper_std_binning
    } else if (scope == "population" && robust) {
        se_mu_raw      <- mug_se_sw_binning
        se_sd_raw      <- sigmag_se_sw_binning
        ci_mu_lo_raw   <- mug_ci_lower_raw_sw_binning
        ci_mu_hi_raw   <- mug_ci_upper_raw_sw_binning
        ci_sd_lo_raw   <- sigmag_ci_lower_raw_sw_binning
        ci_sd_hi_raw   <- sigmag_ci_upper_raw_sw_binning
        ci_mu_lo_std   <- mug_ci_lower_std_sw_binning
        ci_mu_hi_std   <- mug_ci_upper_std_sw_binning
        ci_sd_lo_std   <- sigmag_ci_lower_std_sw_binning
        ci_sd_hi_std   <- sigmag_ci_upper_std_sw_binning
    } else if (scope == "sample" && !robust) {
        se_mu_raw      <- mug_se_total
        se_sd_raw      <- sigmag_se_total
        ci_mu_lo_raw   <- mug_ci_lower_raw_total
        ci_mu_hi_raw   <- mug_ci_upper_raw_total
        ci_sd_lo_raw   <- sigmag_ci_lower_raw_total
        ci_sd_hi_raw   <- sigmag_ci_upper_raw_total
        ci_mu_lo_std   <- mug_ci_lower_std_total
        ci_mu_hi_std   <- mug_ci_upper_std_total
        ci_sd_lo_std   <- sigmag_ci_lower_std_total
        ci_sd_hi_std   <- sigmag_ci_upper_std_total
    } else {
        se_mu_raw      <- mug_se_sandwich
        se_sd_raw      <- sigmag_se_sandwich
        ci_mu_lo_raw   <- mug_ci_lower_raw_sandwich
        ci_mu_hi_raw   <- mug_ci_upper_raw_sandwich
        ci_sd_lo_raw   <- sigmag_ci_lower_raw_sandwich
        ci_sd_hi_raw   <- sigmag_ci_upper_raw_sandwich
        ci_mu_lo_std   <- mug_ci_lower_std_sandwich
        ci_mu_hi_std   <- mug_ci_upper_std_sandwich
        ci_sd_lo_std   <- sigmag_ci_lower_std_sandwich
        ci_sd_hi_std   <- sigmag_ci_upper_std_sandwich
    }

    est_std <- list()
    est_std[[mean_nm]]      <- mug_std
    est_std[[sd_nm]]        <- sigmag_std
    est_std$mean_se         <- se_mu_raw * scale_mult
    est_std$sd_se           <- se_sd_raw * scale_mult
    est_std$mean_ci_lower   <- ci_mu_lo_std
    est_std$mean_ci_upper   <- ci_mu_hi_std
    est_std$sd_ci_lower     <- ci_sd_lo_std
    est_std$sd_ci_upper     <- ci_sd_hi_std
    est_std$scope           <- scope
    est_std$estimator       <- estimator
    est_std$cutpoints_known <- cutpoints_known
    ## The scale est_std is reported on: the pooled distribution is
    ## placed at this mean and SD.  (0, 1) is the usual standardized
    ## scale; other values put the estimates on whatever external scale
    ## the caller anchored them to.
    est_std$pooled_mean     <- pooled_mean
    est_std$pooled_sd       <- pooled_sd
    est_std$conf.level      <- conf.level
    est_std$cutpoints                 <- cuts_std
    est_std$icc                       <- icc_raw   ## ICC is scale-invariant

    ## Record the tolerances actually in force.  Both are parameter-change
    ## criteria in SD units, but they govern different loops and are easy
    ## to confuse, so a fitted object should be able to report what it
    ## used rather than leaving the caller to recall which defaults
    ## applied.  Set here, after every branch above has finished writing
    ## to iter_info (the EM branch replaces the list wholesale).
    iter_info$tol               <- tol
    iter_info$tol_cutpoints_std <- tol_cutpoints_std

    ## Tell the caller which kind of uncertainty they are holding.  The
    ## two are not interchangeable and the field names no longer say
    ## which is which, so state it once per call rather than leaving it
    ## to be inferred from the argument list.
    message(sprintf(
        paste0("fast_hetop: scope = \"%s\". The reported SEs ",
               "and %g%% CIs reflect\n  %s,\n  estimated from %s."),
        scope, 100 * conf.level,
        if (scope == "population") {
            paste0("binning error only -- the coarsening from observing bin ",
                   "counts\n  rather than individual scores; each group's ",
                   "units are treated as its\n  whole population, so there ",
                   "is no sampling error")
        } else {
            paste0("sampling and binning error combined -- both drawing ",
                   "these units\n  rather than others, and observing only ",
                   "bin counts rather than\n  individual scores")
        },
        paste0("the model-implied Fisher information, which assumes\n",
               "  within-group normality holds")))
    if (!cutpoints_known) {
        message(sprintf(
            paste0("  Cutpoints were estimated from the pooled bin ",
                   "proportions; est_std is\n  reported on the scale where ",
                   "the pooled distribution has mean %g and SD %g%s."),
            pooled_mean, pooled_sd,
            if (identical(pooled_mean, 0) && identical(pooled_sd, 1)) {
                " (the default\n  standardized scale)"
            } else ""))
    }

    ## `est_raw` is meaningful only when cutpoints were supplied by the
    ## caller (then it lives on the test-score scale).  When cutpoints
    ## are estimated from the data, the "raw" scale is whatever
    ## arbitrary location and scale the EM converges to, which has no
    ## external interpretation, so only `est_std` is returned.
    if (scale_known) {
        est_raw <- list()
        est_raw[[mean_nm]]      <- mug
        est_raw[[sd_nm]]        <- sigmag
        est_raw$mean_se         <- se_mu_raw
        est_raw$sd_se           <- se_sd_raw
        est_raw$mean_ci_lower   <- ci_mu_lo_raw
        est_raw$mean_ci_upper   <- ci_mu_hi_raw
        est_raw$sd_ci_lower     <- ci_sd_lo_raw
        est_raw$sd_ci_upper     <- ci_sd_hi_raw
        est_raw$scope           <- scope
        est_raw$estimator       <- estimator
        est_raw$cutpoints_known <- cutpoints_known
        est_raw$conf.level      <- conf.level
        est_raw$cutpoints       <- cutpoints
        est_raw$icc             <- icc_raw

        list(est_raw = est_raw, est_std = est_std,
             gof = gof, iter_info = iter_info)
    } else {
        list(est_std = est_std,
             gof = gof, iter_info = iter_info)
    }
}
