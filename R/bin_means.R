###############################################################################
## bin_means: fast estimator of group means and SDs from binned count
## data, under the assumption that within each group the underlying
## scores are normally distributed.
##
## Arguments:
##   ngk        GxK numeric matrix of bin counts.
##   cutpoints  Optional numeric vector of length K-1.  If supplied,
##              treated as the K-1 known cutpoints on the test-score
##              scale.  If NULL (default), cutpoints are derived from the
##              pooled bin proportions via qnorm(cumsum(props)).
##   iterate    Logical.  If TRUE and cutpoints = NULL, refine the
##              cutpoints by alternating between (a) fitting per-district
##              means and SDs given current cutpoints, and (b) computing
##              new cutpoints by matching the pooled bin proportions to
##              the implied mixture-of-normals CDF.  Iterates until the
##              cutpoints change by less than `tol`.  Default FALSE.
##              Not valid when cutpoints are supplied.
##   tol        Convergence tolerance for iteration (default 1e-6).
##   maxit      Iteration cap (default 100).
##
## Output: a list with components est_raw, est_std, gof, iter_info.
##   est_raw   estimates on the raw (test-score) scale, when cutpoints
##             were supplied; otherwise on the standardized scale.
##   est_std   estimates on the standardized scale where the population-
##             weighted state mean is 0 and the total (within + between)
##             state SD is 1.
##   gof       per-district Pearson chi-square goodness-of-fit of the
##             within-district normality assumption (chisq, df, p).
##   iter_info diagnostic flags including `within`, `eb_shrink`, and
##             (when applicable) shrinkage tuning parameters.
###############################################################################

.bin_means_oneshot <- function(ngk, cutpoints, scale_known) {
    G  <- nrow(ngk)
    K  <- ncol(ngk)
    ng <- rowSums(ngk)
    pg <- ng / sum(ng)
    qgk <- ngk / ng

    pooled_props <- colSums(ngk) / sum(ngk)
    pooled_cum   <- cumsum(pooled_props)

    if (scale_known) {
        z_cuts_implied <- qnorm(pooled_cum[1:(K-1)])
        lmfit          <- stats::lm(cutpoints ~ z_cuts_implied)
        pooled_mean    <- as.numeric(stats::coef(lmfit)[1])
        pooled_sd      <- as.numeric(stats::coef(lmfit)[2])
        cuts_for_moments <- (cutpoints - pooled_mean) / pooled_sd
    } else {
        pooled_mean      <- 0
        pooled_sd        <- 1
        cuts_for_moments <- cutpoints
    }

    c_left  <- c(-Inf, cuts_for_moments)
    c_right <- c(cuts_for_moments, Inf)
    phi_left  <- dnorm(c_left)
    phi_right <- dnorm(c_right)
    Phi_left  <- pnorm(c_left)
    Phi_right <- pnorm(c_right)
    p_k       <- Phi_right - Phi_left
    bin_mean_z <- (phi_left - phi_right) / p_k
    cphi_left  <- ifelse(is.infinite(c_left),  0, c_left  * phi_left)
    cphi_right <- ifelse(is.infinite(c_right), 0, c_right * phi_right)
    bin_2nd_z  <- 1 + (cphi_left - cphi_right) / p_k

    mug_z    <- as.numeric(qgk %*% bin_mean_z)
    mug2_z_g <- as.numeric(qgk %*% bin_2nd_z)
    sigmag_z <- sqrt(pmax(mug2_z_g - mug_z^2, 0))

    if (scale_known) {
        mug    <- pooled_mean + pooled_sd * mug_z
        sigmag <- pooled_sd * sigmag_z
    } else {
        mug    <- mug_z
        sigmag <- sigmag_z
    }

    list(mug = mug, sigmag = sigmag, mug_z = mug_z, sigmag_z = sigmag_z,
         pg = pg, pooled_cum = pooled_cum,
         pooled_mean = pooled_mean, pooled_sd = pooled_sd)
}

.refine_within <- function(qgk, cuts_for_moments, mug_z, sigmag_z,
                           tol = 1e-3, maxit = 100) {
    ## For each group g, iteratively refine (mug, sigmag) so that the
    ## within-bin truncated-normal moments are computed using the
    ## group's own implied normal distribution rather than the pooled
    ## (mean 0, SD 1) normal.
    ##
    ## At each iteration, for each group g and each bin k with
    ## boundaries (c_{k-1}, c_k) on the standardized scale:
    ##   z_left  = (c_{k-1} - mug[g]) / sigmag[g]
    ##   z_right = (c_k     - mug[g]) / sigmag[g]
    ##   E[Y | bin k, group g]   = mug[g] + sigmag[g] *
    ##       (phi(z_left) - phi(z_right)) / (Phi(z_right) - Phi(z_left))
    ##   E[Y^2 | bin k, group g] is computed similarly.
    ## Then update mug, sigmag from the weighted sums of these moments.
    G <- length(mug_z)
    K <- length(cuts_for_moments) + 1
    c_left  <- c(-Inf, cuts_for_moments)
    c_right <- c(cuts_for_moments, Inf)
    for (iter in 1:maxit) {
        new_mug <- numeric(G)
        new_sig <- numeric(G)
        for (g in 1:G) {
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
            ## E[Y^2 | bin k] = mu^2 + 2 mu sigma * (phi_l - phi_r)/p_k
            ##               + sigma^2 * (1 + (z_l*phi_l - z_r*phi_r)/p_k)
            zphi_l <- ifelse(is.infinite(z_left),  0, z_left  * phi_l)
            zphi_r <- ifelse(is.infinite(z_right), 0, z_right * phi_r)
            b2 <- mug_z[g]^2 +
                  2 * mug_z[g] * sg * (phi_l - phi_r) / p_k +
                  sg^2 * (1 + (zphi_l - zphi_r) / p_k)
            new_mug[g] <- sum(qgk[g, ] * bm)
            new_v      <- sum(qgk[g, ] * b2) - new_mug[g]^2
            new_sig[g] <- sqrt(max(new_v, 0))
        }
        delta <- max(abs(new_mug - mug_z), abs(new_sig - sigmag_z),
                     na.rm = TRUE)
        mug_z   <- new_mug
        sigmag_z <- new_sig
        if (delta < tol) {
            attr(mug_z, "iterations") <- iter
            attr(mug_z, "converged")  <- TRUE
            return(list(mug = mug_z, sigmag = sigmag_z,
                        iterations = iter, converged = TRUE))
        }
    }
    list(mug = mug_z, sigmag = sigmag_z,
         iterations = maxit, converged = FALSE)
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

.eb_shrink <- function(log_sig_hat, n_g, mug = NULL, sigmag = NULL,
                       cuts_std = NULL, sampling_var = "continuous") {
    ## Fay-Herriot empirical-Bayes shrinkage of per-district log SDs.
    ## sampling_var = "continuous" uses s_g^2 = 1 / (2 n_g) (ignores binning).
    ## sampling_var = "fisher"     uses the Fisher information for the
    ##                             binned-normal model with known cutpoints.
    ## NA entries in log_sig_hat are treated as unidentified and passed
    ## through to the output; they are excluded from tau^2 / overall-mean
    ## estimation.
    G <- length(log_sig_hat)
    ok <- !is.na(log_sig_hat)
    if (sampling_var == "fisher") {
        if (is.null(mug) || is.null(sigmag) || is.null(cuts_std)) {
            stop("sampling_var = 'fisher' requires mug, sigmag, and cuts_std")
        }
        s2 <- rep(NA_real_, G)
        s2[ok] <- .fisher_s2_log_sigma(mug[ok], sigmag[ok], n_g[ok], cuts_std)
    } else {
        s2 <- 1 / (2 * pmax(n_g, 1))
        s2[!ok] <- NA_real_
    }
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

.eb_shrink_mu <- function(mu_hat, n_g, mug = NULL, sigmag = NULL,
                          cuts_std = NULL, sampling_var = "continuous") {
    ## Fay-Herriot empirical-Bayes shrinkage of per-district means.
    ## sampling_var = "continuous" uses s_g^2 = sigma_g^2 / n_g (ignores binning).
    ## sampling_var = "fisher"     uses the Fisher information for the
    ##                             binned-normal model with known cutpoints.
    if (sampling_var == "fisher") {
        if (is.null(mug) || is.null(sigmag) || is.null(cuts_std)) {
            stop("sampling_var = 'fisher' requires mug, sigmag, and cuts_std")
        }
        s2 <- .fisher_s2_mu(mug, sigmag, n_g, cuts_std)
    } else {
        if (is.null(sigmag)) {
            stop("continuous-data sampling variance for mu requires sigmag")
        }
        s2 <- sigmag^2 / pmax(n_g, 1)
    }
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

bin_means <- function(ngk, cutpoints = NULL, eb_shrink = FALSE,
                     iterate = FALSE, tol = 1e-3, maxit = 100) {
    ## Internal defaults retained as locals for the existing downstream
    ## code paths: always do per-district MLE, always use the
    ## Fisher-information sampling variance for EB shrinkage.
    within       <- TRUE
    sampling_var <- "fisher"

    ## iterate = FALSE is the default because the non-iterative
    ## estimates already match HETOP MLE in both the Texas data and the
    ## simulation; iteration has little room to help and did not
    ## meaningfully change the per-district estimates in either case.
    ## We expose the option for users who want it but do not iterate
    ## by default.

    if (!is.numeric(ngk))  stop("ngk must be a GxK numeric matrix of category counts")
    if (!is.matrix(ngk))   stop("ngk must be a GxK numeric matrix of category counts")
    if (any(is.na(ngk)))   stop("ngk cannot contain missing values")
    if (any(ngk < 0))      stop("ngk must contain only non-negative values")
    if (any(rowSums(ngk) <= 0)) stop("ngk contains at least one row with insufficient data")

    G <- nrow(ngk)
    K <- ncol(ngk)
    if (K <= 2) stop("Function requires K >= 3 categories")

    if (iterate && !is.null(cutpoints)) {
        stop("iterate = TRUE is not valid when cutpoints are supplied; ",
             "iteration only refines the data-derived cutpoints")
    }

    pooled_props <- colSums(ngk) / sum(ngk)
    pooled_cum   <- cumsum(pooled_props)
    if (any(pooled_props == 0)) {
        stop("At least one pooled bin has zero count; bin_means requires all pooled bins to be non-empty")
    }

    if (is.null(cutpoints)) {
        scale_known <- FALSE
        cutpoints   <- qnorm(pooled_cum[1:(K-1)])
    } else {
        scale_known <- TRUE
        if (!is.numeric(cutpoints))         stop("cutpoints must be a numeric vector of length K-1")
        if (length(cutpoints) != (K-1))     stop("cutpoints must have length K-1 = ", K-1)
        if (any(is.na(cutpoints)))          stop("cutpoints cannot contain missing values")
        if (any(diff(cutpoints) <= 0))      stop("cutpoints must be strictly increasing")
    }

    ## One-shot fit.
    fit <- .bin_means_oneshot(ngk, cutpoints, scale_known)

    iter_info <- list(iterated = FALSE, within = FALSE, eb_shrink = FALSE)

    ## Optional EM-style refinement when cutpoints are unknown.
    ## Alternates between (a) refitting per-district (mu_g, sigma_g)
    ## given current cutpoints and (b) refitting cutpoints to match the
    ## pooled bin proportions under the implied mixture-of-normals CDF,
    ## using the current per-district fits.  Iterates until cutpoints
    ## stop changing.
    if (iterate) {
        em_converged  <- FALSE
        em_iterations <- 0
        for (em_iter in seq_len(maxit)) {
            old_cuts <- cutpoints
            new_cuts <- .refine_cutpoints(fit$mug_z, fit$sigmag_z, fit$pg,
                                          pooled_cum, init_cuts = cutpoints,
                                          tol = tol, maxit = maxit)
            attributes(new_cuts) <- NULL
            cutpoints <- new_cuts
            fit <- .bin_means_oneshot(ngk, cutpoints, scale_known = FALSE)
            em_iterations <- em_iter
            if (max(abs(cutpoints - old_cuts), na.rm = TRUE) < tol) {
                em_converged <- TRUE
                break
            }
        }
        iter_info <- list(
            iterated   = TRUE,
            iterations = em_iterations,
            converged  = em_converged
        )
    }

    ## Optional within-district refinement: refit each group's
    ## (mu_g, sigma_g) using its own truncated-normal bin-conditional
    ## moments, iterating to a fixed point.  Districts with fewer than
    ## three populated bins do not have likelihood-identified (mu_g,
    ## sigma_g) and are returned as NA.
    if (within) {
        ng_v <- rowSums(ngk)
        qgk  <- ngk / ng_v
        if (scale_known) {
            cuts_z <- (cutpoints - fit$pooled_mean) / fit$pooled_sd
        } else {
            cuts_z <- cutpoints
        }
        bins_used <- rowSums(ngk > 0)
        unidentified <- bins_used < 3
        ref <- .refine_within(qgk, cuts_z, fit$mug_z, fit$sigmag_z,
                              tol = tol, maxit = maxit)
        fit$mug_z    <- ref$mug
        fit$sigmag_z <- ref$sigmag
        fit$mug_z[unidentified]    <- NA_real_
        fit$sigmag_z[unidentified] <- NA_real_
        if (scale_known) {
            fit$mug    <- fit$pooled_mean + fit$pooled_sd * fit$mug_z
            fit$sigmag <- fit$pooled_sd * fit$sigmag_z
        } else {
            fit$mug    <- fit$mug_z
            fit$sigmag <- fit$sigmag_z
        }
        iter_info$within            <- TRUE
        iter_info$within_iterations <- ref$iterations
        iter_info$within_converged  <- ref$converged
        iter_info$within_unidentified <- which(unidentified)
    }

    ## Optional empirical-Bayes shrinkage of the per-district SDs and means.
    if (eb_shrink) {
        ng_v <- rowSums(ngk)
        log_sig <- ifelse(is.na(fit$sigmag), NA_real_,
                          log(pmax(fit$sigmag, .Machine$double.eps)))
        if (scale_known) {
            cuts_std <- (cutpoints - fit$pooled_mean) / fit$pooled_sd
        } else {
            cuts_std <- cutpoints
        }
        ## SD shrinkage
        sh <- .eb_shrink(log_sig, ng_v,
                         mug = fit$mug_z, sigmag = fit$sigmag_z,
                         cuts_std = cuts_std,
                         sampling_var = sampling_var)
        fit$sigmag    <- exp(sh$log_sig)
        if (scale_known) {
            fit$sigmag_z <- fit$sigmag / fit$pooled_sd
        } else {
            fit$sigmag_z <- fit$sigmag
        }
        ## Mean shrinkage via per-district MAP under the actual bin
        ## likelihood with shrunk sigma_g and a Gaussian prior on mu.
        ## This avoids the Fay-Herriot Gaussian approximation, which can
        ## over-shrink sparse-count districts whose unshrunk mu_hat is
        ## extreme but well-identified by the bin counts.
        sh_mu <- .map_shrink_mu(ngk, fit$mug_z,
                                sigmag = fit$sigmag_z,
                                cuts_std = cuts_std)
        fit$mug_z <- sh_mu$mu
        if (scale_known) {
            fit$mug <- fit$pooled_mean + fit$pooled_sd * fit$mug_z
        } else {
            fit$mug <- fit$mug_z
        }
        iter_info$eb_shrink     <- TRUE
        iter_info$eb_sampling_var <- sampling_var
        iter_info$eb_tau2       <- sh$tau2
        iter_info$eb_log_mean   <- sh$overall_log_mean
        iter_info$eb_weights    <- sh$w
        iter_info$eb_tau2_mu    <- sh_mu$tau2
        iter_info$eb_mean_mu    <- sh_mu$overall_mean
    }

    pg    <- fit$pg
    mug   <- fit$mug
    sigmag <- fit$sigmag

    ## Per-district Pearson chi-square goodness-of-fit test of the
    ## within-district normality assumption.  For each district g, the
    ## statistic compares observed bin counts n_{gk} to the expected
    ## counts under the fitted normal N(mug, sigmag).  Under the null,
    ## the statistic has a chi-square distribution with
    ##     df = K_populated - 1 - 2 = K_populated - 3
    ## where K_populated is the number of bins with positive count.
    ## The test is unidentified (df = 0) when only three bins are
    ## populated, and returns NA in that case.
    G_ <- nrow(ngk)
    chisq_g <- rep(NA_real_, G_)
    df_g    <- rep(NA_real_, G_)
    pval_g  <- rep(NA_real_, G_)
    ng_v <- rowSums(ngk)
    for (g in seq_len(G_)) {
        if (is.na(mug[g]) || is.na(sigmag[g]) || sigmag[g] <= 0) next
        z       <- c(-Inf, (cutpoints - mug[g]) / sigmag[g], Inf)
        p_k     <- diff(stats::pnorm(z))
        exp_k   <- ng_v[g] * p_k
        obs_k   <- ngk[g, ]
        keep    <- obs_k > 0 & exp_k > 0
        kk      <- sum(keep)
        d       <- kk - 3
        if (d < 1) next
        chisq_g[g] <- sum((obs_k[keep] - exp_k[keep])^2 / exp_k[keep])
        df_g[g]    <- d
        pval_g[g]  <- stats::pchisq(chisq_g[g], df = d, lower.tail = FALSE)
    }
    gof <- data.frame(chisq = chisq_g, df = df_g, p = pval_g)

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

    ## Standardized scale: population-weighted state mean = 0, total
    ## (within + between) state SD = 1.
    raw_total_sd <- sqrt(between_var + within_var)
    mug_std      <- (mug - raw_mean) / raw_total_sd
    sigmag_std   <- sigmag / raw_total_sd
    cuts_std     <- (cutpoints - raw_mean) / raw_total_sd

    ## Suffix the group-level estimate names with _mle (per-district
    ## MLE) or _eb (empirical-Bayes shrunk) to make the output
    ## self-describing.
    suffix <- if (isTRUE(iter_info$eb_shrink)) "_eb" else "_mle"
    mean_nm <- paste0("group_mean", suffix)
    sd_nm   <- paste0("group_sd",   suffix)

    est_raw <- list()
    est_raw[[mean_nm]] <- mug
    est_raw[[sd_nm]]   <- sigmag
    est_raw$cutpoints  <- cutpoints
    est_raw$icc        <- icc_raw

    est_std <- list()
    est_std[[mean_nm]] <- mug_std
    est_std[[sd_nm]]   <- sigmag_std
    est_std$cutpoints  <- cuts_std
    est_std$icc        <- icc_raw   ## ICC is scale-invariant

    list(est_raw = est_raw, est_std = est_std,
         gof = gof, iter_info = iter_info)
}
