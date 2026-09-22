## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>"
)

## -----------------------------------------------------------------------------
library(binest)
data(tx_g6_math_2018)
dim(tx_g6_math_2018)
head(tx_g6_math_2018, 3)

## -----------------------------------------------------------------------------
ngk  <- with(tx_g6_math_2018,
             cbind(unsatisfactory, approaches, meets, masters))
cuts <- c(1536, 1653, 1772)
truth <- tx_g6_math_2018$reported_mean

## -----------------------------------------------------------------------------
t0 <- Sys.time()
fit_bm_known <- fast_hetop(ngk, cutpoints_known = TRUE, cutpoints = cuts,
                           scope = "population")
t_bm_known   <- as.numeric(Sys.time() - t0, units = "secs")
cor(fit_bm_known$est_raw$mean, truth)

## -----------------------------------------------------------------------------
fit_bm_sample <- fast_hetop(ngk, cutpoints_known = TRUE, cutpoints = cuts,
                            scope = "sample")
head(cbind(population = fit_bm_known$est_raw$mean_se,
           sample     = fit_bm_sample$est_raw$mean_se))

## ----fig.width = 4.5, fig.height = 4.5----------------------------------------
plot(truth, fit_bm_known$est_raw$mean,
     pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.3),
     xlab = "True district mean",
     ylab = "Estimated mean (test-score scale)",
     main = "fast_hetop (known cuts)")
abline(0, 1, col = "red", lty = 2)

## -----------------------------------------------------------------------------
t0 <- Sys.time()
fit_bm_null <- fast_hetop(ngk, scope = "population")
t_bm_null   <- as.numeric(Sys.time() - t0, units = "secs")

## ----fig.width = 4.5, fig.height = 4.5----------------------------------------
plot(truth, fit_bm_null$est_std$mean,
     pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.3),
     xlab = "True district mean",
     ylab = "Estimated mean (standardized)",
     main = "fast_hetop (cuts from data)")

## -----------------------------------------------------------------------------
fit_bm_eb <- fast_hetop(ngk, cutpoints_known = TRUE, cutpoints = cuts,
                        scope = "sample", estimator = "EB_shrunk")

## ----fig.width = 4.5, fig.height = 4.5----------------------------------------
plot(truth, fit_bm_eb$est_raw$mean,
     pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.3),
     xlab = "True district mean",
     ylab = "EB-shrunk estimated mean",
     main = "fast_hetop (EB_shrunk, known cuts)")

## -----------------------------------------------------------------------------
set.seed(1)
sub <- sample(nrow(ngk), 50)
t0 <- Sys.time()
fit_mle <- mle_hetop(ngk[sub, ], iterlim = 200)
t_mle <- as.numeric(Sys.time() - t0, units = "secs")
cor(fit_mle$est_star$mug, truth[sub])

## ----fig.width = 4.5, fig.height = 4.5----------------------------------------
plot(truth[sub], fit_mle$est_star$mug,
     pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.3),
     xlab = "True district mean",
     ylab = "Estimated mean (standardized)",
     main = "HETOP MLE (50-district subsample)")

## ----eval = FALSE-------------------------------------------------------------
# t0 <- Sys.time()
# fit_fh <- fh_hetop(
#   ngk       = ngk,
#   p         = c(10, 10),
#   m         = c(100, 100),
#   gridL     = c(-5.0, log(0.10)),
#   gridU     = c( 5.0, log(5.0)),
#   n.iter    = 2000,
#   n.burnin  = 1000,
#   seed      = 3142
# )
# t_fh <- as.numeric(Sys.time() - t0, units = "secs")
# cor(fit_fh$fh_hetop_extras$est_star_mug$theta_pm, truth)

## ----fig.width = 4.5, fig.height = 4.5----------------------------------------
fh_means_file <- system.file("extdata",
                             "fh_hetop_means_tx_g6_math_2018.rds",
                             package = "binest")
fh_means <- readRDS(fh_means_file)

plot(truth, fh_means,
     pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.3),
     xlab = "True district mean",
     ylab = "Estimated mean (standardized)",
     main = "HETOP Bayes (full data, cached)")

