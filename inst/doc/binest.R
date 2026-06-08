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
fit_bm_known <- bin_means(ngk, cutpoints = cuts)
t_bm_known   <- as.numeric(Sys.time() - t0, units = "secs")
cor(fit_bm_known$est_raw$group_mean_mle, truth)

## ----fig.width = 4.5, fig.height = 4.5----------------------------------------
plot(truth, fit_bm_known$est_raw$group_mean_mle,
     pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.3),
     xlab = "True district mean",
     ylab = "Estimated mean (test-score scale)",
     main = "bin_means (known cuts)")
abline(0, 1, col = "red", lty = 2)

## -----------------------------------------------------------------------------
t0 <- Sys.time()
fit_bm_null <- bin_means(ngk)
t_bm_null   <- as.numeric(Sys.time() - t0, units = "secs")

## ----fig.width = 4.5, fig.height = 4.5----------------------------------------
plot(truth, fit_bm_null$est_std$group_mean_mle,
     pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.3),
     xlab = "True district mean",
     ylab = "Estimated mean (standardized)",
     main = "bin_means (cuts from data)")

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
# state_props <- colSums(ngk) / sum(ngk)
# state_cum   <- cumsum(state_props)
# 
# t0 <- Sys.time()
# fit_fh <- fh_hetop(
#   ngk       = ngk,
#   fixedcuts = qnorm(state_cum[1:2]),
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

