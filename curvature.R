## ============================================================
##  Simulation for curvature–calibrated confidence ellipsoids
## ============================================================

## Packages (install first if needed)
library(numDeriv)   # Jacobians / Hessians
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
library(MASS)   # for ginv if we decide to use it (optional)
set.seed(123)
## -----------------------------
##  Utility: ellipsoid volume
## -----------------------------
ellipsoid_volume <- function(V, c) {
  # Volume of { theta : (theta - mu)^T V^{-1} (theta - mu) <= c }
  p      <- ncol(V)
  detV   <- det(V)
  unit_V <- pi^(p/2) / gamma(p/2 + 1)  # volume of unit p-ball
  unit_V * (c^(p/2)) * sqrt(detV)
}

## ---------------------------------------
##  Design x for each model and n
## ---------------------------------------
design_x <- function(model_name, n) {
  switch(model_name,
         "Michaelis-Menten"  = seq(0.1, 10, length.out = n),
         "Exponential decay" = seq(0,    5,  length.out = n),
         "Logistic growth"   = seq(-2,   2,  length.out = n),
         seq(0, 1, length.out = n))
}

## ---------------------------------------
##  Models: mean functions + setups
## ---------------------------------------
models <- list(
  MM = list(
    name   = "Michaelis-Menten",
    p      = 2,
    f      = function(x, theta) {
      V <- theta[1]; K <- theta[2]
      V * x / (K + x)
    },
    formula      = y ~ V * x / (K + x),
    start        = list(V = 3, K = 1),
    par_names    = c("V", "K"),
    theta_configs = list(
      low    = c(5, 10),   # relatively linear
      medium = c(5,  2),
      high   = c(5,  0.5)  # high curvature
    )
  ),
  EXP = list(
    name   = "Exponential decay",
    p      = 2,
    f      = function(x, theta) {
      A <- theta[1]; B <- theta[2]
      A * exp(-B * x)
    },
    formula      = y ~ A * exp(-B * x),
    start        = list(A = 2, B = 0.2),
    par_names    = c("A", "B"),
    theta_configs = list(
      low    = c(3, 0.1),
      medium = c(3, 0.3),
      high   = c(3, 0.7)
    )
  ),
  LOG = list(
    name   = "Logistic growth",
    p      = 3,
    f      = function(x, theta) {
      L <- theta[1]; k <- theta[2]; x0 <- theta[3]
      L / (1 + exp(-k * (x - x0)))
    },
    formula      = y ~ L / (1 + exp(-k * (x - x0))),
    start        = list(L = 1, k = 1, x0 = 0),
    par_names    = c("L", "k", "x0"),
    theta_configs = list(
      low    = c(1, 1, 0),
      medium = c(1, 2, 0),
      high   = c(1, 3, 0)
    )
  )
)

## ------------------------------------------------
## Curvature factor kappa(θ): numerical proxy
## ------------------------------------------------
# We approximate curvature using Hessians of f(x_i, θ) for a
# small subset of design points and average their Frobenius norms.
curvature_factor <- function(theta, model, x, n_points = 10) {
  m <- length(x)
  # pick up to n_points roughly equally spaced indices
  idx <- seq(1, m, length.out = min(n_points, m))
  idx <- unique(round(idx))
  idx <- pmax(1, pmin(m, idx))  # clamp to [1, m]
  
  k2 <- 0
  for (i in idx) {
    Hi <- numDeriv::hessian(function(th) model$f(x[i], th), theta)
    k2 <- k2 + norm(Hi, type = "F")^2
  }
  sqrt(k2 / length(idx))
}



## ------------------------------------------------
## Jacobian J(θ) via numDeriv::jacobian
## ------------------------------------------------
compute_jacobian <- function(theta, model, x) {
  numDeriv::jacobian(function(th) model$f(x, th), theta)
}

## ------------------------------------------------
## Inflation function δ(κ, n)
## ------------------------------------------------
delta_fun <- function(kappa, n, max_delta = 0.05) {
  # Simple monotone function: δ ∝ κ / sqrt(n)
  pmin(max_delta, 2 * kappa / sqrt(n))
}

## ------------------------------------------------
## Run simulation for one configuration
## ------------------------------------------------
## ------------------------------------------------
## Run simulation for one configuration
## ------------------------------------------------
run_one_config <- function(model, theta0, n, M, alpha, regime_label,
                           sigma = 0.5, track_paths = FALSE) {
  x <- design_x(model$name, n)
  p <- length(theta0)
  
  # curvature at the true parameter (design-level diagnostic)
  kappa_true <- curvature_factor(theta0, model, x)
  
  coverage_wald <- coverage_corr <- coverage_lr <- rep(NA, M)
  vol_wald      <- vol_corr      <- vol_lr      <- rep(NA, M)
  fail_fit      <- 0L
  
  for (m in seq_len(M)) {
    ## simulate data
    y <- model$f(x, theta0) + rnorm(n, sd = sigma)
    dat <- data.frame(x = x, y = y)
    
    ## fit nonlinear model
    fit <- tryCatch(
      nls(model$formula,
          data    = dat,
          start   = model$start,
          control = nls.control(warnOnly = TRUE, maxiter = 100)),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      fail_fit <- fail_fit + 1L
      next
    }
    
    est <- coef(fit)[model$par_names]
    est <- as.numeric(est)
    
    ## Jacobian and XtX
    Jhat <- compute_jacobian(est, model, x)
    XtX  <- t(Jhat) %*% Jhat
    
    # Check for singular / ill-conditioned XtX
    if (any(!is.finite(XtX))) {
      fail_fit <- fail_fit + 1L
      next
    }
    if (rcond(XtX) < 1e-10) {
      # Too ill-conditioned; skip this replicate
      fail_fit <- fail_fit + 1L
      next
    }
    
    sig2_hat <- summary(fit)$sigma^2
    V_hat    <- sig2_hat * solve(XtX)
    
    ## Wald region
    Q      <- t(est - theta0) %*% solve(V_hat) %*% (est - theta0)
    c_wald <- qchisq(1 - alpha, df = p)
    coverage_wald[m] <- as.numeric(Q <= c_wald)
    vol_wald[m]      <- ellipsoid_volume(V_hat, c_wald)
    
    ## Corrected region
    kappa_hat  <- curvature_factor(est, model, x)
    delta      <- delta_fun(kappa_hat, n)
    alpha_star <- pmax(alpha - delta, 1e-3)
    c_star     <- qchisq(1 - alpha_star, df = p)
    coverage_corr[m] <- as.numeric(Q <= c_star)
    vol_corr[m]      <- ellipsoid_volume(V_hat, c_star)
    
    ## Likelihood-ratio region via profile CIs
    prof <- tryCatch(stats::profile(fit), error = function(e) NULL)
    if (!is.null(prof)) {
      ci_lr <- tryCatch(confint(prof, level = 1 - alpha),
                        error = function(e) NULL)
      if (!is.null(ci_lr)) {
        ci_sub <- ci_lr[model$par_names, , drop = FALSE]
        inside <- all(theta0 >= ci_sub[, 1] & theta0 <= ci_sub[, 2])
        coverage_lr[m] <- as.numeric(inside)
        # simple hyper-rectangle volume as a proxy for LR region size
        vol_lr[m]      <- prod(ci_sub[, 2] - ci_sub[, 1])
      }
    }
  }
  
  summary_tbl <- tibble(
    Model       = model$name,
    Regime      = regime_label,
    n           = n,
    p           = p,
    kappa_true  = kappa_true,
    Method      = c("Wald", "Corrected", "LR"),
    Coverage    = c(mean(coverage_wald, na.rm = TRUE),
                    mean(coverage_corr, na.rm = TRUE),
                    mean(coverage_lr,   na.rm = TRUE)),
    Volume      = c(mean(vol_wald, na.rm = TRUE),
                    mean(vol_corr, na.rm = TRUE),
                    mean(vol_lr,   na.rm = TRUE)),
    Failed_fits = fail_fit / M
  )
  
  if (!track_paths) {
    return(summary_tbl)
  }
  
  # Return also the per-replicate coverage indicators
  coverage_paths <- tibble(
    Iteration = 1:M,
    Wald      = coverage_wald,
    Corrected = coverage_corr,
    LR        = coverage_lr
  )
  
  list(summary = summary_tbl, paths = coverage_paths)
}


  


## ------------------------------------------------
## Master simulation over models / regimes / n
## ------------------------------------------------
run_full_simulation <- function(
    models, sample_sizes = c(30, 60, 120),
    M = 10000, alpha = 0.05, sigma = 0.5
) {
  res_list <- list()
  idx      <- 1L
  
  for (mname in names(models)) {
    model <- models[[mname]]
    for (regime in names(model$theta_configs)) {
      theta0 <- model$theta_configs[[regime]]
      for (n in sample_sizes) {
        cat("Running:", model$name,
            "| regime =", regime,
            "| n =", n, "\n")
        res_list[[idx]] <- run_one_config(
          model        = model,
          theta0       = theta0,
          n            = n,
          M            = M,
          alpha        = alpha,
          regime_label = regime,
          sigma        = sigma
        )
        idx <- idx + 1L
      }
    }
  }
  bind_rows(res_list)
}

## ------------------------------------------------
##   RUN THE SIMULATION
## ------------------------------------------------
set.seed(123)
sim_summary <- run_full_simulation(models,
                                   sample_sizes = c(30, 60),
                                   M = 10000,
                                   alpha = 0.05,
                                   sigma = 0.5)

## Look at first few rows
print(head(sim_summary))

## ------------------------------------------------
##  LaTeX table for the paper
## ------------------------------------------------
latex_table <- kable(
  sim_summary,
  format   = "latex",
  booktabs = TRUE,
  digits   = 3,
  caption  = "Empirical coverage, average region volume and failure rate by model, curvature regime, sample size and method.",
  label    = "tab:sim_results"
)

cat(latex_table)


plot_convergence <- function(paths, nominal = 0.95, title = "") {
  cov_long <- paths %>%
    mutate(
      Wald      = cumsum(tidyr::replace_na(Wald, 0)) / seq_len(n()),
      Corrected = cumsum(tidyr::replace_na(Corrected, 0)) / seq_len(n()),
      LR        = cumsum(tidyr::replace_na(LR, 0)) / seq_len(n())
    ) %>%
    tidyr::pivot_longer(
      cols      = c("Wald", "Corrected", "LR"),
      names_to  = "Method",
      values_to = "Coverage"
    )
  
  ggplot(cov_long, aes(x = Iteration, y = Coverage, color = Method)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = nominal, linetype = "dashed") +
    labs(
      title = title,
      x     = "Simulation replicate",
      y     = "Running empirical coverage"
    ) +
    theme_minimal()
}



set.seed(123)

res_exp_med60 <- run_one_config(
  model        = models$EXP,
  theta0       = models$EXP$theta_configs$medium,
  n            = 60,
  M            = 10000,
  alpha        = 0.05,
  regime_label = "medium",
  sigma        = 0.5,
  track_paths  = TRUE
)

res_exp_med60$summary
plot_convergence(
  res_exp_med60$paths,
  nominal = 0.95,
  title = "Exponential decay, medium curvature, n = 60"
)


set.seed(123)

res_log_high60 <- run_one_config(
  model        = models$LOG,
  theta0       = models$LOG$theta_configs$high,
  n            = 60,
  M            = 10000,
  alpha        = 0.05,
  regime_label = "high",
  sigma        = 0.5,
  track_paths  = TRUE
)

res_log_high60$summary
plot_convergence(
  res_log_high60$paths,
  nominal = 0.95,
  title = "Logistic growth, high curvature, n = 60"
)
