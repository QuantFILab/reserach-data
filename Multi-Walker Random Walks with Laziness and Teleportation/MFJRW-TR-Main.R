############################################################
## High-quality PNG saving (Windows Downloads folder)
############################################################
out_dir <- "C:/Users/USER/Downloads"  # <-- your target directory
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## Use ragg if available (best quality); otherwise fallback to default ggsave device
save_png_hq <- function(plot, filename,
                        width = 11, height = 7, units = "in",
                        dpi = 600, bg = "white") {
  full_path <- file.path(out_dir, filename)
  
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(
      filename = full_path,
      plot = plot,
      device = ragg::agg_png,
      width = width, height = height, units = units,
      dpi = dpi, bg = bg
    )
  } else {
    ggplot2::ggsave(
      filename = full_path,
      plot = plot,
      width = width, height = height, units = units,
      dpi = dpi, bg = bg
    )
  }
  message("Saved: ", full_path)
}

bump_fonts <- function(p, base = 16) {
  p + ggplot2::theme(
    text = ggplot2::element_text(size = base),
    axis.title = ggplot2::element_text(size = base + 1),
    axis.text  = ggplot2::element_text(size = base),
    strip.text = ggplot2::element_text(size = base),
    legend.title = ggplot2::element_text(size = base),
    legend.text  = ggplot2::element_text(size = base),
    plot.title   = ggplot2::element_text(size = base + 3, face = "bold"),
    plot.subtitle= ggplot2::element_text(size = base)
  )
}

############################################################
## Packages (NO viridis)
############################################################
pkgs <- c("ggplot2", "dplyr", "tidyr", "scales", "purrr")
need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need) > 0) install.packages(need)
library(ggplot2); library(dplyr); library(tidyr)
library(scales); library(purrr)

############################################################
## Global styling helpers (color-blind friendly; easy to see)
############################################################
# Okabe–Ito palette
OI <- c(
  black  = "#000000",
  blue   = "#0072B2",
  orange = "#D55E00",
  green  = "#009E73",
  yellow = "#E69F00",
  sky    = "#56B4E9",
  pink   = "#CC79A7",
  grey   = "#999999"
)

ctc_named_colors <- function(lvls) {
  lvls <- as.character(lvls)
  
  preferred <- c(
    "Monte Carlo" = OI["black"],
    "Theory"      = OI["blue"]
  )
  
  out <- setNames(rep(NA_character_, length(lvls)), lvls)
  hit <- intersect(names(preferred), lvls)
  out[hit] <- preferred[hit]
  
  remaining <- names(out)[is.na(out)]
  base <- c(OI["blue"], OI["orange"], OI["green"], OI["yellow"], OI["pink"], OI["sky"], OI["grey"])
  base <- base[!base %in% out]
  
  if (length(remaining) > length(base)) {
    base <- c(base, hcl.colors(length(remaining) - length(base), "Dark 3"))
  }
  out[remaining] <- base[seq_along(remaining)]
  out
}

ctc_d_colors <- c("50" = OI["blue"], "100" = OI["orange"], "200" = OI["green"])

ctc_rule_colors <- c(
  "Baseline (none)" = OI["blue"],
  "Blocking (sync)" = OI["orange"],
  "Bouncing (sync)" = OI["green"],
  "Pushing (seq)"   = OI["pink"]
)

ctc_fill_scale <- function(...) {
  scale_fill_gradientn(colours = hcl.colors(256, "YlOrRd"), ...)
}

############################################################
## Shared utilities
############################################################
mod_d <- function(x, d) ((x %% d) + d) %% d

################################################################################
################################################################################
##  SECTION 1 — Two-walker displacement-chain validation (SMALL CASE)
################################################################################
################################################################################

set.seed(123)

d     <- 50
alpha <- 0.05
s     <- 0.30

## forward-jump increment law (choose ONE fixed law)
mu_type <- "unit"  # "unit" | "uniformJ" | "geomJ"
J   <- 4
rho <- 0.60

## Arriving-distribution checks (panels)
n_grid        <- c(1, 5, 10, 25)
delta_arrive  <- c(1, floor(d/4), floor(d/2))

## Meeting-time checks
delta_meet <- 1:floor(d/2)

## Monte Carlo trials
N_arrive <- 120000   # per (delta,n)
N_meet   <- 50000    # per delta

## cap
n_max <- ceiling(200 * max(d^2, 1 / max(alpha, 1e-12)))

make_mu_forward_unit <- function() c(`1` = 1.0)
make_mu_forward_uniformJ <- function(J = 4) {
  inc <- 1:J
  p <- rep(1/J, J)
  stats::setNames(p, as.character(inc))
}
make_mu_forward_geomJ <- function(J = 6, rho = 0.6) {
  j <- 1:J
  w <- (1 - rho) * rho^(j - 1)
  w <- w / sum(w)
  stats::setNames(w, as.character(j))
}
sample_mu <- function(n, mu) {
  inc  <- as.integer(names(mu))
  prob <- as.numeric(mu)
  sample(inc, size = n, replace = TRUE, prob = prob)
}

mu <- switch(mu_type,
             unit     = make_mu_forward_unit(),
             uniformJ = make_mu_forward_uniformJ(J),
             geomJ    = make_mu_forward_geomJ(J, rho),
             make_mu_forward_unit()
)

update_positions <- function(x, d, alpha, s, mu) {
  N <- length(x)
  
  tele <- runif(N) < alpha
  if (any(tele)) {
    x[tele] <- sample.int(d, sum(tele), replace = TRUE) - 1L
  }
  
  not_tele <- !tele
  if (any(not_tele)) {
    lazy <- runif(N) < s
    lazy <- lazy & not_tele
    
    inc_idx <- not_tele & !lazy
    if (any(inc_idx)) {
      inc <- sample_mu(sum(inc_idx), mu)
      x[inc_idx] <- mod_d(x[inc_idx] + inc, d)
    }
  }
  x
}

simulate_displacement_at_n <- function(N, d, alpha, s, mu, delta, n) {
  x1 <- rep.int(as.integer(mod_d(delta, d)), N)
  x2 <- rep.int(0L, N)
  
  if (n >= 1) {
    for (t in seq_len(n)) {
      x1 <- update_positions(x1, d, alpha, s, mu)
      x2 <- update_positions(x2, d, alpha, s, mu)
    }
  }
  mod_d(x1 - x2, d)
}

simulate_meeting_time <- function(N, d, alpha, s, mu, delta, n_max) {
  x1 <- rep.int(as.integer(mod_d(delta, d)), N)
  x2 <- rep.int(0L, N)
  
  T <- rep.int(n_max, N)
  active <- rep.int(TRUE, N)
  
  for (n in seq_len(n_max)) {
    idx <- which(active)
    if (length(idx) == 0) break
    
    x1[idx] <- update_positions(x1[idx], d, alpha, s, mu)
    x2[idx] <- update_positions(x2[idx], d, alpha, s, mu)
    
    hit <- (x1[idx] == x2[idx])
    if (any(hit)) {
      hit_ids <- idx[hit]
      T[hit_ids] <- n
      active[hit_ids] <- FALSE
    }
  }
  
  list(T = T, censored = active)
}

mu_hat <- function(d, mu) {
  j <- as.integer(names(mu))
  p <- as.numeric(mu)
  k <- 0:(d-1)
  E <- exp(1i * 2*pi/d * outer(k, j, "*"))
  as.vector(E %*% p)
}
lambda_k <- function(d, alpha, s, mu) {
  mh <- mu_hat(d, mu)
  lam <- complex(length = d)
  lam[1] <- 1 + 0i
  if (d >= 2) lam[2:d] <- (1 - alpha) * (s + (1 - s) * mh[2:d])
  lam
}
gamma_k <- function(d, alpha, s, mu) {
  lam <- lambda_k(d, alpha, s, mu)
  Mod(lam)^2
}
arrive_diff_theory <- function(d, alpha, s, mu, delta, n) {
  G <- gamma_k(d, alpha, s, mu)
  r <- 0:(d-1)
  k <- 1:(d-1)
  
  phase <- exp(-1i * 2*pi/d * outer(r - delta, k, "*"))
  series <- phase %*% (G[k+1]^n)
  
  a <- (1/d) + (1/d) * series
  a <- pmax(Re(a), 0)
  a <- a / sum(a)
  as.numeric(a)
}
mean_meeting_theory <- function(d, alpha, s, mu, delta) {
  G <- gamma_k(d, alpha, s, mu)
  k <- 1:(d-1)
  omega_term <- exp(1i * 2*pi/d * (delta * k))
  term <- (1 - omega_term) * (G[k+1] / (1 - G[k+1]))
  Re(d + sum(term))
}
var_meeting_theory <- function(d, alpha, s, mu, delta) {
  G <- gamma_k(d, alpha, s, mu)
  k <- 1:(d-1)
  omega_term <- exp(1i * 2*pi/d * (delta * k))
  mu_delta <- mean_meeting_theory(d, alpha, s, mu, delta)
  
  A <- 2 * sum((1 - omega_term) * (G[k+1] / (1 - G[k+1])^2))
  B <- (d - 1) + sum((1 + omega_term) * (G[k+1] / (1 - G[k+1])))
  Re(A + mu_delta * B)
}

############################################################
## PLOT 1) Arriving distribution panels (MC=dots, Theory=solid)
############################################################
arrive_rows <- list()
idx <- 1
for (delta in delta_arrive) {
  for (n in n_grid) {
    Dn <- simulate_displacement_at_n(N_arrive, d, alpha, s, mu, delta, n)
    emp <- tabulate(Dn + 1, nbins = d) / length(Dn)
    th  <- arrive_diff_theory(d, alpha, s, mu, delta, n)
    
    arrive_rows[[idx]] <- tibble(
      r = 0:(d-1),
      prob_emp = emp,
      prob_th  = th,
      n = factor(n, levels = n_grid),
      delta = factor(delta, levels = delta_arrive)
    )
    idx <- idx + 1
  }
}
arrive_df <- bind_rows(arrive_rows)

arrive_long <- arrive_df %>%
  pivot_longer(c(prob_emp, prob_th), names_to = "source", values_to = "prob") %>%
  mutate(source = recode(source, prob_emp = "Monte Carlo", prob_th = "Theory"),
         source = factor(source, levels = c("Monte Carlo","Theory")))

p1 <- ggplot(arrive_long, aes(x = r, y = prob)) +
  geom_line(
    data = subset(arrive_long, source == "Theory"),
    aes(color = source),
    linewidth = 0.9, alpha = 0.80, linetype = "solid"
  ) +
  geom_point(
    data = subset(arrive_long, source == "Monte Carlo"),
    aes(color = source),
    size = 0.95, alpha = 0.55
  ) +
  facet_grid(delta ~ n,
             labeller = labeller(delta = function(x) paste0("delta=", x),
                                 n = function(x) paste0("n=", x))) +
  scale_color_manual(values = ctc_named_colors(levels(arrive_long$source))) +
  labs(x = "r (displacement)", y = "probability", color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
p1 <- bump_fonts(p1, base = 16)
print(p1)
save_png_hq(p1, "two_walker_arriving_panels.png", width = 12, height = 8, dpi = 600)

############################################################
## PLOT 2) Meeting-time mean/variance vs delta (MC=dots, Theory=solid)
############################################################
meet_rows <- vector("list", length(delta_meet))
for (ii in seq_along(delta_meet)) {
  delta <- delta_meet[ii]
  
  out <- simulate_meeting_time(N_meet, d, alpha, s, mu, delta, n_max)
  T <- out$T
  cens_rate <- mean(out$censored)
  
  meet_rows[[ii]] <- tibble(
    delta = delta,
    cens_rate = cens_rate,
    mc_mean = mean(T),
    mc_var  = var(T),
    th_mean = mean_meeting_theory(d, alpha, s, mu, delta),
    th_var  = var_meeting_theory(d, alpha, s, mu, delta)
  )
}
meet_stats <- bind_rows(meet_rows)
cat(sprintf("Max censoring rate over deltas: %.4f\n", max(meet_stats$cens_rate)))

meet_long <- meet_stats %>%
  select(delta, mc_mean, th_mean, mc_var, th_var) %>%
  pivot_longer(-delta, names_to = "key", values_to = "value") %>%
  mutate(
    stat = ifelse(grepl("mean", key), "Mean  E[T]", "Variance  Var(T)"),
    source = ifelse(grepl("^mc_", key), "Monte Carlo", "Theory"),
    source = factor(source, levels = c("Monte Carlo","Theory"))
  )

p2 <- ggplot(meet_long, aes(x = delta, y = value)) +
  geom_line(
    data = subset(meet_long, source == "Theory"),
    aes(color = source),
    linewidth = 0.95, alpha = 0.80, linetype = "solid"
  ) +
  geom_point(
    data = subset(meet_long, source == "Monte Carlo"),
    aes(color = source),
    size = 1.15, alpha = 0.55
  ) +
  facet_wrap(~ stat, scales = "free_y", ncol = 2) +
  scale_color_manual(values = ctc_named_colors(levels(meet_long$source))) +
  labs(x = "initial displacement delta", y = "value", color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
p2 <- bump_fonts(p2, base = 16)
print(p2)
save_png_hq(p2, "two_walker_meeting_moments.png", width = 11, height = 6, dpi = 600)

############################################################
## PLOT 3) Heatmap with value labels (fill easy to see)
############################################################
plot_small_heatmap_labeled_facet <- function(M, V, title,
                                             alpha_tile = 0.90, digits = 1) {
  make_df <- function(A, stat_name) {
    as.data.frame(A) %>%
      mutate(i = row_number()) %>%
      pivot_longer(-i, names_to = "j", values_to = "val") %>%
      mutate(j = as.integer(gsub("V", "", j)),
             stat = stat_name)
  }
  
  df <- bind_rows(make_df(M, "Mean  E[T]"),
                  make_df(V, "Variance  Var(T)"))
  
  ggplot(df, aes(x = j, y = i, fill = val)) +
    geom_tile(alpha = alpha_tile, color = NA) +
    geom_text(aes(label = format(round(val, digits), nsmall = digits)),
              size = 4.2, alpha = 0.92, color = OI["black"]) +
    facet_wrap(~ stat, nrow = 1) +
    coord_equal() +
    scale_y_reverse(breaks = 1:nrow(M)) +
    scale_x_continuous(breaks = 1:ncol(M)) +
    ctc_fill_scale() +
    labs(title = title, x = "walker j", y = "walker i", fill = "value") +
    theme_minimal(base_size = 12)
}

delta0 <- 1
out0 <- simulate_meeting_time(N_meet, d, alpha, s, mu, delta0, n_max)
m0 <- mean(out0$T); v0 <- var(out0$T)

Mhat <- matrix(c(0, m0, m0, 0), nrow = 2, byrow = TRUE)
Vhat <- matrix(c(0, v0, v0, 0), nrow = 2, byrow = TRUE)

p3 <- plot_small_heatmap_labeled_facet(
  Mhat, Vhat,
  title ="",
  digits = 1
)
p3 <- bump_fonts(p3, base = 16)
print(p3)
save_png_hq(p3, "two_walker_heatmap_mean_var.png", width = 11, height = 4.5, dpi = 600)

################################################################################
################################################################################
##  SECTION 2 — Many-walker first-collision validation (SMALL GRID)
################################################################################
################################################################################

set.seed(123)

alpha <- 0.05
s     <- 0.30

d_grid <- c(50, 100, 200)
L_grid <- c(4, 8, 16)

N_coll <- 30000
N_pn   <- 60000

n_tail_max_cap <- 400
n_tail_step    <- 2

d_ex <- 50
L_ex <- 8
n_ex_grid <- c(1, 2, 5, 10, 25, 50, 100)

mu <- c(`1` = 1.0)

sample_mu <- function(n, mu) {
  inc  <- as.integer(names(mu))
  prob <- as.numeric(mu)
  sample(inc, size = n, replace = TRUE, prob = prob)
}

pi_collision <- function(d, L) {
  if (L > d) return(1)
  log_no_coll <- sum(log(d - (0:(L-1)))) - L * log(d)
  1 - exp(log_no_coll)
}

has_collision_rows <- function(X) {
  apply(X, 1, function(r) anyDuplicated(r) > 0)
}

update_positions_mat <- function(X, d, alpha, s, mu) {
  N <- nrow(X); L <- ncol(X)
  
  tele <- matrix(runif(N*L) < alpha, nrow = N, ncol = L)
  if (any(tele)) {
    X[tele] <- sample.int(d, sum(tele), replace = TRUE) - 1L
  }
  
  not_tele <- !tele
  if (any(not_tele)) {
    lazy <- matrix(runif(N*L) < s, nrow = N, ncol = L)
    lazy <- lazy & not_tele
    
    inc_idx <- not_tele & !lazy
    if (any(inc_idx)) {
      if (length(mu) == 1 && names(mu)[1] == "1") {
        X[inc_idx] <- mod_d(X[inc_idx] + 1L, d)
      } else {
        inc <- sample_mu(sum(inc_idx), mu)
        X[inc_idx] <- mod_d(X[inc_idx] + inc, d)
      }
    }
  }
  X
}

init_positions <- function(N, d, L, init = c("clustered_shifted", "clustered_fixed", "iid_with_replacement")) {
  init <- match.arg(init)
  if (init == "iid_with_replacement") {
    X <- matrix(sample.int(d, N*L, replace = TRUE) - 1L, nrow = N, ncol = L)
    return(X)
  }
  base <- 0:(L-1)
  if (init == "clustered_fixed") {
    X <- matrix(rep(base, each = N), nrow = N, ncol = L, byrow = FALSE)
    X <- apply(X, 2, function(col) mod_d(col, d))
    return(X)
  }
  shift <- sample.int(d, N, replace = TRUE) - 1L
  X <- sapply(base, function(b) mod_d(shift + b, d))
  X
}

simulate_Tcoll_many <- function(N_rep, d, L, alpha, s, mu,
                                init = "clustered_shifted",
                                n_max = NULL) {
  
  if (is.null(n_max)) {
    pC <- pi_collision(d, L)
    n_max <- min(5000, ceiling(30 / max(pC, 1e-6)))
  }
  
  X <- init_positions(N_rep, d, L, init = init)
  
  T <- rep.int(n_max, N_rep)
  active <- rep.int(TRUE, N_rep)
  
  for (n in seq_len(n_max)) {
    idx <- which(active)
    if (length(idx) == 0) break
    
    X[idx, ] <- update_positions_mat(X[idx, , drop = FALSE], d, alpha, s, mu)
    
    coll <- has_collision_rows(X[idx, , drop = FALSE])
    if (any(coll)) {
      hit_ids <- idx[coll]
      T[hit_ids] <- n
      active[hit_ids] <- FALSE
    }
  }
  
  list(T = T, censored = active, n_max = n_max)
}

## Pairwise theory for union bound plot 3
mu_hat <- function(d, mu) {
  j <- as.integer(names(mu))
  p <- as.numeric(mu)
  k <- 0:(d-1)
  E <- exp(1i * 2*pi/d * outer(k, j, "*"))
  as.vector(E %*% p)
}
lambda_k <- function(d, alpha, s, mu) {
  mh <- mu_hat(d, mu)
  lam <- complex(length = d)
  lam[1] <- 1 + 0i
  if (d >= 2) lam[2:d] <- (1 - alpha) * (s + (1 - s) * mh[2:d])
  lam
}
gamma_k <- function(d, alpha, s, mu) {
  lam <- lambda_k(d, alpha, s, mu)
  Mod(lam)^2
}
p_pair_meet_at_time_n <- function(d, alpha, s, mu, delta, n) {
  G <- gamma_k(d, alpha, s, mu)
  k <- 1:(d-1)
  phase <- exp(1i * 2*pi/d * (delta * k))
  p0 <- (1/d) + (1/d) * sum(phase * (G[k+1]^n))
  Re(p0)
}
union_upper_bound_collision_prob <- function(d, L, alpha, s, mu, n) {
  ub <- 0
  for (i in 1:(L-1)) {
    for (j in (i+1):L) {
      delta_ij <- mod_d((i-1) - (j-1), d)
      ub <- ub + p_pair_meet_at_time_n(d, alpha, s, mu, delta_ij, n)
    }
  }
  ub
}

############################################################
## PLOT 1) Tail: MC dots vs Geometric approx solid (facets)
############################################################
grid <- expand.grid(d = d_grid, L = L_grid) %>%
  dplyr::filter(L <= d) %>%
  as_tibble()

tail_rows <- list()
tail_idx <- 1

cat("Simulating T_coll for (d,L) grid...\n")
for (rr in seq_len(nrow(grid))) {
  d0 <- grid$d[rr]
  L0 <- grid$L[rr]
  pC <- pi_collision(d0, L0)
  
  out <- simulate_Tcoll_many(N_coll, d0, L0, alpha, s, mu, init = "clustered_shifted")
  T <- out$T
  
  n_plot_max <- min(n_tail_max_cap, max(20, ceiling(quantile(T, 0.95))))
  n_eval <- seq(0, n_plot_max, by = n_tail_step)
  
  S_mc <- sapply(n_eval, function(n) mean(T > n))
  S_th <- (1 - pC)^n_eval
  
  tail_rows[[tail_idx]] <- tibble(
    d = d0, L = L0,
    n = n_eval,
    S_mc = S_mc,
    S_th = S_th,
    piC = pC
  )
  tail_idx <- tail_idx + 1
  
  cat(sprintf("  done (d=%d, L=%d): pi(C)=%.4g, n_max used=%d\n", d0, L0, pC, out$n_max))
}

tail_df <- bind_rows(tail_rows) %>%
  pivot_longer(c(S_mc, S_th), names_to = "source", values_to = "S") %>%
  mutate(source = recode(source, S_mc = "Monte Carlo", S_th = "Geometric approx"),
         source = factor(source, levels = c("Monte Carlo","Geometric approx")))

p1 <- ggplot(tail_df, aes(x = n, y = S)) +
  geom_line(
    data = subset(tail_df, source == "Geometric approx"),
    aes(color = source),
    linewidth = 0.9, alpha = 0.80, linetype = "solid"
  ) +
  geom_point(
    data = subset(tail_df, source == "Monte Carlo"),
    aes(color = source),
    size = 1.0, alpha = 0.55
  ) +
  facet_grid(d ~ L, labeller = labeller(d = function(x) paste0("d=", x),
                                        L = function(x) paste0("L=", x))) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(values = ctc_named_colors(levels(tail_df$source))) +
  labs(x = "n", y = "P(T_coll > n)", color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
p1 <- bump_fonts(p1, base = 16)
print(p1)
save_png_hq(p1, "mw_tail_facets.png", width = 12, height = 8, dpi = 600)

############################################################
## PLOT 2) Mean E[T_coll] vs L: MC dots; 1/pi(C) solid; heuristic dashed
############################################################
mean_rows <- list()
idx <- 1
for (rr in seq_len(nrow(grid))) {
  d0 <- grid$d[rr]
  L0 <- grid$L[rr]
  out <- simulate_Tcoll_many(N_coll, d0, L0, alpha, s, mu, init = "clustered_shifted")
  T <- out$T
  pC <- pi_collision(d0, L0)
  mean_rows[[idx]] <- tibble(
    d = d0, L = L0,
    mc_mean = mean(T),
    th_mean = 1 / pC,
    heur_mean = 2*d0 / (L0*(L0-1))
  )
  idx <- idx + 1
}
mean_df <- bind_rows(mean_rows)

mean_long <- mean_df %>%
  pivot_longer(c(mc_mean, th_mean, heur_mean), names_to = "source", values_to = "value") %>%
  mutate(source = recode(source,
                         mc_mean = "Monte Carlo",
                         th_mean = "1/pi(C)",
                         heur_mean = "2d/(L(L-1))"),
         source = factor(source, levels = c("Monte Carlo","1/pi(C)","2d/(L(L-1))")),
         d = factor(d))

p2 <- ggplot(mean_long, aes(x = L, y = value, group = interaction(d, source))) +
  geom_line(
    data = subset(mean_long, source != "Monte Carlo"),
    aes(color = d, linetype = source),
    linewidth = 0.9, alpha = 0.80
  ) +
  geom_point(
    data = subset(mean_long, source == "Monte Carlo"),
    aes(color = d),
    size = 2.0, alpha = 0.55
  ) +
  scale_color_manual(values = ctc_d_colors) +
  scale_linetype_manual(values = c("1/pi(C)" = "solid", "2d/(L(L-1))" = "dashed")) +
  labs(x = "L", y = "E[T_coll]", color = "d", linetype = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
p2 <- bump_fonts(p2, base = 16)
print(p2)
save_png_hq(p2, "mw_mean_scaling.png", width = 11, height = 6.5, dpi = 600)

############################################################
## PLOT 3) P(R_n in C): MC dots vs UNION upper bound dashed (+ pi(C) dashed)
############################################################
cat(sprintf("\nPlot 3: estimating P(R_n in C) for (d=%d, L=%d) clustered start...\n", d_ex, L_ex))

p_emp <- numeric(length(n_ex_grid))
for (tt in seq_along(n_ex_grid)) {
  n0 <- n_ex_grid[tt]
  X <- init_positions(N_pn, d_ex, L_ex, init = "clustered_fixed")
  for (k in seq_len(n0)) {
    X <- update_positions_mat(X, d_ex, alpha, s, mu)
  }
  p_emp[tt] <- mean(has_collision_rows(X))
  cat(sprintf("  n=%d done\n", n0))
}

p_ub <- sapply(n_ex_grid, function(n0)
  union_upper_bound_collision_prob(d_ex, L_ex, alpha, s, mu, n0)
)

p_pi <- pi_collision(d_ex, L_ex)

df3 <- tibble(
  n = n_ex_grid,
  `Monte Carlo` = p_emp,
  `Union bound` = p_ub
) %>%
  pivot_longer(-n, names_to = "series", values_to = "prob") %>%
  mutate(series = factor(series, levels = c("Monte Carlo","Union bound")))

p3 <- ggplot(df3, aes(x = n, y = prob, color = series)) +
  geom_line(
    data = subset(df3, series == "Union bound"),
    linewidth = 0.9, alpha = 0.80, linetype = "dashed"
  ) +
  geom_point(
    data = subset(df3, series == "Monte Carlo"),
    size = 2.0, alpha = 0.55
  ) +
  geom_hline(yintercept = p_pi, linetype = "dashed", linewidth = 0.9, alpha = 0.65) +
  annotate("text", x = max(n_ex_grid), y = p_pi, label = "pi(C)",
           hjust = 1.1, vjust = -0.4, size = 3.5, alpha = 0.8) +
  scale_color_manual(values = ctc_named_colors(levels(df3$series))) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "n", y = "Probability", color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
p3 <- bump_fonts(p3, base = 16)
print(p3)
save_png_hq(p3, "mw_collision_prob_union_bound.png", width = 11, height = 6.5, dpi = 600)

################################################################################
################################################################################
##  SECTION 3 — Coalescence-time validation (FORWARD-JUMP on cycle Z_d)
################################################################################
################################################################################

set.seed(123)

alpha <- 0.05
s     <- 0.30
mu <- c(`1` = 1.0)

d_grid <- c(50, 100, 200)
L_grid <- c(4, 8, 16)

d_all   <- 50
L_all   <- c(2, 3, 4)

N_coal_grid   <- 12000
N_tau_profile <- 15000
N_all         <- 3000

nmax_coal <- function(d) ceiling(500 * d)
nmax_all  <- 2e6

falling_fact <- function(d, k) {
  if (k == 0) return(1)
  if (k > d) return(0)
  prod(d:(d - k + 1))
}
pi_collision_k <- function(d, k) {
  if (k <= 1) return(0)
  1 - falling_fact(d, k) / (d^k)
}
pi_all_equal <- function(d, L) {
  d^(-(L - 1))
}

step_many <- function(pos, d, alpha, s, mu) {
  K <- length(pos)
  y <- pos
  
  tele <- runif(K) < alpha
  if (any(tele)) y[tele] <- sample.int(d, sum(tele), replace = TRUE) - 1L
  
  not_tele <- !tele
  if (any(not_tele)) {
    lazy <- (runif(K) < s) & not_tele
    inc_idx <- not_tele & !lazy
    if (any(inc_idx)) {
      inc <- sample_mu(sum(inc_idx), mu)
      y[inc_idx] <- mod_d(pos[inc_idx] + inc, d)
    }
  }
  y
}

simulate_T_all_one <- function(d, L, alpha, s, mu, n_max) {
  pos <- sample.int(d, L, replace = FALSE) - 1L
  for (n in seq_len(n_max)) {
    pos <- step_many(pos, d, alpha, s, mu)
    if (all(pos == pos[1])) return(list(Tall = n, censored = FALSE))
  }
  list(Tall = n_max, censored = TRUE)
}
simulate_T_all_many <- function(N, d, L, alpha, s, mu, n_max) {
  Tall <- integer(N)
  cens <- logical(N)
  for (m in seq_len(N)) {
    out <- simulate_T_all_one(d, L, alpha, s, mu, n_max)
    Tall[m] <- out$Tall
    cens[m] <- out$censored
  }
  list(Tall = Tall, censored = cens)
}

simulate_coalescing_one <- function(d, L, alpha, s, mu, n_max,
                                    init = c("iid_no_collision","clustered_shifted")) {
  init <- match.arg(init)
  
  cid <- seq_len(L)
  
  if (init == "iid_no_collision") {
    pos <- sample.int(d, L, replace = FALSE) - 1L
  } else {
    U <- sample.int(d, 1) - 1L
    pos <- mod_d(U + (0:(L-1)), d)
  }
  
  tau <- rep(NA_integer_, L)
  names(tau) <- as.character(seq_len(L))
  
  active <- sort(unique(cid))
  reps <- sapply(active, function(a) pos[which(cid == a)[1]])
  grp <- split(active, reps)
  if (any(lengths(grp) >= 2)) {
    for (g in grp) {
      if (length(g) >= 2) {
        new_id <- min(g)
        for (old in g) cid[cid == old] <- new_id
      }
    }
  }
  
  t_last <- 0L
  
  for (n in seq_len(n_max)) {
    active <- sort(unique(cid))
    K_old <- length(active)
    if (K_old == 1) return(list(Tcoal = n-1L, tau = tau, censored = FALSE))
    
    pos_active <- sapply(active, function(a) pos[which(cid == a)[1]])
    pos_active_new <- step_many(pos_active, d, alpha, s, mu)
    
    for (idx in seq_along(active)) {
      a <- active[idx]
      pos[cid == a] <- pos_active_new[idx]
    }
    
    reps2 <- sapply(active, function(a) pos[which(cid == a)[1]])
    grp2 <- split(active, reps2)
    if (any(lengths(grp2) >= 2)) {
      for (g in grp2) {
        if (length(g) >= 2) {
          new_id <- min(g)
          for (old in g) cid[cid == old] <- new_id
        }
      }
    }
    
    K_new <- length(unique(cid))
    
    if (K_new < K_old) {
      dt <- n - t_last
      tau[as.character(K_old)] <- dt
      if (K_old - K_new >= 2) {
        for (kk in (K_old-1):(K_new+1)) tau[as.character(kk)] <- 0L
      }
      t_last <- n
    }
  }
  
  list(Tcoal = n_max, tau = tau, censored = TRUE)
}

simulate_coalescing_many <- function(N, d, L, alpha, s, mu, n_max,
                                     init = "iid_no_collision") {
  Tcoal <- integer(N)
  cens  <- logical(N)
  tau_mat <- matrix(NA_integer_, nrow = N, ncol = L,
                    dimnames = list(NULL, as.character(seq_len(L))))
  for (m in seq_len(N)) {
    out <- simulate_coalescing_one(d, L, alpha, s, mu, n_max, init = init)
    Tcoal[m] <- out$Tcoal
    cens[m]  <- out$censored
    tau_mat[m, ] <- out$tau
  }
  list(Tcoal = Tcoal, censored = cens, tau_mat = tau_mat)
}

coal_grid <- expand.grid(d = d_grid, L = L_grid) %>% as_tibble()

coal_summ <- pmap_dfr(coal_grid, function(d, L) {
  out <- simulate_coalescing_many(N_coal_grid, d, L, alpha, s, mu, nmax_coal(d), init = "iid_no_collision")
  tibble(
    d = d, L = L,
    mc_mean = mean(out$Tcoal),
    mc_var  = var(out$Tcoal),
    cens_rate = mean(out$censored),
    th_mean_geo = sum(sapply(2:L, function(k) 1 / pi_collision_k(d, k)))
  )
})

d0 <- 50
L0 <- 16
out_prof <- simulate_coalescing_many(N_tau_profile, d0, L0, alpha, s, mu, nmax_coal(d0),
                                     init = "iid_no_collision")

tau_df <- as.data.frame(out_prof$tau_mat) %>%
  mutate(rep = row_number()) %>%
  pivot_longer(-rep, names_to = "k", values_to = "tau") %>%
  mutate(k = as.integer(k)) %>%
  filter(k >= 2) %>%
  group_by(k) %>%
  summarise(mc_tau_mean = mean(tau, na.rm = TRUE),
            mc_tau_var  = var(tau,  na.rm = TRUE),
            .groups = "drop") %>%
  mutate(th_tau_geo = 1 / sapply(k, function(kk) pi_collision_k(d0, kk)))

all_summ <- map_dfr(L_all, function(L) {
  out <- simulate_T_all_many(N_all, d_all, L, alpha, s, mu, nmax_all)
  tibble(
    d = d_all, L = L,
    mc_mean = mean(out$Tall),
    mc_var  = var(out$Tall),
    cens_rate = mean(out$censored),
    th_mean_geo = 1 / pi_all_equal(d_all, L)
  )
})

## Plot B1: E[T_coal] vs L
coal_long <- coal_summ %>%
  select(d, L, mc_mean, th_mean_geo) %>%
  pivot_longer(c(mc_mean, th_mean_geo), names_to = "source", values_to = "value") %>%
  mutate(source = recode(source, mc_mean = "Monte Carlo", th_mean_geo = "Heuristic  sum 1/pi_k(C_k)"),
         source = factor(source, levels = c("Monte Carlo","Heuristic  sum 1/pi_k(C_k)")))

p_coal_mean <- ggplot(coal_long, aes(x = L, y = value, color = source)) +
  geom_line(data = subset(coal_long, source != "Monte Carlo"),
            linewidth = 0.9, alpha = 0.80, linetype = "solid") +
  geom_point(data = subset(coal_long, source == "Monte Carlo"),
             size = 1.8, alpha = 0.55) +
  facet_wrap(~ d, scales = "free_y",
             labeller = labeller(d = function(x) paste0("d=", x))) +
  scale_color_manual(values = ctc_named_colors(levels(coal_long$source))) +
  labs(x = "number of walkers L", y = "mean time", color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
p_coal_mean <- bump_fonts(p_coal_mean, base = 16)
print(p_coal_mean)
save_png_hq(p_coal_mean, "coal_mean_vs_L.png", width = 11, height = 7, dpi = 600)

## Plot B2: merge-time profile E[tau_k]
tau_long <- tau_df %>%
  select(k, mc_tau_mean, th_tau_geo) %>%
  pivot_longer(c(mc_tau_mean, th_tau_geo), names_to = "source", values_to = "value") %>%
  mutate(source = recode(source, mc_tau_mean = "Monte Carlo", th_tau_geo = "Heuristic  1/pi_k(C_k)"),
         source = factor(source, levels = c("Monte Carlo","Heuristic  1/pi_k(C_k)")))

p_tau <- ggplot(tau_long, aes(x = k, y = value, color = source)) +
  geom_line(data = subset(tau_long, source != "Monte Carlo"),
            linewidth = 0.9, alpha = 0.80, linetype = "solid") +
  geom_point(data = subset(tau_long, source == "Monte Carlo"),
             size = 1.8, alpha = 0.55) +
  scale_x_reverse(breaks = sort(unique(tau_long$k), decreasing = TRUE)) +
  scale_color_manual(values = ctc_named_colors(levels(tau_long$source))) +
  labs(x = "number of clusters k (decreases to 1)", y = "mean merge time", color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
p_tau <- bump_fonts(p_tau, base = 16)
print(p_tau)
save_png_hq(p_tau, "coal_merge_profile_tau_k.png", width = 11, height = 6.5, dpi = 600)

## Plot A: all-equal diagonal time E[T_all]
all_long <- all_summ %>%
  select(L, mc_mean, th_mean_geo) %>%
  pivot_longer(c(mc_mean, th_mean_geo), names_to = "source", values_to = "value") %>%
  mutate(source = recode(source, mc_mean = "Monte Carlo", th_mean_geo = "Heuristic  d^(L-1)"),
         source = factor(source, levels = c("Monte Carlo","Heuristic  d^(L-1)")))

p_all <- ggplot(all_long, aes(x = L, y = value, color = source)) +
  geom_line(data = subset(all_long, source != "Monte Carlo"),
            linewidth = 0.9, alpha = 0.80, linetype = "solid") +
  geom_point(data = subset(all_long, source == "Monte Carlo"),
             size = 2.0, alpha = 0.55) +
  scale_y_log10(labels = label_number()) +
  scale_color_manual(values = ctc_named_colors(levels(all_long$source))) +
  labs(x = "number of walkers L", y = "mean time (log10 scale)", color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
p_all <- bump_fonts(p_all, base = 16)
print(p_all)
save_png_hq(p_all, "coal_all_equal_Tall.png", width = 11, height = 6.5, dpi = 600)

cat("\nCensoring rates (T_coal grid):\n")
print(coal_summ %>% select(d, L, cens_rate))
cat("\nCensoring rates (T_all):\n")
print(all_summ %>% select(d, L, cens_rate))

################################################################################
################################################################################
##  SECTION 4 — Teleportation & Collision Frequency — Validation Code
##  (3 figures: p1, p2, p3)
################################################################################
################################################################################

set.seed(123)

s_fixed <- 0.30

pi_collision <- function(d, L) {
  if (L > d) return(1)
  log_falling <- lgamma(d + 1) - lgamma(d - L + 1)
  1 - exp(log_falling - L * log(d))
}
pair_deltas_clustered <- function(L, d) {
  ij <- t(combn(L, 2))
  mod_d(ij[,1] - ij[,2], d)
}

mu_hat_unit <- function(d) {
  k <- 0:(d-1)
  exp(1i * 2*pi/d * k)
}
lambda_k_unit <- function(d, alpha, s) {
  mh <- mu_hat_unit(d)
  lam <- complex(length = d)
  lam[1] <- 1 + 0i
  if (d >= 2) lam[2:d] <- (1 - alpha) * (s + (1 - s) * mh[2:d])
  lam
}
gamma_k_unit <- function(d, alpha, s) {
  lam <- lambda_k_unit(d, alpha, s)
  Mod(lam)^2
}
p_pair_equal_theory <- function(d, alpha, s, delta, n) {
  G <- gamma_k_unit(d, alpha, s)
  k <- 1:(d-1)
  omega_term <- exp(1i * 2*pi/d * (delta * k))
  val <- (1/d) + (1/d) * sum(omega_term * (G[k+1]^n))
  Re(val)
}
union_bound_collision_time_n <- function(d, L, alpha, s, n) {
  deltas <- pair_deltas_clustered(L, d)
  ub <- sum(vapply(deltas, function(delta) p_pair_equal_theory(d, alpha, s, delta, n), numeric(1)))
  pmin(1, ub)
}

tv_bound_independent <- function(L, alpha, n) {
  rho <- (1 - alpha)^n
  1 - (1 - rho)^L
}
tv_bound_shared_clock <- function(alpha, n) (1 - alpha)^n
lower_bound_collision_prob <- function(piC, tvb) pmax(0, piC - tvb)

init_positions <- function(N, d, L, init = c("iid_uniform", "clustered")) {
  init <- match.arg(init)
  if (init == "iid_uniform") {
    matrix(sample.int(d, N * L, replace = TRUE) - 1L, nrow = N, ncol = L)
  } else {
    base <- mod_d(0:(L-1), d)
    matrix(rep(base, each = N), nrow = N, ncol = L, byrow = FALSE)
  }
}
collision_indicator <- function(pos_mat) {
  apply(pos_mat, 1, function(r) anyDuplicated(r) > 0L)
}
step_batch <- function(pos, d, alpha, s, coupling = c("independent", "shared_clock")) {
  coupling <- match.arg(coupling)
  if (is.na(alpha)) stop("alpha is NA; check your parameter grid.")
  
  N <- nrow(pos); L <- ncol(pos)
  
  if (coupling == "independent") {
    tele <- matrix(runif(N * L) < alpha, nrow = N, ncol = L)
    if (any(tele)) pos[tele] <- sample.int(d, sum(tele), replace = TRUE) - 1L
    
    not_tele <- !tele
    lazy <- matrix(runif(N * L) < s, nrow = N, ncol = L) & not_tele
    inc  <- not_tele & !lazy
    if (any(inc)) pos[inc] <- mod_d(pos[inc] + 1L, d)
    return(pos)
  }
  
  B <- runif(N) < alpha
  if (any(B)) {
    idx <- which(B)
    pos[idx, ] <- matrix(sample.int(d, length(idx) * L, replace = TRUE) - 1L,
                         nrow = length(idx), ncol = L)
  }
  if (any(!B)) {
    idx <- which(!B)
    lazy <- matrix(runif(length(idx) * L) < s, nrow = length(idx), ncol = L)
    inc  <- !lazy
    if (any(inc)) {
      sub <- pos[idx, , drop = FALSE]
      sub[inc] <- mod_d(sub[inc] + 1L, d)
      pos[idx, ] <- sub
    }
  }
  pos
}

simulate_collision_prob_timegrid <- function(N, d, L, alpha, s,
                                             coupling = c("independent", "shared_clock"),
                                             init = c("clustered", "iid_uniform"),
                                             times = c(0,1,2,5,10,20,50,100,200,500,1000)) {
  coupling <- match.arg(coupling)
  init <- match.arg(init)
  
  times <- sort(unique(times))
  tmax <- max(times)
  pos <- init_positions(N, d, L, init = init)
  
  out <- tibble(n = times, p_mc = NA_real_)
  if (0 %in% times) out$p_mc[out$n == 0] <- mean(collision_indicator(pos))
  
  if (tmax >= 1) {
    for (t in 1:tmax) {
      pos <- step_batch(pos, d, alpha, s, coupling = coupling)
      if (t %in% times) out$p_mc[out$n == t] <- mean(collision_indicator(pos))
    }
  }
  out
}

simulate_Tcoll <- function(N, d, L, alpha, s,
                           coupling = c("independent", "shared_clock"),
                           init = c("clustered", "iid_uniform"),
                           n_max = 200000L) {
  coupling <- match.arg(coupling)
  init <- match.arg(init)
  
  pos <- init_positions(N, d, L, init = init)
  T <- rep(NA_integer_, N)
  
  hit0 <- collision_indicator(pos)
  T[hit0] <- 0L
  
  for (n in 1:n_max) {
    alive <- which(is.na(T))
    if (length(alive) == 0) break
    
    sub <- pos[alive, , drop = FALSE]
    sub <- step_batch(sub, d, alpha, s, coupling = coupling)
    pos[alive, ] <- sub
    
    hit <- collision_indicator(sub)
    if (any(hit)) T[alive[hit]] <- n
  }
  
  T[is.na(T)] <- n_max
  T
}

############################################################
## PLOT 1 — P(R_n in C) vs n (MC dots, theory solid, bounds dashed)
############################################################
d0 <- 50
L0 <- 8
s0 <- s_fixed

alpha_list <- tibble(
  alpha_label = c("0 (no teleport)", "1/(10 d^2) low", "1/d^2 cross", "0.05 high", "0.20 very high"),
  alpha = c(0, 1/(10*d0^2), 1/(d0^2), 0.05, 0.20)
)

time_grid <- c(0,1,2,5,10,20,50,100,200,500,1000,1500,2000,2500,3000)
N_prob <- 8000

grid1 <- tidyr::crossing(alpha_list, coupling = c("independent", "shared_clock"))

prob_panels <- purrr::pmap_dfr(
  grid1,
  function(alpha_label, alpha, coupling) {
    sim <- simulate_collision_prob_timegrid(
      N = N_prob, d = d0, L = L0, alpha = alpha, s = s0,
      coupling = coupling,
      init = "clustered",
      times = time_grid
    )
    piC <- pi_collision(d0, L0)
    
    tvb <- if (coupling == "independent") {
      tv_bound_independent(L0, alpha, sim$n)
    } else {
      tv_bound_shared_clock(alpha, sim$n)
    }
    lb <- lower_bound_collision_prob(piC, tvb)
    ub <- vapply(sim$n, function(nn) union_bound_collision_time_n(d0, L0, alpha, s0, nn), numeric(1))
    
    tibble(
      alpha_label = alpha_label,
      alpha = alpha,
      coupling = coupling,
      n = sim$n,
      p_mc = sim$p_mc,
      piC = piC,
      lb = lb,
      ub = ub
    )
  }
) %>%
  mutate(
    coupling = recode(coupling,
                      independent = "Independent teleports",
                      shared_clock = "Shared teleport clock")
  )

prob_long <- prob_panels %>%
  pivot_longer(c(p_mc, piC, lb, ub), names_to = "curve", values_to = "value") %>%
  mutate(
    curve = recode(curve,
                   p_mc = "Monte Carlo",
                   piC  = "Stationary  pi(C)",
                   lb   = "Lower bound  pi(C) - TV",
                   ub   = "Union upper bound"),
    curve_type = case_when(
      curve == "Monte Carlo" ~ "Monte Carlo",
      grepl("bound", curve, ignore.case = TRUE) ~ "Bound",
      TRUE ~ "Theory"
    ),
    curve = factor(curve, levels = c("Monte Carlo","Stationary  pi(C)","Lower bound  pi(C) - TV","Union upper bound"))
  )

p1 <- ggplot(prob_long, aes(x = n, y = value, color = curve)) +
  geom_line(
    data = subset(prob_long, curve_type != "Monte Carlo"),
    aes(linetype = curve_type),
    linewidth = 0.9, alpha = 0.80
  ) +
  geom_point(
    data = subset(prob_long, curve_type == "Monte Carlo"),
    size = 1.1, alpha = 0.55
  ) +
  facet_grid(coupling ~ alpha_label) +
  scale_linetype_manual(values = c("Theory" = "solid", "Bound" = "dashed")) +
  scale_color_manual(values = ctc_named_colors(levels(prob_long$curve))) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "time n", y = "probability", color = "", linetype = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        strip.text.x = element_text(size = 9),
        strip.text.y = element_text(size = 10))
p1 <- bump_fonts(p1, base = 16)
print(p1)
save_png_hq(p1, "teleport_collision_prob_p1.png", width = 16, height = 9, dpi = 600)

############################################################
## PLOT 2 — Tail of T_coll (MC dots; theory solid; bounds dashed)
############################################################
N_tail <- 30000
alpha_tail <- tibble(
  alpha_label = c("1/(10 d^2) low", "1/d^2 cross", "0.05 high"),
  alpha = c(1/(10*d0^2), 1/(d0^2), 0.05)
)

n_max_tail <- 20000L

tail_df <- purrr::pmap_dfr(
  alpha_tail,
  function(alpha_label, alpha) {
    T <- simulate_Tcoll(
      N = N_tail, d = d0, L = L0, alpha = alpha, s = s0,
      coupling = "independent", init = "clustered",
      n_max = n_max_tail
    )
    
    piC <- pi_collision(d0, L0)
    n_grid <- 0:max(T)
    
    geom_tail <- (1 - piC)^n_grid
    p_coll <- (alpha^L0) * piC
    one_step_bound <- exp(-p_coll * n_grid)
    
    rhs <- 1 - (1 - piC/2)^(1/L0)
    t_block <- ceiling(log(rhs) / log(1 - alpha))
    t_block <- max(t_block, 1L)
    block_bound <- exp(-floor(n_grid / t_block) * (piC/2))
    
    # MC survival on grid (as points)
    S_emp <- sapply(n_grid, function(n) mean(T > n))
    
    tibble(
      alpha_label = alpha_label,
      n = n_grid,
      S_emp = S_emp,
      S_geom = geom_tail,
      S_one = one_step_bound,
      S_block = block_bound
    )
  }
)

tail_long <- tail_df %>%
  pivot_longer(-c(alpha_label, n), names_to = "curve", values_to = "S") %>%
  mutate(
    curve = recode(curve,
                   S_emp  = "Monte Carlo",
                   S_geom = "Geom approx  (1-pi(C))^n",
                   S_one  = "One-step bound  exp(-alpha^L pi(C) n)",
                   S_block= "Block bound  exp(-floor(n/t)*pi(C)/2)"),
    curve_type = case_when(
      curve == "Monte Carlo" ~ "Monte Carlo",
      grepl("bound", curve, ignore.case = TRUE) ~ "Bound",
      TRUE ~ "Theory"
    ),
    curve = factor(curve, levels = c("Monte Carlo",
                                     "Geom approx  (1-pi(C))^n",
                                     "One-step bound  exp(-alpha^L pi(C) n)",
                                     "Block bound  exp(-floor(n/t)*pi(C)/2)"))
  )

p2 <- ggplot(tail_long, aes(x = n, y = S, color = curve)) +
  geom_line(
    data = subset(tail_long, curve_type != "Monte Carlo"),
    aes(linetype = curve_type),
    linewidth = 0.9, alpha = 0.78
  ) +
  geom_point(
    data = subset(tail_long, curve_type == "Monte Carlo"),
    size = 0.8, alpha = 0.55
  ) +
  facet_wrap(~ alpha_label, ncol = 3) +
  scale_linetype_manual(values = c("Theory" = "solid", "Bound" = "dashed")) +
  scale_color_manual(values = ctc_named_colors(levels(tail_long$curve))) +
  scale_y_log10(labels = percent_format(accuracy = 0.1)) +
  labs(x = "n", y = "P(T_coll > n)", color = "", linetype = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
p2 <- bump_fonts(p2, base = 16)
print(p2)
save_png_hq(p2, "teleport_collision_tail_p2.png", width = 14, height = 6.5, dpi = 600)

############################################################
## PLOT 3 — Mean(T_coll) across regimes (MC dots; theory solid)
############################################################
d_vec <- c(50, 100, 200)
L_vec <- c(4, 8, 16)
N_mean <- 8000

mean_grid <- tidyr::crossing(d = d_vec, L = L_vec) %>%
  filter(L <= d) %>%
  mutate(
    alpha_low   = 1/(10*d^2),
    alpha_cross = 1/(d^2),
    alpha_high  = 0.05
  ) %>%
  pivot_longer(c(alpha_low, alpha_cross, alpha_high),
               names_to = "regime", values_to = "alpha") %>%
  mutate(
    regime = recode(regime,
                    alpha_low = "low  1/(10 d^2)",
                    alpha_cross = "cross  1/d^2",
                    alpha_high = "high  0.05")
  )

mean_stats <- mean_grid %>%
  group_by(d, L, regime, alpha) %>%
  group_modify(~{
    d0 <- .y$d
    L0 <- .y$L
    a0 <- .y$alpha
    
    T <- simulate_Tcoll(
      N = N_mean, d = d0, L = L0, alpha = a0, s = s_fixed,
      coupling = "independent", init = "clustered",
      n_max = 200000L
    )
    
    piC <- pi_collision(d0, L0)
    
    tibble(
      mc_mean        = mean(T),
      approx_occ     = 1 / piC,
      approx_sparse  = 2 * d0 / (L0 * (L0 - 1))
    )
  }) %>%
  ungroup()


mean_long <- mean_stats %>%
  pivot_longer(c(mc_mean, approx_occ, approx_sparse),
               names_to = "curve", values_to = "value") %>%
  mutate(
    curve = recode(curve,
                   mc_mean = "Monte Carlo mean",
                   approx_occ = "Occ. approx  1/pi(C)",
                   approx_sparse = "Sparse approx  2d/(L(L-1))"),
    curve = factor(curve, levels = c("Monte Carlo mean","Occ. approx  1/pi(C)","Sparse approx  2d/(L(L-1))"))
  )

p3 <- ggplot(mean_long, aes(x = alpha, y = value, color = curve)) +
  geom_line(
    data = subset(mean_long, curve != "Monte Carlo mean"),
    linewidth = 0.9, alpha = 0.78, linetype = "solid"
  ) +
  geom_point(
    data = subset(mean_long, curve == "Monte Carlo mean"),
    size = 1.8, alpha = 0.55
  ) +
  facet_grid(d ~ L, labeller = labeller(d = function(x) paste0("d=", x),
                                        L = function(x) paste0("L=", x))) +
  scale_x_log10(labels = scientific_format(digits = 2)) +
  scale_y_log10(labels = comma_format()) +
  scale_color_manual(values = ctc_named_colors(levels(mean_long$curve))) +
  labs(x = "teleport rate alpha (log scale)", y = "mean(T_coll) (log scale)", color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
p3 <- bump_fonts(p3, base = 16)
print(p3)
save_png_hq(p3, "teleport_collision_mean_p3.png", width = 14, height = 8, dpi = 600)
################################################################################
################################################################################
##  SECTION 5 — Interaction Rules Validation: Blocking / Bouncing / Pushing
##  (MODIFIED TO BE COMPUTABLE: hard caps + adaptive reps)
################################################################################
################################################################################

## Packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

## --- Helpers (safe to keep even if you defined them earlier) -----------------
mod_d <- function(x, d) ((x %% d) + d) %% d

# If you already have these in earlier sections, you can delete these defaults.
if (!exists("OI")) {
  OI <- c(grey = "grey50")
}
if (!exists("ctc_rule_colors")) {
  ctc_rule_colors <- c(
    "Baseline (none)"   = "#1b9e77",
    "Blocking (sync)"   = "#d95f02",
    "Bouncing (sync)"   = "#7570b3",
    "Pushing (seq)"     = "#e7298a"
  )
}
if (!exists("bump_fonts")) {
  bump_fonts <- function(p, base = 14) p + theme(text = element_text(size = base))
}
if (!exists("save_png_hq")) {
  save_png_hq <- function(p, filename, width = 12, height = 7, dpi = 300) {
    ggsave(filename, plot = p, width = width, height = height, dpi = dpi, units = "in")
  }
}

set.seed(123)

pi_collision <- function(d, L) {
  if (L > d) return(1)
  falling <- prod(d:(d - L + 1))
  1 - falling / (d^L)
}

make_mu_forward_unit <- function() c(`1` = 1.0)
mu <- make_mu_forward_unit()
s  <- 0.30

sample_mu <- function(n, mu) {
  inc  <- as.integer(names(mu))
  prob <- as.numeric(mu)
  sample(inc, size = n, replace = TRUE, prob = prob)
}

propose_independent <- function(x, d, alpha, s, mu) {
  L <- length(x)
  
  tele <- runif(L) < alpha
  y <- x
  move_type <- rep("stay", L)
  inc_used  <- integer(L)
  
  if (any(tele)) {
    y[tele] <- sample.int(d, sum(tele), replace = TRUE) - 1L
    move_type[tele] <- "teleport"
  }
  
  not_tele <- !tele
  if (any(not_tele)) {
    lazy <- (runif(L) < s) & not_tele
    inc_idx <- not_tele & !lazy
    move_type[lazy] <- "lazy"
    
    if (any(inc_idx)) {
      X <- sample_mu(sum(inc_idx), mu)
      inc_used[inc_idx] <- X
      y[inc_idx] <- mod_d(x[inc_idx] + X, d)
      move_type[inc_idx] <- "increment"
    }
  }
  
  list(y = y, move_type = move_type, inc = inc_used, tele = tele)
}

resolve_interaction <- function(x, prop, d,
                                interaction = c("none","blocking","bouncing","pushing"),
                                timing = c("synchronous","sequential"),
                                teleport_constrained = FALSE) {
  interaction <- match.arg(interaction)
  timing <- match.arg(timing)
  
  y <- prop$y
  move_type <- prop$move_type
  inc <- prop$inc
  tele <- prop$tele
  L <- length(x)
  
  if (interaction == "none") {
    return(list(x_new = y, modified = any(y != x)))
  }
  
  if (teleport_constrained && any(tele)) {
    occ <- unique(x)
    vacant <- setdiff(0:(d-1), occ)
    ntele <- sum(tele)
    if (length(vacant) >= ntele) {
      y[tele] <- sample(vacant, ntele, replace = FALSE)
    } else if (length(vacant) > 0) {
      y[tele] <- sample(vacant, ntele, replace = TRUE)
    } else {
      y[tele] <- sample.int(d, ntele, replace = TRUE) - 1L
    }
  }
  
  occ_other_vec <- function(i) x[-i]
  
  if (timing == "synchronous") {
    x_new <- y
    nontele_idx <- which(!tele)
    
    if (length(nontele_idx) > 0) {
      targets <- y[nontele_idx]
      dup_targets <- targets[duplicated(targets)]
      dup_set <- unique(dup_targets)
      
      for (i in nontele_idx) {
        tgt <- y[i]
        blocked_by_occ <- (tgt != x[i]) && (tgt %in% occ_other_vec(i))
        blocked_by_dup <- (tgt %in% dup_set)
        
        if (interaction == "blocking") {
          if (blocked_by_occ || blocked_by_dup) x_new[i] <- x[i]
        }
        
        if (interaction == "bouncing") {
          if ((blocked_by_occ || blocked_by_dup) && move_type[i] == "increment") {
            refl <- mod_d(x[i] - inc[i], d)
            refl_blocked_occ <- (refl %in% occ_other_vec(i))
            refl_blocked_dup <- (refl %in% dup_set)
            x_new[i] <- if (!refl_blocked_occ && !refl_blocked_dup) refl else x[i]
          } else if (blocked_by_occ || blocked_by_dup) {
            x_new[i] <- x[i]
          }
        }
        
        if (interaction == "pushing") {
          return(resolve_interaction(x, prop, d,
                                     interaction = "pushing",
                                     timing = "sequential",
                                     teleport_constrained = teleport_constrained))
        }
      }
    }
    
    modified <- any(x_new != y)
    return(list(x_new = x_new, modified = modified))
  }
  
  cur <- x
  modified <- FALSE
  
  for (i in seq_len(L)) {
    if (tele[i]) {
      cur[i] <- y[i]
      next
    }
    
    tgt <- y[i]
    if (tgt == cur[i]) next
    
    occupied <- which(cur == tgt)
    
    if (length(occupied) == 0) {
      cur[i] <- tgt
      next
    }
    
    if (interaction == "blocking") {
      modified <- TRUE
      next
    }
    
    if (interaction == "bouncing" && move_type[i] == "increment") {
      refl <- mod_d(cur[i] - inc[i], d)
      if (!any(cur == refl)) {
        cur[i] <- refl
      } else {
        modified <- TRUE
      }
      next
    }
    
    if (interaction == "pushing") {
      victim <- occupied[1]
      push_to <- mod_d(tgt + 1L, d)
      if (!any(cur == push_to)) {
        cur[victim] <- push_to
        cur[i] <- tgt
      } else {
        modified <- TRUE
      }
      next
    }
    
    modified <- TRUE
  }
  
  list(x_new = cur, modified = modified)
}

simulate_one_path <- function(d, L, alpha, s, mu,
                              interaction = c("none","blocking","bouncing","pushing"),
                              timing = c("synchronous","sequential"),
                              init = c("iid_uniform","clustered"),
                              teleport_constrained = FALSE,
                              n_max = 1e6) {
  interaction <- match.arg(interaction)
  timing <- match.arg(timing)
  init <- match.arg(init)
  
  x <- if (init == "iid_uniform") sample.int(d, L, replace = TRUE) - 1L else mod_d(0:(L-1), d)
  
  if (any(duplicated(x))) {
    return(list(Tcoll = 0L, censored = FALSE, mod_rate = 0))
  }
  
  mod_count <- 0L
  
  for (n in seq_len(n_max)) {
    prop <- propose_independent(x, d, alpha, s, mu)
    res  <- resolve_interaction(x, prop, d,
                                interaction = interaction,
                                timing = timing,
                                teleport_constrained = teleport_constrained)
    x <- res$x_new
    if (isTRUE(res$modified)) mod_count <- mod_count + 1L
    
    if (any(duplicated(x))) {
      return(list(Tcoll = n, censored = FALSE, mod_rate = mod_count / n))
    }
  }
  
  list(Tcoll = as.integer(n_max), censored = TRUE, mod_rate = mod_count / n_max)
}

simulate_many <- function(N_rep, d, L, alpha, s, mu,
                          interaction = "none",
                          timing = "synchronous",
                          init = "clustered",
                          teleport_constrained = FALSE,
                          n_max = 1e6,
                          seed = 1) {
  set.seed(seed)
  
  out <- vector("list", N_rep)
  for (m in seq_len(N_rep)) {
    out[[m]] <- simulate_one_path(d, L, alpha, s, mu,
                                  interaction = interaction,
                                  timing = timing,
                                  init = init,
                                  teleport_constrained = teleport_constrained,
                                  n_max = n_max)
  }
  
  tibble(
    d = d, L = L, alpha = alpha, s = s,
    interaction = interaction, timing = timing,
    init = init, teleport_constrained = teleport_constrained,
    Tcoll = vapply(out, `[[`, integer(1), "Tcoll"),
    censored = vapply(out, `[[`, logical(1), "censored"),
    mod_rate = vapply(out, `[[`, numeric(1), "mod_rate")
  )
}

## --- GRID -------------------------------------------------------------------
d_grid <- c(50, 100, 200)
L_grid <- c(4, 8, 16)

alpha_regimes <- function(d) {
  c(low = 1/(10*d^2), cross = 1/(d^2), high = 0.05)
}

rules <- tibble(
  rule_id = c("none_sync", "block_sync", "bounce_sync", "push_seq"),
  interaction = c("none", "blocking", "bouncing", "pushing"),
  timing = c("synchronous", "synchronous", "synchronous", "sequential"),
  label = c("Baseline (none)", "Blocking (sync)", "Bouncing (sync)", "Pushing (seq)")
)

## --- COMPUTABLE SETTINGS (KEY CHANGE) ----------------------------------------
N_rep_map <- function(alpha_regime, interaction) {
  if (alpha_regime == "high") return(2500)
  if (alpha_regime == "cross") return(1200)
  # low
  if (interaction == "none") return(800)
  return(250)
}

n_max_map <- function(alpha_regime, interaction, d, L, alpha) {
  if (alpha_regime == "high")  return(2e4)
  if (alpha_regime == "cross") return(8e4)
  # low
  if (interaction == "none") return(1e5)
  return(2e5)
}

## --- RUN --------------------------------------------------------------------
all_runs <- list()
idx <- 1

for (d in d_grid) {
  for (L in L_grid) {
    a_vec <- alpha_regimes(d)
    
    for (a_name in names(a_vec)) {
      alpha <- unname(a_vec[a_name])
      
      for (rr in seq_len(nrow(rules))) {
        interaction <- rules$interaction[rr]
        timing <- rules$timing[rr]
        
        N_rep <- N_rep_map(a_name, interaction)
        n_max <- n_max_map(a_name, interaction, d, L, alpha)
        
        message(sprintf(
          "[%03d] d=%d L=%d alpha_reg=%s alpha=%.3g rule=%s timing=%s | N_rep=%d n_max=%d",
          idx, d, L, a_name, alpha, interaction, timing, N_rep, n_max
        ))
        
        sim <- simulate_many(
          N_rep = N_rep,
          d = d, L = L, alpha = alpha, s = s, mu = mu,
          interaction = interaction,
          timing = timing,
          init = "clustered",
          teleport_constrained = FALSE,
          n_max = n_max,
          seed = 1000 + 37*idx
        )
        
        sim$alpha_regime <- a_name
        sim$rule_label <- rules$label[rr]
        all_runs[[idx]] <- sim
        idx <- idx + 1
      }
    }
  }
}

df <- bind_rows(all_runs) %>%
  mutate(
    alpha_regime = factor(alpha_regime, levels = c("low","cross","high")),
    rule_label   = factor(rule_label, levels = rules$label),
    d = factor(d), L = factor(L)
  )

summ <- df %>%
  group_by(d, L, alpha_regime, alpha, rule_label) %>%
  summarise(
    mean_T = mean(Tcoll),
    se_T   = sd(Tcoll)/sqrt(n()),
    cens_rate = mean(censored),
    mean_mod = mean(mod_rate),
    .groups = "drop"
  ) %>%
  mutate(
    piC = mapply(function(dd, LL) pi_collision(as.integer(as.character(dd)),
                                               as.integer(as.character(LL))), d, L),
    geom_mean_ref = 1 / piC
  )

## reference per facet (d, L)
ref_panel <- summ %>% distinct(d, L, geom_mean_ref)

pA <- ggplot(summ, aes(x = alpha_regime, y = mean_T, color = rule_label, group = rule_label)) +
  geom_line(alpha = 0.65, linewidth = 0.7) +
  geom_point(alpha = 0.65, size = 1.8) +
  geom_errorbar(aes(ymin = pmax(mean_T - 2*se_T, 1), ymax = mean_T + 2*se_T),
                alpha = 0.25, width = 0.12) +
  geom_hline(
    data = ref_panel,
    aes(yintercept = geom_mean_ref),
    inherit.aes = FALSE,
    linetype = "dashed", alpha = 0.60, linewidth = 0.9,
    color = OI["grey"]
  ) +
  facet_grid(d ~ L, scales = "free_y") +
  scale_color_manual(values = ctc_rule_colors) +
  scale_y_log10(labels = comma_format()) +
  labs(x = "teleport regime", y = "mean T_coll (log scale)", color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

pA <- bump_fonts(pA, base = 16)
print(pA)
save_png_hq(pA, "interaction_rules_pA_mean_Tcoll.png", width = 16, height = 9, dpi = 600)

## Panel B: modification rate (dots)
pB <- ggplot(summ, aes(x = alpha_regime, y = mean_mod, color = rule_label, group = rule_label)) +
  geom_line(alpha = 0.65, linewidth = 0.7) +
  geom_point(alpha = 0.65, size = 1.8) +
  facet_grid(d ~ L, scales = "free_y") +
  scale_color_manual(values = ctc_rule_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "teleport regime", y = "mean modification rate", color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
pB <- bump_fonts(pB, base = 16)
print(pB)
save_png_hq(pB, "interaction_rules_pB_mod_rate.png", width = 16, height = 9, dpi = 600)

## Tail comparison
emp_survival <- function(t) {
  u <- sort(unique(t))
  S <- vapply(u, function(x) mean(t > x), numeric(1))
  tibble(t = u, S = S)
}

pick_d <- "100"
pick_L <- "8"
pick_reg <- "high"

sub <- df %>% filter(d == pick_d, L == pick_L, alpha_regime == pick_reg)

tail_df <- sub %>%
  group_by(rule_label) %>%
  group_modify(~ emp_survival(.x$Tcoll)) %>%
  ungroup() %>%
  mutate(rule_label = factor(rule_label, levels = rules$label))

d_num <- as.integer(as.character(pick_d))
L_num <- as.integer(as.character(pick_L))
piC   <- pi_collision(d_num, L_num)

ref_df <- tibble(
  t = sort(unique(tail_df$t)),
  S = (1 - piC)^sort(unique(tail_df$t))
)

pC <- ggplot(tail_df, aes(x = t, y = S, color = rule_label)) +
  geom_point(alpha = 0.55, size = 1.0) +
  geom_line(data = ref_df, aes(x = t, y = S),
            inherit.aes = FALSE,
            linetype = "dashed", alpha = 0.70, linewidth = 0.9, color = OI["grey"]) +
  scale_color_manual(values = ctc_rule_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "t", y = "P(T_coll > t)", color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
pC <- bump_fonts(pC, base = 16)
print(pC)
save_png_hq(pC, "interaction_rules_pC_tail.png", width = 14, height = 8, dpi = 600)
