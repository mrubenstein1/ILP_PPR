# ══════════════════════════════════════════════════════════════════════════════
# 07__ilp_core.R  —  Shared engine for the three land-acquisition decision models
# ══════════════════════════════════════════════════════════════════════════════
#
# This script defines (but does not run) the machinery shared by all three models:
#   * build_scenario_matrices()  : assemble per-EAU (benefit, hazard, cost) matrices
#   * compute_survival_matrix()  : cumulative survival S[i,t] from per-period hazard
#   * compute_value_vector()     : the ILP objective coefficient V[i,tau]
#   * make_budget()              : per-period budget in EAU-cost units
#   * solve_acquisition_ilp()    : Gurobi wrapper for the binary program
#   * run_greedy/myopic/rolling  : the three policies
#   * run_full_horizon()         : single full-horizon solve (validation reference)
#   * evaluate_policy()          : TRUE-future landscape ducks for a given schedule
#   * make_synthetic_instance()  : small synthetic generator (used by validation)
#
# It is sourced by 08__run_models.R (to run the experiment) and 09__validation.R
# (to check correctness). Source it, then call the functions.
#
# ──────────────────────────────────────────────────────────────────────────────
# OBJECTIVE (landscape-total).  We maximise the time-discounted sum of breeding
# duck pairs on the WHOLE landscape. A parcel acquired at period tau contributes
# its (climate-driven) abundance b[i,t] with certainty for every t >= tau; a
# parcel left in the unprotected pool contributes b[i,t] * S[i,t] in expectation,
# where S is the probability it has not yet been lost to conversion. Total welfare
# decomposes as
#       J(policy) = J_baseline + sum_i V[i, tau_i] ,
# where J_baseline (no acquisition) is a constant and
#       V[i,tau] = sum_{t>=tau} delta^t * b[i,t] * (1 - S[i,t])
# is the incremental welfare of acquiring parcel i at tau (= the expected
# conversion loss it prevents). Maximising sum_i V[i,tau_i] therefore maximises J.
#
# CONVERSION RISK enters ONLY through the (1 - S) weighting in compute_value_vector()
# and the S weighting in evaluate_policy(). Both read the hazard column trans_prob.
# If the conversion metric is later revised (new dataset, different definition),
# only the trans_prob column upstream changes — the model logic here is untouched.
#
# SURVIVAL CLOCK.  S is accumulated GLOBALLY from t = 0 (period 2020) under each
# model's belief about the hazard trajectory. This is the convention that makes the
# rolling-horizon policy the provable welfare optimum and makes the stationary null
# exact (under a stationary scenario the myopic frozen-future belief is correct, so
# myopic and rolling become mathematically identical). It was selected after testing
# the alternative (re-planning survival forward from each decision period), which
# breaks rolling-optimality (myopic can then beat rolling — an order violation).
#
# DISCOUNTING is absolute present value to 2020: delta^t with t the decade index
# (0 = 2020, ..., 8 = 2100). delta is a free parameter for sensitivity analysis.
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(Matrix)   # sparse constraint matrix for Gurobi
})

# ══ PARAMETERS (edit here for sensitivity analyses) ══════════════════════════════####

# Per-decade discount factor. 2100 is worth delta^8 of present value
# (0.95^8 ~= 0.66). Change to test sensitivity of results to discounting.
DELTA <- 0.95

# Budget per period, expressed in "median EAUs": the per-period budget is
# BUDGET_EAUS_PER_PERIOD * median(cost at that period). This scales automatically
# with inflation (cost already carries 2%/yr) and is dimensionless. Default 5;
# sensitivity range 2-10.
BUDGET_EAUS_PER_PERIOD <- 5

# Decision periods: 2020 (baseline) stitched onto the decadal trajectory 2030-2100.
DECISION_YEARS <- seq(2020, 2100, by = 10)   # length 9; index t = (year-2020)/10

# Scenario coding in data_panel (rcp / gcm factor levels). The 2020 baseline row
# (rcp = "baseline", gcm = "baseline") is used as period 0 for EVERY scenario.
# NOTE: verify these strings match the levels in your data_panel.
SCENARIOS <- list(
  rcp45_wet  = c(rcp = "45",         gcm = "wet"),
  rcp45_dry  = c(rcp = "45",         gcm = "dry"),
  rcp85_wet  = c(rcp = "85",         gcm = "wet"),
  rcp85_dry  = c(rcp = "85",         gcm = "dry"),
  stationary = c(rcp = "stationary", gcm = "stationary")
)
# ══════════════════════════════════════════════════════════════════════════════


# ── 1. Build per-EAU (benefit, hazard, cost) matrices for one scenario ────────####
#
# Returns matrices of shape [n_eau x n_t] (rows = EAUs in sorted id order, columns
# = decision periods 2020..2100). Period 0 (2020) is taken from the shared baseline
# row; periods 1..8 (2030..2100) from the scenario's own rows. Cost is
# scenario-invariant, so it is pulled from the same rows for convenience.
#
#   b   = scaled_abundance  (EAU-level breeding pairs)
#   lam = trans_prob        (per-period conversion probability proxy)
#   cost = cost             (USD to acquire EAU at that period)

build_scenario_matrices <- function(data_panel, scenario,
                                     years = DECISION_YEARS) {
  rcp_s <- scenario[["rcp"]]
  gcm_s <- scenario[["gcm"]]
  decade_years <- years[years > 2020]

  eau_ids <- sort(unique(data_panel$eau_id))
  n_eau   <- length(eau_ids)
  n_t     <- length(years)

  # period-0 (2020) baseline slice
  base <- data_panel %>%
    filter(year == 2020, rcp == "baseline", gcm == "baseline") %>%
    select(eau_id, prop_suitable, scaled_abundance, trans_prob, cost)

  # periods 1..8 scenario slice
  scen <- data_panel %>%
    filter(year %in% decade_years, rcp == rcp_s, gcm == gcm_s) %>%
    select(eau_id, year, scaled_abundance, trans_prob, cost)

  # ---- integrity checks: one row per (eau) baseline, one per (eau, year) scenario
  checks <- c(
    "baseline has one row per EAU"            = nrow(base) == n_eau &&
                                                 !any(duplicated(base$eau_id)),
    "scenario has one row per EAU-year"       = nrow(scen) == n_eau * length(decade_years) &&
                                                 !any(duplicated(scen[c("eau_id", "year")])),
    "all EAUs present in scenario slice"      = setequal(unique(scen$eau_id), eau_ids),
    "no missing benefit/hazard/cost"          = !anyNA(base[c("scaled_abundance","trans_prob","cost")]) &&
                                                 !anyNA(scen[c("scaled_abundance","trans_prob","cost")])
  )
  if (!all(checks)) {
    cat("build_scenario_matrices FAILED for rcp =", rcp_s, "gcm =", gcm_s, "\n")
    for (nm in names(checks)) cat(sprintf("  [%s] %s\n", ifelse(checks[[nm]], "PASS", "FAIL"), nm))
    stop("Scenario assembly integrity check failed.")
  }

  # allocate
  b    <- matrix(NA_real_, n_eau, n_t, dimnames = list(eau_ids, years))
  lam  <- matrix(NA_real_, n_eau, n_t, dimnames = list(eau_ids, years))
  cost <- matrix(NA_real_, n_eau, n_t, dimnames = list(eau_ids, years))

  # period 0 (2020) from baseline, aligned by eau order
  bi <- match(eau_ids, base$eau_id)
  b[,    1] <- base$scaled_abundance[bi]
  lam[,  1] <- base$trans_prob[bi]
  cost[, 1] <- base$cost[bi]

  # periods 1..8 from scenario, aligned by (eau, year)
  for (k in seq_along(decade_years)) {
    yr   <- decade_years[k]
    slab <- scen %>% filter(year == yr)
    si   <- match(eau_ids, slab$eau_id)
    b[,    k + 1] <- slab$scaled_abundance[si]
    lam[,  k + 1] <- slab$trans_prob[si]
    cost[, k + 1] <- slab$cost[si]
  }

  list(b = b, lam = lam, cost = cost, eau_ids = eau_ids, years = years)
}


# ── 2. Core quantities: survival, acquisition value, budget ──────────────────####

# Cumulative GLOBAL survival from t = 0. S[,1] = 1 (no loss possible in the first
# period); S[,t] = S[,t-1] * (1 - lam[,t-1]). The terminal-period hazard never
# enters S (there is no period after it), matching the upstream convention that
# sets trans_prob = 0 at 2100.
compute_survival_matrix <- function(lam_mat) {
  n_t <- ncol(lam_mat)
  S <- matrix(1, nrow(lam_mat), n_t, dimnames = dimnames(lam_mat))
  if (n_t >= 2) for (t in 2:n_t) S[, t] <- S[, t - 1] * (1 - lam_mat[, t - 1])
  S
}

# Acquisition-value vector for ONE parcel given its benefit and hazard trajectories.
# Returns V[tau] for tau = period 0..(n_t-1) (R index 1..n_t):
#     V[tau] = sum_{t >= tau} delta^t * b[t] * (1 - S[t])
# computed with a reverse cumulative sum. Survival is global from t = 0.
compute_value_vector <- function(b_vec, lam_vec, delta = DELTA) {
  n_t <- length(b_vec)
  S <- numeric(n_t); S[1] <- 1
  if (n_t >= 2) for (t in 2:n_t) S[t] <- S[t - 1] * (1 - lam_vec[t - 1])
  disc    <- delta^(0:(n_t - 1))
  contrib <- disc * b_vec * (1 - S)        # per-period loss-prevented value
  rev(cumsum(rev(contrib)))                # V[tau] = sum_{t>=tau} contrib[t]
}

# Per-period budget vector (length n_t) in cost units: n_eaus * median period cost.
make_budget <- function(cost_mat, n_eaus = BUDGET_EAUS_PER_PERIOD) {
  apply(cost_mat, 2, function(col) n_eaus * median(col))
}


# ── 3. ILP solver (Gurobi) ───────────────────────────────────────────────────####
#
#   maximise   sum_{i in avail, t in periods}  V[i,t] * y[i,t]
#   subject to sum_t y[i,t] <= 1                 (each parcel acquired at most once)
#              sum_i cost[i,t] * y[i,t] <= B[t]  (per-period budget, no rollover)
#              y in {0,1}      (or [0,1] if relax = TRUE, for the LP bound)
#
# Returns list(picks = data.frame(eau_idx, period), objval, x, grid). `eau_idx` and
# `period` are 1-based indices into the rows/cols of the supplied matrices.
#
# Requires the gurobi R package (ships with a Gurobi install; activate a free
# academic licence first). To swap in another MILP solver, replace ONLY the block
# marked "SOLVER CALL" — the model is assembled in the standard (A, obj, sense,
# rhs, vtype) form that most R MILP interfaces accept.

solve_acquisition_ilp <- function(V_mat, cost_mat, budget, avail, periods,
                                   relax = FALSE) {
  if (length(avail) == 0 || length(periods) == 0)
    return(list(picks = data.frame(eau_idx = integer(0), period = integer(0)),
                objval = 0, x = numeric(0), grid = NULL))

  # one binary variable per (available EAU, allowed period)
  grid  <- expand.grid(eau_idx = avail, period = periods)
  nvar  <- nrow(grid)
  obj   <- V_mat[cbind(grid$eau_idx, grid$period)]
  costs <- cost_mat[cbind(grid$eau_idx, grid$period)]

  n_av <- length(avail); n_pe <- length(periods)
  # constraint rows: 1..n_av = acquire-once (per EAU); next n_pe = budget (per period)
  row_once <- match(grid$eau_idx, avail)
  row_budg <- n_av + match(grid$period, periods)


  # ---- SOLVER CALL (Gurobi) -------------------------------------------------
  if (!requireNamespace("gurobi", quietly = TRUE))
    stop("The 'gurobi' R package is required. Install it from your Gurobi ",
         "installation (install.packages('<gurobi>/R/gurobi_x.y-z.tar.gz', ",
         "repos = NULL)) and activate an academic licence.")
  
  # Scale costs to avoid numerical issues from large USD values (billions range)
  cost_scale <- 1e6
  costs_s    <- costs / cost_scale
  
  A <- Matrix::sparseMatrix(
    i    = c(row_once, row_budg),
    j    = c(seq_len(nvar), seq_len(nvar)),
    x    = c(rep(1, nvar), costs_s),
    dims = c(n_av + n_pe, nvar)
  )
  sense <- rep("<=", n_av + n_pe)
  rhs   <- c(rep(1, n_av), budget[periods] / cost_scale)
  
  model <- list(
    A          = A,
    obj        = obj,
    sense      = sense,
    rhs        = rhs,
    lb         = rep(0, nvar),
    ub         = rep(1, nvar),
    vtype      = if (relax) "C" else "B",
    modelsense = "max"
  )
  res <- gurobi::gurobi(model, params = list(
    OutputFlag   = 1
  ))
  # ---------------------------------------------------------------------------

  x <- res$x
  chosen <- if (relax) integer(0) else which(x > 0.5)
  picks  <- grid[chosen, , drop = FALSE]
  list(picks = picks, objval = res$objval, x = x, grid = grid)
}


# ── 4. The three policies (+ full-horizon reference) ─────────────────────────####
#
# Each returns an integer vector `acquired` of length n_eau giving the period index
# (1..n_t) at which each EAU is acquired, or NA if never acquired.

# ROLLING HORIZON — re-solve each period using the TRUE projected future (global V),
# implement only the current period, advance. Under deterministic foresight this
# reproduces the full-horizon optimum (see run_full_horizon); the loop form is kept
# for methodological fidelity and so stochastic extensions can slot in.
run_rolling_horizon <- function(b_mat, lam_mat, cost_mat, budget, delta = DELTA) {
  n_eau <- nrow(b_mat); n_t <- ncol(b_mat)
  V <- t(vapply(seq_len(n_eau),
                function(i) compute_value_vector(b_mat[i, ], lam_mat[i, ], delta),
                numeric(n_t)))
  acquired <- rep(NA_integer_, n_eau)
  for (tnow in seq_len(n_t)) {
    avail <- which(is.na(acquired))
    if (length(avail) == 0) break
    sol <- solve_acquisition_ilp(V, cost_mat, budget, avail, periods = tnow:n_t)
    now <- sol$picks$eau_idx[sol$picks$period == tnow]
    acquired[now] <- tnow
  }
  acquired
}

# MYOPIC ILP — at each period freeze benefit and hazard at the CURRENT period's
# values for the whole remaining horizon (a fabricated stationary future), solve,
# implement only the current period, advance and observe the true next period.
# Survival is global under the frozen belief, so under a stationary scenario this
# equals the rolling policy exactly.
run_myopic_ilp <- function(b_mat, lam_mat, cost_mat, budget, delta = DELTA) {
  n_eau <- nrow(b_mat); n_t <- ncol(b_mat)
  acquired <- rep(NA_integer_, n_eau)
  for (tnow in seq_len(n_t)) {
    avail <- which(is.na(acquired))
    if (length(avail) == 0) break
    Vf <- matrix(0, n_eau, n_t)
    for (i in avail) {
      b_frozen   <- rep(b_mat[i, tnow],   n_t)   # believe current benefit persists
      lam_frozen <- rep(lam_mat[i, tnow], n_t)   # believe current hazard persists
      Vf[i, ] <- compute_value_vector(b_frozen, lam_frozen, delta)
    }
    sol <- solve_acquisition_ilp(Vf, cost_mat, budget, avail, periods = tnow:n_t)
    now <- sol$picks$eau_idx[sol$picks$period == tnow]
    acquired[now] <- tnow
  }
  acquired
}

# GREEDY HEURISTIC — status quo. Each period rank available EAUs by current
# benefit-cost ratio (scaled_abundance / cost) and acquire down the list while the
# budget allows. No foresight, no survival weighting, no future periods.
run_greedy <- function(b_mat, cost_mat, budget) {
  n_eau <- nrow(b_mat); n_t <- ncol(b_mat)
  acquired <- rep(NA_integer_, n_eau)
  for (tnow in seq_len(n_t)) {
    avail <- which(is.na(acquired))
    if (length(avail) == 0) break
    bcr <- b_mat[avail, tnow] / cost_mat[avail, tnow]
    ord <- avail[order(bcr, decreasing = TRUE)]
    budget_left <- budget[tnow]
    for (i in ord) {
      if (cost_mat[i, tnow] <= budget_left + 1e-9) {
        acquired[i] <- tnow
        budget_left <- budget_left - cost_mat[i, tnow]
      }
    }
  }
  acquired
}

# FULL-HORIZON single solve (validation reference for the rolling policy). Solves
# the whole multi-period program once with the true global V.
run_full_horizon <- function(b_mat, lam_mat, cost_mat, budget, delta = DELTA) {
  n_eau <- nrow(b_mat); n_t <- ncol(b_mat)
  V <- t(vapply(seq_len(n_eau),
                function(i) compute_value_vector(b_mat[i, ], lam_mat[i, ], delta),
                numeric(n_t)))
  sol <- solve_acquisition_ilp(V, cost_mat, budget,
                               avail = seq_len(n_eau), periods = seq_len(n_t))
  acquired <- rep(NA_integer_, n_eau)
  if (nrow(sol$picks) > 0) acquired[sol$picks$eau_idx] <- sol$picks$period
  acquired
}


# ── 5. Evaluation against the TRUE future (landscape-total welfare) ──────────####
#
# Scores ANY acquisition schedule against the TRUE benefit/hazard trajectories of a
# scenario. All three models are evaluated this way: the myopic schedule is scored
# on reality even though it was chosen under a fabricated future.
#
#   J = sum_i sum_t delta^t * b[i,t] * w[i,t],   w = 1 if protected by t else S[i,t]
#
# Also returns the per-period landscape duck trajectory (discounted and raw counts),
# the number of EAUs acquired, and the total acquisition spend.

evaluate_policy <- function(acquired, b_mat, lam_mat, cost_mat, delta = DELTA) {
  n_eau <- nrow(b_mat); n_t <- ncol(b_mat)
  S    <- compute_survival_matrix(lam_mat)
  disc <- delta^(0:(n_t - 1))

  # protection weight: 1 from acquisition period onward, else expected survival
  W <- S
  for (i in seq_len(n_eau)) {
    tau <- acquired[i]
    if (!is.na(tau)) W[i, tau:n_t] <- 1
  }

  weighted <- b_mat * W
  contrib  <- sweep(weighted, 2, disc, `*`)

  spend <- 0
  acq_idx <- which(!is.na(acquired))
  if (length(acq_idx) > 0)
    spend <- sum(cost_mat[cbind(acq_idx, acquired[acq_idx])])

  list(
    J                    = sum(contrib),
    per_period_discounted = colSums(contrib),
    per_period_ducks      = colSums(weighted),
    n_acquired            = length(acq_idx),
    total_spend           = spend,
    acquired              = acquired
  )
}


# ── 6. Small synthetic-instance generator (used by 09__validation.R) ─────────####
#
# Produces a tiny, fully-controlled instance for brute-force and ordering checks.
# Decoupled benefit shapes and hazards, optional stationarity, severity-scaled.
make_synthetic_instance <- function(n_eau, n_t, seed = 1,
                                     stationary = FALSE, severity = 1) {
  set.seed(seed)
  base  <- runif(n_eau, 3, 20)
  cost  <- outer(runif(n_eau, 0.5, 3.0), 1.02^(seq_len(n_t) - 1) * 10) # ~2%/yr
  shape <- sample(0:3, n_eau, replace = TRUE)

  b <- matrix(0, n_eau, n_t)
  for (i in seq_len(n_eau)) {
    x <- seq(0, 1, length.out = n_t)
    f <- switch(as.character(shape[i]),
      "0" = 1 + severity * 0.8 * x,                              # rising
      "1" = 1 - severity * 0.4 * x,                              # falling
      "2" = 1 + severity * 0.7 * sin(pi * x),                    # peak mid
      "3" = 1 + severity * 0.6 * pmax(x - 0.5, 0) * 2)           # late surge
    b[i, ] <- pmax(base[i] * (if (stationary) rep(1, n_t) else f), 0.05)
  }

  base_lam <- runif(n_eau, 0.02, 0.18)
  lam <- pmin(matrix(base_lam, n_eau, n_t) * (if (stationary) 1 else severity), 0.6)
  lam[, n_t] <- 0   # terminal period: no subsequent decision
  list(b = b, lam = lam, cost = cost)
}


# Brute-force optimum over a tiny instance: enumerate every assignment of each
# parcel to a period-or-never and return the best feasible objective. Only tractable
# for very small n_eau (used to certify the ILP formulation).
brute_force_optimum <- function(V_mat, cost_mat, budget) {
  n_eau <- nrow(V_mat); n_t <- ncol(V_mat)
  choices <- c(seq_len(n_t), 0L)   # 0 = never
  best <- -Inf
  grid <- expand.grid(rep(list(choices), n_eau))
  for (r in seq_len(nrow(grid))) {
    assign <- as.integer(grid[r, ])
    spend <- numeric(n_t); obj <- 0
    for (i in seq_len(n_eau)) if (assign[i] > 0) {
      spend[assign[i]] <- spend[assign[i]] + cost_mat[i, assign[i]]
      obj <- obj + V_mat[i, assign[i]]
    }
    if (all(spend <= budget + 1e-9) && obj > best) best <- obj
  }
  best
}

cat("07__ilp_core.R loaded: DELTA =", DELTA,
    "| BUDGET_EAUS_PER_PERIOD =", BUDGET_EAUS_PER_PERIOD,
    "| scenarios:", paste(names(SCENARIOS), collapse = ", "), "\n")
