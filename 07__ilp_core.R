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

# ── Solver control ────────────────────────────────────────────────────────────
# The myopic frozen-belief solves reach a near-optimal incumbent within seconds,
# then spend the remaining time only CERTIFYING optimality against the LP bound
# across a large band of co-valued schedules (a benign symmetry of the flattened
# frozen objective). Because each policy implements only the current period and
# re-optimises, certifying the discarded tail of the schedule does not change the
# enacted trajectory. We therefore cap each solve's runtime. Empirically (diagnostic
# run, script 14) only 1-2 of the 9 per-period climate solves reach a 60s cap, each
# at <=0.05% optimality gap; the resulting realised value_added is within <0.5% of
# the converged value, far below the ~13% rolling-vs-myopic effect. Truncation only
# ever biases the myopic incumbent slightly DOWN, so it cannot manufacture a
# rolling advantage — the comparison is conservative under the cap.
SOLVER_TIME_LIMIT <- 60      # seconds per solve (cap on the optimality-proof phase)
SOLVER_MIP_GAP    <- 1e-4    # Gurobi default; the time cap usually binds first
SOLVER_THREADS    <- 0L      # 0 = all cores (fast). Set 1L for bit-reproducible runs.
SOLVER_SEED       <- 1L      # affects only solver tie-breaking
SOLVER_OUTPUT     <- 0L      # 0 = quiet; 1 = full Gurobi search log

# ── Budget deployment (spend-down) ────────────────────────────────────────────
# Real managers deploy the whole per-period budget (it cannot be rolled over) and
# do not leave money idle merely because the next parcel looks low-value under their
# current belief. Each ILP solve therefore optimises acquisition value FIRST and then,
# as a STRICTLY LOWER-PRIORITY tiebreak among value-co-optimal plans, deploys the
# IMPLEMENTED period's budget as fully as possible. This removes a MIPGap-sensitive
# artefact — myopic leaving ~22% of its budget unspent because its frozen belief
# assigns the marginal parcels near-zero value, so the solver is indifferent to them —
# and replaces it with an explicit, defensible behavioural rule.
#
#   "off"   — value only. Reproduces the pre-spend-down numbers; used for clean
#             formulation validation (P1/P2/P5 in 09).
#   "spend" — secondary tiebreak maximises DOLLARS deployed in the implemented period
#             ("don't leave money on the table").
#   "count" — secondary tiebreak maximises NUMBER of parcels acquired (cheaper-leaning;
#             "more lottery tickets" on currently-low parcels that may rise in value).
#
# Greedy already deploys its budget by construction (it buys down its ranked list until
# nothing affordable remains), so this knob affects only the ILP policies — materially
# for myopic, negligibly for rolling (which already nearly exhausts the budget).
SPEND_DOWN_MODE   <- "spend"   # "off" | "spend" | "count"

# Max fraction of the primary (value) objective the spend tiebreak may sacrifice.
# Kept tiny so the value-optimal core is PINNED and only genuinely marginal parcels
# (value within this tolerance of the optimum) are added to fill the budget. With the
# core pinned, the tiebreak is pure ADDITION of non-negative-value parcels in the
# implemented period, so realised (true) J can only weakly improve — see header note.
SPEND_DOWN_RELTOL <- 1e-4
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

  # ── Stationary null correction ──────────────────────────────────────────────
  # The block above fills period-0 (2020) from the SHARED cross-RCP baseline row for
  # EVERY scenario. For the stationary null that anchor is wrong: the baseline 2020
  # hazard is the mean of the RCP4.5/8.5 2020->2030 transitions (> epsilon), so the
  # stationary hazard trajectory is NOT flat at 2020. That single non-stationary step
  # makes the myopic frozen belief (current values persist) disagree with reality and
  # breaks the intended null in which myopic must equal rolling when nothing changes.
  # Fix: set the stationary 2020 hazard to the scenario's own background floor epsilon
  # (recovered as its 2030 hazard, which script 05 sets to epsilon for all parcels).
  # Climate scenarios are untouched. NOTE: this flattens the HAZARD only; whether the
  # benefit trajectory is also flat (the other precondition of the null) is certified
  # separately by P3 in 09__validation.R.
  if (identical(as.character(scenario[["rcp"]]), "stationary") && n_t >= 2) {
    lam[, 1] <- lam[, 2]
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
                                   relax = FALSE,
                                   spend_down = SPEND_DOWN_MODE,
                                   implement_periods = periods) {
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
    sense      = sense,
    rhs        = rhs,
    lb         = rep(0, nvar),
    ub         = rep(1, nvar),
    vtype      = if (relax) "C" else "B",
    modelsense = "max"
  )

  # Primary objective is ALWAYS the acquisition value V. When spend-down is active
  # (integer solves only), attach a strictly lower-priority secondary objective that,
  # among value-co-optimal plans, maximises deployment of the IMPLEMENTED period's
  # budget. Targeting only `implement_periods` (the period(s) actually enacted — a
  # single period in the rolling/myopic loop) is deliberate: it forces the budget we
  # are about to spend to be deployed now, and avoids a perverse "defer to a pricier
  # later period" tie that a horizon-wide spend objective would create under cost
  # inflation. SPEND_DOWN_RELTOL pins the value-optimal core; only marginal parcels
  # are added.
  use_spend_down <- !relax && !identical(spend_down, "off")
  if (use_spend_down) {
    in_impl <- grid$period %in% implement_periods
    sec <- numeric(nvar)
    if (identical(spend_down, "spend")) {
      sec[in_impl] <- costs_s[in_impl]   # maximise dollars deployed (scaled units)
    } else if (identical(spend_down, "count")) {
      sec[in_impl] <- 1                  # maximise number of parcels acquired
    } else {
      stop("spend_down must be one of 'off', 'spend', 'count'.")
    }
    model$multiobj <- list(
      list(objn = obj, priority = 2L, weight = 1, reltol = SPEND_DOWN_RELTOL, abstol = 0),
      list(objn = sec, priority = 1L, weight = 1, reltol = 0,                 abstol = 0)
    )
  } else {
    model$obj <- obj
  }

  res <- gurobi::gurobi(model, params = list(
    OutputFlag = SOLVER_OUTPUT,
    TimeLimit  = SOLVER_TIME_LIMIT,
    MIPGap     = SOLVER_MIP_GAP,
    Threads    = SOLVER_THREADS,
    Seed       = SOLVER_SEED
  ))

  # Optional instrumentation: if a driver has created a global .SOLVE_LOG list, append
  # this solve's status / optimality gap / runtime. Pure side-effect — it does not
  # touch the model, the solution, or any returned value, and is inert when .SOLVE_LOG
  # does not exist (so solving is unaffected outside an instrumented run).
  if (exists(".SOLVE_LOG", envir = .GlobalEnv)) {
    .gx  <- function(v, d = NA_real_) if (is.null(v)) d else v
    .tag <- if (exists(".SOLVE_TAG", envir = .GlobalEnv))
              get(".SOLVE_TAG", envir = .GlobalEnv) else NA_character_
    .lg  <- get(".SOLVE_LOG", envir = .GlobalEnv)
    .lg[[length(.lg) + 1L]] <- data.frame(
      tag = .tag, n_avail = length(avail), n_periods = length(periods), nvar = nvar,
      status = .gx(res$status, "UNKNOWN"), gap_pct = 100 * .gx(res$mipgap, NA_real_),
      runtime_s = round(.gx(res$runtime, 0), 1), stringsAsFactors = FALSE)
    assign(".SOLVE_LOG", .lg, envir = .GlobalEnv)
  }
  # ---------------------------------------------------------------------------

  x <- res$x
  chosen <- if (relax) integer(0) else which(x > 0.5)
  picks  <- grid[chosen, , drop = FALSE]
  # Under multi-objective, res$objval is a vector (one per objective, in list order);
  # element 1 is the primary acquisition value. Return that scalar so downstream
  # comparisons (validation P1/P5, reporting) are unaffected by the spend tiebreak.
  prim_obj <- if (length(res$objval) > 1L) res$objval[[1L]] else res$objval
  list(picks = picks, objval = prim_obj, x = x, grid = grid)
}


# ── 4. The three policies (+ full-horizon reference) ─────────────────────────####
#
# Each returns an integer vector `acquired` of length n_eau giving the period index
# (1..n_t) at which each EAU is acquired, or NA if never acquired.

# ROLLING HORIZON — re-solve each period using the TRUE projected future (global V),
# implement only the current period, advance. Under deterministic foresight this
# reproduces the full-horizon optimum (see run_full_horizon); the loop form is kept
# for methodological fidelity and so stochastic extensions can slot in.
run_rolling_horizon <- function(b_mat, lam_mat, cost_mat, budget, delta = DELTA,
                                spend_down = SPEND_DOWN_MODE) {
  n_eau <- nrow(b_mat); n_t <- ncol(b_mat)
  V <- t(vapply(seq_len(n_eau),
                function(i) compute_value_vector(b_mat[i, ], lam_mat[i, ], delta),
                numeric(n_t)))
  acquired <- rep(NA_integer_, n_eau)
  for (tnow in seq_len(n_t)) {
    avail <- which(is.na(acquired))
    if (length(avail) == 0) break
    sol <- solve_acquisition_ilp(V, cost_mat, budget, avail, periods = tnow:n_t,
                                 spend_down = spend_down, implement_periods = tnow)
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
run_myopic_ilp <- function(b_mat, lam_mat, cost_mat, budget, delta = DELTA,
                           spend_down = SPEND_DOWN_MODE) {
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
    sol <- solve_acquisition_ilp(Vf, cost_mat, budget, avail, periods = tnow:n_t,
                                 spend_down = spend_down, implement_periods = tnow)
    now <- sol$picks$eau_idx[sol$picks$period == tnow]
    acquired[now] <- tnow
  }
  acquired
}

# GREEDY HEURISTIC — status quo. Each period rank available EAUs by current
# benefit-cost ratio (scaled_abundance / cost) and acquire down the list while the
# budget allows. No foresight, no survival weighting, no future periods. (Greedy
# already deploys its budget by construction — it buys until nothing affordable
# remains — so the SPEND_DOWN_MODE knob does not apply to it.)
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
run_full_horizon <- function(b_mat, lam_mat, cost_mat, budget, delta = DELTA,
                             spend_down = SPEND_DOWN_MODE) {
  n_eau <- nrow(b_mat); n_t <- ncol(b_mat)
  V <- t(vapply(seq_len(n_eau),
                function(i) compute_value_vector(b_mat[i, ], lam_mat[i, ], delta),
                numeric(n_t)))
  # Single full-horizon solve: every period is "enacted", so implement_periods is the
  # whole horizon (the default). NB validation calls this with spend_down = "off" so the
  # rolling-vs-full equivalence (P2) is checked on the pure value optimum, independent of
  # any spend tiebreak (which a single solve and the per-period loop break differently).
  sol <- solve_acquisition_ilp(V, cost_mat, budget,
                               avail = seq_len(n_eau), periods = seq_len(n_t),
                               spend_down = spend_down)
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
