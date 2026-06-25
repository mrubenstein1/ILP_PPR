# ══════════════════════════════════════════════════════════════════════════════
# 14__run_models_capped.R
#   The headline read: realized value_added gaps for greedy/myopic/rolling across
#   all five scenarios, with a per-solve solver-convergence log.
# ══════════════════════════════════════════════════════════════════════════════
#
# WHY THIS SCRIPT
# ------------------------------------------------------------------------------
# This is 08__run_models.R's logic, run with (a) OutputFlag = 0 and (b) a per-solve
# TimeLimit, plus a global log that records the status/gap/runtime of EVERY solver
# call. Two outputs come out of one run:
#
#   1. RESULTS table  — J, value_added, and gap-vs-rolling on both total J and
#      value_added, per (scenario, model). THIS IS THE HEADLINE: whether myopic and
#      rolling separate on the climate scenarios, and by how much vs the ~0.3%
#      solver-noise floor we measured in Q4.
#
#   2. CONVERGENCE summary — per (scenario, model): how many solves hit the cap and
#      the worst gap. This tells us whether the cap is even BINDING on the climate
#      scenarios. If climate myopic solves close OPTIMAL inside the cap, the cap
#      isn't distorting those numbers and no separate climate Q4 ladder is needed.
#      If they hit the cap at non-trivial gaps, that flags exactly which
#      (scenario, period) to probe further.
#
# Script 07 is untouched: solve_acquisition_ilp is overridden LOCALLY with a
# faithful superset (identical model assembly + 1e6 cost scaling) that is quiet,
# capped, and logs diagnostics. The override keeps the 07 call signature, so the
# run_greedy / run_myopic_ilp / run_rolling_horizon loops run unmodified.
#
# READ THE CAVEATS at the bottom before treating any single number as reportable.
# ══════════════════════════════════════════════════════════════════════════════

source("07__ilp_core.R")

# ══ PARAMETERS ═══════════════════════════════════════════════════════════════════
TIME_LIMIT <- 60     # seconds per INDIVIDUAL solve (myopic/rolling call ~9 each).
                     # 60 = the incumbent-saturation point we measured on stationary.
                     # For a faster first pass, drop to 30 (≈1.5% realized error there).
THREADS    <- 0      # 0 = all cores (production-style, best incumbent per wall-second).
                     # Set to 1 for a deterministic, reproducible run (slower).
SEED       <- 1      # only perturbs Gurobi tie-breaking; mild reproducibility help.
OUT_SUFFIX <- "_capped"
# ══════════════════════════════════════════════════════════════════════════════

OUT_DIR <- "output_data"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ── Global solve log + a context tag the override stamps onto each solve ─────────
.solve_log <- new.env()
.solve_log$rows <- list()
CURRENT_TAG <- "untagged"

# ── Overridden solver: faithful superset of 07 (quiet, capped, logged) ──────────
#    Same (A, obj, sense, rhs, vtype) assembly and 1e6 cost scaling as 07; same
#    return shape so the policy loops are unaffected.
solve_acquisition_ilp <- function(V_mat, cost_mat, budget, avail, periods,
                                   relax = FALSE) {
  if (length(avail) == 0 || length(periods) == 0)
    return(list(picks = data.frame(eau_idx = integer(0), period = integer(0)),
                objval = 0, x = numeric(0), grid = NULL))

  grid  <- expand.grid(eau_idx = avail, period = periods)
  nvar  <- nrow(grid)
  obj   <- V_mat[cbind(grid$eau_idx, grid$period)]
  costs <- cost_mat[cbind(grid$eau_idx, grid$period)]

  n_av <- length(avail); n_pe <- length(periods)
  row_once <- match(grid$eau_idx, avail)
  row_budg <- n_av + match(grid$period, periods)

  if (!requireNamespace("gurobi", quietly = TRUE))
    stop("The 'gurobi' R package is required (activate an academic licence).")

  cost_scale <- 1e6
  costs_s    <- costs / cost_scale

  A <- Matrix::sparseMatrix(
    i = c(row_once, row_budg), j = c(seq_len(nvar), seq_len(nvar)),
    x = c(rep(1, nvar), costs_s), dims = c(n_av + n_pe, nvar)
  )
  sense <- rep("<=", n_av + n_pe)
  rhs   <- c(rep(1, n_av), budget[periods] / cost_scale)

  model <- list(A = A, obj = obj, sense = sense, rhs = rhs,
                lb = rep(0, nvar), ub = rep(1, nvar),
                vtype = if (relax) "C" else "B", modelsense = "max")

  res <- gurobi::gurobi(model, params = list(
    OutputFlag = 0, TimeLimit = TIME_LIMIT, Threads = THREADS, Seed = SEED
  ))

  gx <- function(v, d = NA_real_) if (is.null(v)) d else v

  # ---- log this solve ----
  .solve_log$rows[[length(.solve_log$rows) + 1L]] <- data.frame(
    tag        = CURRENT_TAG,
    n_avail    = length(avail),
    n_periods  = length(periods),
    nvar       = nvar,
    status     = gx(res$status, "UNKNOWN"),
    gap_pct    = 100 * gx(res$mipgap, NA_real_),
    runtime_s  = round(gx(res$runtime, 0), 1),
    stringsAsFactors = FALSE
  )

  x <- res$x
  chosen <- if (relax) integer(0) else which(x > 0.5)
  list(picks = grid[chosen, , drop = FALSE], objval = res$objval, x = x, grid = grid)
}


# ══ DRIVER — 08 logic, with tagging before each policy ══════════════════════════
data_panel <- readRDS("input_data/data_panel.rds")
MODELS <- c("greedy", "myopic", "rolling")

results_rows <- list()
traj_rows    <- list()

for (sc_name in names(SCENARIOS)) {
  cat(sprintf("\n── Scenario: %-11s ──  (this can take a few minutes)\n", sc_name))
  sc  <- build_scenario_matrices(data_panel, SCENARIOS[[sc_name]])
  bud <- make_budget(sc$cost)

  J_baseline <- evaluate_policy(rep(NA_integer_, nrow(sc$b)),
                                sc$b, sc$lam, sc$cost)$J

  CURRENT_TAG <<- paste0(sc_name, "/greedy")
  sched_greedy  <- run_greedy(sc$b, sc$cost, bud)
  CURRENT_TAG <<- paste0(sc_name, "/myopic")
  sched_myopic  <- run_myopic_ilp(sc$b, sc$lam, sc$cost, bud)
  CURRENT_TAG <<- paste0(sc_name, "/rolling")
  sched_rolling <- run_rolling_horizon(sc$b, sc$lam, sc$cost, bud)

  schedules <- list(greedy = sched_greedy, myopic = sched_myopic, rolling = sched_rolling)
  evals     <- lapply(schedules, evaluate_policy,
                      b_mat = sc$b, lam_mat = sc$lam, cost_mat = sc$cost)
  J_rolling <- evals$rolling$J

  for (m in MODELS) {
    e <- evals[[m]]
    results_rows[[length(results_rows) + 1]] <- data.frame(
      scenario            = sc_name,
      model               = m,
      J                   = e$J,
      J_baseline          = J_baseline,
      value_added         = e$J - J_baseline,
      n_acquired          = e$n_acquired,
      total_spend         = e$total_spend,
      gap_J_pct           = 100 * (J_rolling - e$J) / J_rolling,
      gap_value_added_pct = if (abs(J_rolling - J_baseline) < 1e-9) NA_real_
                            else 100 * (J_rolling - e$J) / (J_rolling - J_baseline),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %-8s  J = %14.1f   value_added = %12.1f   n = %3d\n",
                m, e$J, e$J - J_baseline, e$n_acquired))
  }
}

results <- do.call(rbind, results_rows); rownames(results) <- NULL


# ══ RESULTS table (the headline) ════════════════════════════════════════════════
cat("\n\n══ RESULTS: realized landscape welfare and gap vs rolling ══\n")
print(
  results %>%
    select(scenario, model, J, value_added, gap_J_pct, gap_value_added_pct,
           n_acquired, total_spend) %>%
    mutate(across(c(J, value_added, total_spend), ~ round(.x, 1)),
           across(c(gap_J_pct, gap_value_added_pct), ~ round(.x, 4))),
  row.names = FALSE
)

cat("\n══ Ordering check (rolling should be >= myopic; watch for inversions) ══\n")
ord <- results %>%
  select(scenario, model, value_added) %>%
  tidyr::pivot_wider(names_from = model, values_from = value_added) %>%
  mutate(rolling_minus_myopic = round(rolling - myopic, 2),
         rolling_minus_greedy = round(rolling - greedy, 2),
         inversion = rolling < myopic - 1e-6)
print(as.data.frame(ord), row.names = FALSE)


# ══ CONVERGENCE summary (is the cap binding? where?) ════════════════════════════
safe_max <- function(x) { x <- x[is.finite(x)]; if (length(x)) max(x) else NA_real_ }

log_df <- do.call(rbind, .solve_log$rows)
conv <- log_df %>%
  group_by(tag) %>%
  summarise(
    n_solves        = dplyr::n(),
    n_hit_timelimit = sum(status == "TIME_LIMIT"),
    worst_gap_pct   = round(safe_max(gap_pct), 3),
    total_runtime_s = round(sum(runtime_s), 1),
    .groups = "drop"
  )

cat("\n\n══ SOLVER CONVERGENCE per (scenario / model) ══\n")
cat("  n_hit_timelimit = 0  => every solve closed OPTIMAL inside the cap (trustworthy).\n")
cat("  n_hit_timelimit > 0  => the cap is binding; read worst_gap_pct next to the\n")
cat("                          value_added gaps above to judge if it distorts them.\n\n")
print(as.data.frame(conv), row.names = FALSE)

# Per-iteration detail, so you can see WHICH period of myopic stalls on each scenario
cat("\n── Per-solve detail (myopic only; n_periods shows the re-solve iteration) ──\n")
myo_detail <- log_df %>%
  filter(grepl("/myopic$", tag)) %>%
  mutate(gap_pct = round(gap_pct, 3)) %>%
  select(tag, n_avail, n_periods, status, gap_pct, runtime_s)
print(as.data.frame(myo_detail), row.names = FALSE)


# ══ SAVE ════════════════════════════════════════════════════════════════════════
write.csv(results, file.path(OUT_DIR, paste0("model_results", OUT_SUFFIX, ".csv")),
          row.names = FALSE)
write.csv(conv,    file.path(OUT_DIR, paste0("solver_convergence", OUT_SUFFIX, ".csv")),
          row.names = FALSE)
write.csv(log_df,  file.path(OUT_DIR, paste0("solve_log", OUT_SUFFIX, ".csv")),
          row.names = FALSE)
saveRDS(list(results = results, convergence = conv, solve_log = log_df,
             params = list(time_limit = TIME_LIMIT, threads = THREADS, seed = SEED,
                           delta = DELTA, budget_eaus_per_period = BUDGET_EAUS_PER_PERIOD)),
        file.path(OUT_DIR, paste0("model_results", OUT_SUFFIX, ".rds")))

cat(sprintf("\nSaved:\n  %s\n  %s\n  %s\n  %s\n",
            file.path(OUT_DIR, paste0("model_results", OUT_SUFFIX, ".csv")),
            file.path(OUT_DIR, paste0("solver_convergence", OUT_SUFFIX, ".csv")),
            file.path(OUT_DIR, paste0("solve_log", OUT_SUFFIX, ".csv")),
            file.path(OUT_DIR, paste0("model_results", OUT_SUFFIX, ".rds"))))


# ══════════════════════════════════════════════════════════════════════════════
# HOW TO READ THIS — and the caveats that matter
# ──────────────────────────────────────────────────────────────────────────────
# 1. HEADLINE = gap_value_added_pct for myopic (and greedy) on the CLIMATE scenarios
#    (rcp45/85 × wet/dry). That is the foresight signal. Compare it to ~0.3%: the
#    realized solver-noise floor a 60s cap left on the (worst-case) stationary t=1
#    solve in Q4.
#       - signal comfortably > ~0.3%  -> a 60s cap is fine; report the residual gap
#         as a small uncertainty band and move on.
#       - signal <= ~1%               -> the cap competes with the effect; converge
#         the solves (longer cap) or denoise before reporting.
#
# 2. CAP-BINDING CHECK = the convergence table. If the climate scenarios show
#    n_hit_timelimit = 0 for myopic, the cap is NOT binding there: those value_added
#    numbers are at the true optimum and a separate climate Q4 ladder is unnecessary.
#    (This is the Q3 prediction — real, parcel-varying climate hazards give the
#    solver value spread to discriminate on, so it closes. The per-solve detail
#    table confirms it iteration by iteration.)
#
# 3. KNOWN DISTORTIONS, stated honestly:
#    - The cap pushes MYOPIC's realized value slightly DOWN (it returns a marginally
#      sub-optimal frozen schedule), which can INFLATE the apparent rolling-over-
#      myopic gap. So a measured gap is, if anything, an upper bound on the truth:
#      a tiny measured gap is therefore strong evidence the true gap is tiny.
#    - With THREADS = 0 + a TimeLimit, a capped solve is not bit-reproducible run to
#      run; expect ~0.3%-scale wobble on any capped myopic number. For a REPORTABLE
#      figure, either set THREADS = 1 (deterministic) or confirm stability across a
#      couple of runs. For this diagnostic read, one THREADS = 0 run is fine.
#    - The STATIONARY row is the null; its myopic t>=2 solves optimize a ~zero
#      objective (epsilon hazard) and will stall at the cap. That is expected and
#      does not affect the climate headline — but it is why the stationary
#      convergence row will look the worst.
# ══════════════════════════════════════════════════════════════════════════════
