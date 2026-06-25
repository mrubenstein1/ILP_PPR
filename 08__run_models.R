# ══════════════════════════════════════════════════════════════════════════════
# 08__run_models.R  —  Run greedy / myopic / rolling across all scenarios
# ══════════════════════════════════════════════════════════════════════════════
#
# Sources the core engine, then for every scenario builds the (benefit, hazard,
# cost) matrices, runs all three policies, and scores each against that scenario's
# TRUE future. Saves a tidy results table and per-period landscape trajectories.
#
# RUN 09__validation.R FIRST (with Gurobi active) to confirm the formulation passes
# its correctness checks on your machine. The model logic here was validated in a
# separate solver harness; this script implements that verified logic in R/Gurobi.
#
# ── TWO COMPLEMENTARY OUTCOME METRICS ─────────────────────────────────────────
# Because we acquire only a tiny fraction of 841 EAUs, the do-nothing landscape
# total J_baseline dominates J, so gaps measured on the TOTAL J look small. The
# quantity a manager actually controls is the conservation value created by
# acquisition,
#       value_added = J - J_baseline = sum_i V[i, tau_i].
# We report gaps on BOTH: total J (the true landscape outcome, the headline) and
# value_added (where model choice actually shows up). Under the current low-
# conversion data both gaps are modest and grow with conversion magnitude; the
# value_added view makes the model ordering visible where the total-J view masks it.
# ══════════════════════════════════════════════════════════════════════════════

source("07__ilp_core.R")

# ── Reproducibility ───────────────────────────────────────────────────────────####
# A capped solve under multi-threaded search is not bit-reproducible run-to-run
# (parallel workers may return different members of a co-optimal tie). For the
# numbers that go into the thesis, run deterministically (single-thread, fixed seed)
# so the figures are exactly reproducible. Set FALSE for a faster, multi-threaded
# pass (results then carry ~0.3%-scale run-to-run wobble — fine for exploration,
# and far below the ~13% effect, but pin this TRUE for anything you report).
REPRODUCIBLE <- TRUE
if (REPRODUCIBLE) SOLVER_THREADS <- 1L

# ── 0. Load the data panel ────────────────────────────────────────────────────####
data_panel <- readRDS("input_data/data_panel.rds")

OUT_DIR <- "output_data"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

MODELS <- c("greedy", "myopic", "rolling")

# Solve-log accumulator: the solver in 07 appends one row per call when this exists.
.SOLVE_LOG <- list()
.SOLVE_TAG <- NA_character_


# ── 1. Run every model on every scenario ──────────────────────────────────────####

results_rows  <- list()   # one row per (scenario, model)
traj_rows     <- list()   # long: (scenario, model, year, ducks, discounted)

for (sc_name in names(SCENARIOS)) {
  cat("\n── Scenario:", sc_name, "──────────────────────────────\n")
  sc  <- build_scenario_matrices(data_panel, SCENARIOS[[sc_name]])
  bud <- make_budget(sc$cost)

  # do-nothing reference (all NA = nothing acquired)
  base_eval <- evaluate_policy(rep(NA_integer_, nrow(sc$b)),
                               sc$b, sc$lam, sc$cost)
  J_baseline <- base_eval$J

  # run each policy (tag each so the solve log attributes solves to scenario/model;
  # greedy uses no solver, so it produces no log rows)
  .SOLVE_TAG <- paste0(sc_name, "/greedy")
  sched_greedy  <- run_greedy(sc$b, sc$cost, bud)
  .SOLVE_TAG <- paste0(sc_name, "/myopic")
  sched_myopic  <- run_myopic_ilp(sc$b, sc$lam, sc$cost, bud)
  .SOLVE_TAG <- paste0(sc_name, "/rolling")
  sched_rolling <- run_rolling_horizon(sc$b, sc$lam, sc$cost, bud)
  schedules <- list(greedy = sched_greedy, myopic = sched_myopic, rolling = sched_rolling)

  # evaluate each against the scenario's TRUE future
  evals <- lapply(schedules, evaluate_policy,
                  b_mat = sc$b, lam_mat = sc$lam, cost_mat = sc$cost)
  J_rolling <- evals$rolling$J

  for (m in MODELS) {
    e <- evals[[m]]
    results_rows[[length(results_rows) + 1]] <- data.frame(
      scenario             = sc_name,
      model                = m,
      J                    = e$J,
      J_baseline           = J_baseline,
      value_added          = e$J - J_baseline,
      n_acquired           = e$n_acquired,
      total_spend          = e$total_spend,
      gap_J_pct            = 100 * (J_rolling - e$J) / J_rolling,
      gap_value_added_pct  = if (abs(J_rolling - J_baseline) < 1e-9) NA_real_
                             else 100 * (J_rolling - e$J) / (J_rolling - J_baseline),
      stringsAsFactors = FALSE
    )
    traj_rows[[length(traj_rows) + 1]] <- data.frame(
      scenario   = sc_name,
      model      = m,
      year       = sc$years,
      ducks      = e$per_period_ducks,
      discounted = e$per_period_discounted,
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %-8s  J = %12.1f   value_added = %10.1f   n = %3d   spend = %12.0f\n",
                m, e$J, e$J - J_baseline, e$n_acquired, e$total_spend))
  }
}

results      <- do.call(rbind, results_rows)
trajectories <- do.call(rbind, traj_rows)
rownames(results) <- NULL
rownames(trajectories) <- NULL


# ── 2. Summary tables ─────────────────────────────────────────────────────────####

cat("\n\n══ RESULTS: total landscape welfare J and gap vs rolling ══\n")
print(
  results %>%
    select(scenario, model, J, value_added, gap_J_pct, gap_value_added_pct,
           n_acquired, total_spend) %>%
    mutate(across(c(J, value_added, total_spend), ~ round(.x, 1)),
           across(c(gap_J_pct, gap_value_added_pct), ~ round(.x, 3))),
  row.names = FALSE
)

cat("\n══ Model ordering check (rolling should be >= myopic >= ... ) ══\n")
ord <- results %>%
  select(scenario, model, J) %>%
  tidyr::pivot_wider(names_from = model, values_from = J) %>%
  mutate(rolling_ge_myopic = rolling >= myopic - 1e-6,
         myopic_minus_greedy = myopic - greedy)
print(as.data.frame(ord), row.names = FALSE)


# ── 2b. Solver convergence summary ────────────────────────────────────────────####
# Shows, per (scenario / model), how many of the per-period solves reached the time
# cap and the worst optimality gap. n_hit_timelimit = 0 means every solve closed
# OPTIMAL inside the cap (numbers at the true optimum). Where the cap binds, read
# worst_gap_pct next to the value_added gaps above: it should be tiny (<0.05%) on the
# climate scenarios, confirming the cap does not distort the headline.
log_df <- if (length(.SOLVE_LOG)) do.call(rbind, .SOLVE_LOG) else NULL
if (!is.null(log_df)) {
  safe_max <- function(x) { x <- x[is.finite(x)]; if (length(x)) max(x) else NA_real_ }
  convergence <- log_df %>%
    group_by(tag) %>%
    summarise(n_solves        = dplyr::n(),
              n_hit_timelimit = sum(status == "TIME_LIMIT"),
              worst_gap_pct   = round(safe_max(gap_pct), 3),
              total_runtime_s = round(sum(runtime_s), 1),
              .groups = "drop") %>%
    as.data.frame()
  cat("\n══ Solver convergence (per scenario / model) ══\n")
  print(convergence, row.names = FALSE)
} else {
  convergence <- NULL
}


# ── 3. Save outputs ───────────────────────────────────────────────────────────####

saveRDS(list(results = results, trajectories = trajectories,
             convergence = convergence, solve_log = log_df,
             params = list(delta = DELTA,
                           budget_eaus_per_period = BUDGET_EAUS_PER_PERIOD,
                           solver_time_limit = SOLVER_TIME_LIMIT,
                           solver_mip_gap = SOLVER_MIP_GAP,
                           solver_threads = SOLVER_THREADS,
                           reproducible = REPRODUCIBLE)),
        file.path(OUT_DIR, "model_results.rds"))
write.csv(results,      file.path(OUT_DIR, "model_results.csv"),      row.names = FALSE)
write.csv(trajectories, file.path(OUT_DIR, "model_trajectories.csv"), row.names = FALSE)
if (!is.null(convergence))
  write.csv(convergence, file.path(OUT_DIR, "solver_convergence.csv"), row.names = FALSE)

cat("\nSaved:",
    "\n  ", file.path(OUT_DIR, "model_results.rds"),
    "\n  ", file.path(OUT_DIR, "model_results.csv"),
    "\n  ", file.path(OUT_DIR, "model_trajectories.csv"),
    if (!is.null(convergence)) paste0("\n   ", file.path(OUT_DIR, "solver_convergence.csv")) else "",
    "\n")
