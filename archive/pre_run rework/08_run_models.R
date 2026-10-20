# ══════════════════════════════════════════════════════════════════════════════
# 08__run_models.R  —  Run greedy / myopic / rolling across all scenarios
# ══════════════════════════════════════════════════════════════════════════════
#
# Sources the core engine, then for every scenario builds the (benefit, hazard,
# cost) matrices, runs all three policies, and scores each against that scenario's
# TRUE future. 

#Saves a tidy results table, per-period landscape trajectories, AND
# the per-parcel acquisition schedule (joined to EAU coordinates, with a per-parcel
# `prevented_loss` = TRUE-future discounted value_added contribution) so downstream
# reporting/mapping never has to re-solve. 

#The schedule is written both as output_data/acquisition_schedule_spatial.csv and 
#embedded in model_results.rds.

# RUN 09__validation.R FIRST (with Gurobi active) to confirm the formulation passes
# its correctness checks.

# ── TWO COMPLEMENTARY OUTCOME METRICS ─────────────────────────────────────────

# This script reports gaps on both total J (the true landscape outcome) and
# value_added (difference between welfare achieved under do-nothing baseline). 
# ══════════════════════════════════════════════════════════════════════════════

source("07_ilp_core.R")
if (!exists(".SETUP_DONE")) source("00_setup.R")

# ── Reproducibility ───────────────────────────────────────────────────────────####
# A capped solve under multi-threaded search is not bit-reproducible run-to-run
# (parallel workers may return different members of a co-optimal tie). For the
# numbers that go into the thesis, run deterministically (single-thread, fixed seed)
# so the figures are exactly reproducible. Set FALSE for a faster, multi-threaded
# pass (results then carry ~0.3%-scale run-to-run wobble — fine for exploration)

REPRODUCIBLE <- TRUE
if (REPRODUCIBLE) SOLVER_THREADS <- 1L

cat(sprintf("Run config: spend_down = '%s' | budget = %d median-EAUs/period | delta = %.2f | reproducible = %s\n",
            SPEND_DOWN_MODE, BUDGET_EAUS_PER_PERIOD, DELTA, REPRODUCIBLE))
if (identical(SPEND_DOWN_MODE, "off"))
  cat("  (spend_down off — reproduces the pre-spend-down numbers.)\n")

# ── 0. Load the data panel ────────────────────────────────────────────────────####
data_panel <- readRDS(file.path(DIR_DERIVED, "data_panel.rds"))

# EAU -> (wmd_id, x_coord, y_coord) lookup, used to make the persisted schedule a
# complete, mappable artifact 
EAU_LOOKUP <- file.path(DIR_DERIVED, "eau_wmd_lookup.rds")

OUT_DIR <- DIR_OUT
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

MODELS <- c("greedy", "myopic", "rolling")

# Solve-log accumulator: the solver in 07 appends one row per call when this exists.
.SOLVE_LOG <- list()
.SOLVE_TAG <- NA_character_


# ── 1. Run every model on every scenario ──────────────────────────────────────####

results_rows  <- list()   # one row per (scenario, model)
traj_rows     <- list()   # long: (scenario, model, year, ducks, discounted)
sched_rows    <- list()   # long: (scenario, model, eau_id, acquired_period/year, prevented_loss)

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

  # TRUE-future per-parcel value matrix V_true[i, tau] = discounted prevented loss if EAU i
  # is acquired at period tau, scored on the real b/lam (global survival from t=0). This is
  # the same global V the rolling policy optimises; we index it at each policy's REALISED
  # acquisition period to get per-parcel prevented_loss (sums to value_added — checked below).
  V_true <- t(vapply(seq_len(nrow(sc$b)),
                     function(i) compute_value_vector(sc$b[i, ], sc$lam[i, ], DELTA),
                     numeric(ncol(sc$b))))

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
    # per-parcel acquisition schedule: acquired[i] is the period (1..n_t) EAU i was
    # bought, or NA. Row order is sc$eau_ids (= matrix row order); map period -> year.
    # prevented_loss[i] = TRUE-future discounted prevented loss at the realised period
    # (V_true indexed at acquired[i]; 0 if never acquired). By construction these sum to the
    # model's value_added (= J - J_baseline); assert it as an internal correctness check.
    acq_vec <- schedules[[m]]
    pl  <- numeric(length(acq_vec))
    got <- which(!is.na(acq_vec))
    if (length(got)) pl[got] <- V_true[cbind(got, acq_vec[got])]
    va <- e$J - J_baseline
    stopifnot(abs(sum(pl) - va) <= 1e-6 * max(1, abs(va)))
    sched_rows[[length(sched_rows) + 1]] <- data.frame(
      scenario        = sc_name,
      model           = m,
      eau_id          = sc$eau_ids,
      acquired_period = acq_vec,
      acquired_year   = ifelse(is.na(acq_vec), NA_integer_, sc$years[acq_vec]),
      prevented_loss  = pl,
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

# ── 1b. Assemble the per-parcel schedule and attach EAU coordinates ───────────####
# Option A: join wmd_id / x_coord / y_coord here so the persisted schedule is a
# complete, mappable artifact and downstream reporting stays a pure CSV reader.
schedule <- do.call(rbind, sched_rows)
rownames(schedule) <- NULL
if (file.exists(EAU_LOOKUP)) {
  eau_lookup <- readRDS(EAU_LOOKUP)              # eau_id, wmd_id, x_coord, y_coord, ...
  schedule <- merge(schedule,
                    eau_lookup[, c("eau_id", "wmd_id", "x_coord", "y_coord")],
                    by = "eau_id", all.x = TRUE, sort = FALSE)
} else {
  warning("EAU lookup '", EAU_LOOKUP, "' not found; writing schedule WITHOUT coordinates ",
          "(maps will need them joined later).")
  schedule$wmd_id <- NA_character_
  schedule$x_coord <- NA_real_
  schedule$y_coord <- NA_real_
}
schedule <- schedule[order(schedule$scenario, schedule$model, schedule$eau_id),
                     c("scenario", "model", "eau_id", "wmd_id",
                       "x_coord", "y_coord", "acquired_period", "acquired_year",
                       "prevented_loss")]
rownames(schedule) <- NULL


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
             schedule = schedule,
             convergence = convergence, solve_log = log_df,
             params = list(delta = DELTA,
                           budget_eaus_per_period = BUDGET_EAUS_PER_PERIOD,
                           solver_time_limit = SOLVER_TIME_LIMIT,
                           solver_mip_gap = SOLVER_MIP_GAP,
                           solver_threads = SOLVER_THREADS,
                           spend_down_mode = SPEND_DOWN_MODE,
                           reproducible = REPRODUCIBLE)),
        file.path(OUT_DIR, "model_results.rds"))
write.csv(results,      file.path(OUT_DIR, "model_results.csv"),      row.names = FALSE)
write.csv(trajectories, file.path(OUT_DIR, "model_trajectories.csv"), row.names = FALSE)
write.csv(schedule,     file.path(OUT_DIR, "acquisition_schedule_spatial.csv"), row.names = FALSE)
if (!is.null(convergence))
  write.csv(convergence, file.path(OUT_DIR, "solver_convergence.csv"), row.names = FALSE)

cat("\nSaved:",
    "\n  ", file.path(OUT_DIR, "model_results.rds"),
    "\n  ", file.path(OUT_DIR, "model_results.csv"),
    "\n  ", file.path(OUT_DIR, "model_trajectories.csv"),
    "\n  ", file.path(OUT_DIR, "acquisition_schedule_spatial.csv"),
    if (!is.null(convergence)) paste0("\n   ", file.path(OUT_DIR, "solver_convergence.csv")) else "",
    "\n")
