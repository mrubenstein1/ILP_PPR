# ══════════════════════════════════════════════════════════════════════════════
# 10_extract_schedules.R  —  Export the per-parcel acquisition schedule for mapping
# ══════════════════════════════════════════════════════════════════════════════
#
# WHY THIS EXISTS.  08__run_models.R computes each policy's acquisition schedule
# (which EAU is bought in which period) but saves only the SUMMARY rows and the
# per-period trajectories — the schedules themselves are discarded. Maps need the
# schedules, so this script re-runs the three policies on every scenario and writes a
# tidy, spatially-keyed table:
#
#     output_data/acquisition_schedule_spatial.csv
#     columns: scenario, model, eau_id, wmd_id, x_coord, y_coord,
#              acquired_period (1..9 or NA), acquired_year (2020..2100 or NA)
#
# It re-solves, so it costs roughly one 08 run. It is set REPRODUCIBLE by default
# (single thread + fixed seed) so the schedule is stable and re-runnable; the spatial
# PATTERN is what the maps use and is robust to solver tie-breaking either way. The
# aggregate counts/spend it prints should match model_results.csv to within the
# ~0.3%-scale wobble noted in the diagnostic.
#
# HOW TO RUN.  Set CORE below to your spend-down core filename, then:  source(this).
# ══════════════════════════════════════════════════════════════════════════════

# ── set me ────────────────────────────────────────────────────────────────────
CORE    <- "07__ilp_core.R"   # <-- your spend-down core (e.g. "07_ilp_core_spenddown.R")
OUT_DIR <- "output_data"      # where the CSV is written
# ──────────────────────────────────────────────────────────────────────────────

source(CORE)
stopifnot("spend_down" %in% names(formals(solve_acquisition_ilp)))  # confirm spend-down core
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Reproducible schedule: single thread + fixed seed (overridable).
SOLVER_THREADS <- 1L
set.seed(1)

data_panel <- readRDS("input_data/data_panel.rds")
eau_lookup <- readRDS("input_data/eau_wmd_lookup.rds")   # eau_id, wmd_id, x_coord, y_coord, ...

# Matrix row order in 07 is sort(unique(eau_id)); the schedule index maps to this.
eau_ids <- sort(unique(data_panel$eau_id))
MODELS  <- c("greedy", "myopic", "rolling")

sched_rows <- list()
for (sc_name in names(SCENARIOS)) {
  cat("\n── Scenario:", sc_name, "──\n")
  sc  <- build_scenario_matrices(data_panel, SCENARIOS[[sc_name]])
  bud <- make_budget(sc$cost)
  stopifnot(length(eau_ids) == nrow(sc$b))

  schedules <- list(
    greedy  = run_greedy(sc$b, sc$cost, bud),
    myopic  = run_myopic_ilp(sc$b, sc$lam, sc$cost, bud),
    rolling = run_rolling_horizon(sc$b, sc$lam, sc$cost, bud)
  )

  for (m in MODELS) {
    acq <- schedules[[m]]                       # length n_eau; period 1..9 or NA
    cat(sprintf("   %-8s acquired %3d parcels\n", m, sum(!is.na(acq))))
    sched_rows[[length(sched_rows) + 1]] <- data.frame(
      scenario        = sc_name,
      model           = m,
      eau_id          = eau_ids,
      acquired_period = acq,
      acquired_year   = ifelse(is.na(acq), NA_integer_, sc$years[acq]),
      stringsAsFactors = FALSE
    )
  }
}

schedule <- do.call(rbind, sched_rows)

# Attach coordinates (centroids of the equal-area raster cells) for mapping.
schedule <- merge(
  schedule,
  eau_lookup[, c("eau_id", "wmd_id", "x_coord", "y_coord")],
  by = "eau_id", all.x = TRUE, sort = FALSE
)
schedule <- schedule[order(schedule$scenario, schedule$model, schedule$eau_id),
                     c("scenario", "model", "eau_id", "wmd_id",
                       "x_coord", "y_coord", "acquired_period", "acquired_year")]

out <- file.path(OUT_DIR, "acquisition_schedule_spatial.csv")
write.csv(schedule, out, row.names = FALSE)
cat("\nWrote", out, "—", nrow(schedule), "rows (",
    length(unique(schedule$eau_id)), "EAUs ×",
    length(MODELS), "models ×", length(SCENARIOS), "scenarios ).\n")
cat("Acquired-parcel counts per scenario/model:\n")
print(aggregate(acquired_period ~ scenario + model, schedule,
                FUN = function(z) sum(!is.na(z)), na.action = na.pass))
