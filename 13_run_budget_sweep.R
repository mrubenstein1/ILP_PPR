# ══════════════════════════════════════════════════════════════════════════════
# 13_run_budget_sweep.R — budget sensitivity sweep (land cover set in 00_config.R)
# ══════════════════════════════════════════════════════════════════════════════
#
# Re-runs the decision models (08) and the per-run assets (10) at several budget
# levels, each in its own output folder, runs/<LC_SCENARIO>_bXX/. Land-cover
# processing (01-06) is NOT rerun: every budget run reads the derived data of the
# land-cover run, runs/<LC_SCENARIO>/derived_data/, because the budget does not
# enter scripts 01-06.
#
# PREREQUISITE: the land-cover run must already exist. From scratch, first run
#   Rscript 00_run_all.R
# with PPR_BUDGET unset (this is the canonical run).
#
# USAGE (from the project root):
#   Rscript 13_run_budget_sweep.R              run the budgets in BUDGETS below
#   Rscript 13_run_budget_sweep.R 2 10         run only these budgets
#   Rscript 13_run_budget_sweep.R --dry-run    show what would run; run nothing
#   Rscript 13_run_budget_sweep.R --force      rerun budgets even if complete
#
# On a Mac, stop the machine sleeping during the multi-hour sweep:
#   caffeinate -i Rscript 13_run_budget_sweep.R
#
# Each budget runs  PPR_BUDGET=<b> Rscript 00_run_all.R 08  (i.e. 08 then 10), so a
# sweep run is identical to running that command by hand. Finished budgets are
# skipped, so an interrupted sweep restarts with the same command.
#
# Console output of each run goes to logs/budget_sweep/<RUN_ID>.log.
# Follow a run live from another terminal with:
#   tail -f logs/budget_sweep/<RUN_ID>.log
#
# If budget 5 is in the sweep, its results are compared with the canonical run
# (runs/<LC_SCENARIO>/output_data/model_results.csv) at the end. Small
# differences are expected (time-capped solves are not bit-reproducible); larger
# ones are flagged. Base R only: no packages needed to run this driver.
# ══════════════════════════════════════════════════════════════════════════════

# ── Settings ──────────────────────────────────────────────────────────────────
# Budget levels, in median-priced EAUs per period. Budget 5 runs first so the
# check against the canonical run is available after the first run.
BUDGETS <- c(5, 2, 10)

# Flag a b05-vs-canonical difference larger than this fraction of value_added
# (climate scenarios only; the stationary null is too small for a relative test).
CANON_TOL <- 0.005

LOG_DIR <- file.path("logs", "budget_sweep")

# ── Command-line arguments ────────────────────────────────────────────────────
args    <- commandArgs(trailingOnly = TRUE)
opts    <- args[grepl("^--", args)]
nums    <- args[!grepl("^--", args)]
DRY_RUN <- "--dry-run" %in% opts
FORCE   <- "--force"   %in% opts

unknown <- setdiff(opts, c("--dry-run", "--force"))
if (length(unknown))
  stop("Unknown option(s): ", paste(unknown, collapse = ", "),
       ". Use --dry-run and/or --force.", call. = FALSE)

if (length(nums)) {
  BUDGETS <- suppressWarnings(as.numeric(nums))
  if (anyNA(BUDGETS) || any(BUDGETS <= 0))
    stop("Budgets must be positive numbers, e.g. ",
         "'Rscript 13_run_budget_sweep.R 2 10'.", call. = FALSE)
}

# ── Checks ────────────────────────────────────────────────────────────────────
if (!file.exists("00_run_all.R") || !file.exists("00_config.R"))
  stop("00_run_all.R / 00_config.R not found. Run this from the project root.",
       call. = FALSE)

if (nzchar(Sys.getenv("PPR_BUDGET"))) {
  message("Note: PPR_BUDGET was set in this shell (", Sys.getenv("PPR_BUDGET"),
          "); ignoring it. Budgets come from this script.")
  Sys.unsetenv("PPR_BUDGET")
}

# Resolve the run configuration for a budget (NULL = canonical) by sourcing
# 00_config.R, so run IDs and paths come from one place only.
config_for <- function(budget = NULL) {
  if (is.null(budget)) Sys.unsetenv("PPR_BUDGET")
  else Sys.setenv(PPR_BUDGET = format(budget))
  e <- new.env()
  sys.source("00_config.R", envir = e)
  e
}

canonical  <- config_for(NULL)
panel_path <- file.path(canonical$DIR_DERIVED, "data_panel.rds")
if (!file.exists(panel_path))
  stop("Derived data not found at ", panel_path, ".\n",
       "  Run the land-cover pipeline first (PPR_BUDGET unset): Rscript 00_run_all.R",
       call. = FALSE)

if (!DRY_RUN && !requireNamespace("gurobi", quietly = TRUE))
  stop("The 'gurobi' R package is not available; script 08 needs it.",
       call. = FALSE)

# ── Sweep ─────────────────────────────────────────────────────────────────────
rscript <- file.path(R.home("bin"), "Rscript")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
t_sweep <- Sys.time()
log_rows <- list()

message("Budget sweep | land cover: ", canonical$LC_LABEL,
        " | derived data: ", canonical$DIR_DERIVED,
        "\nBudgets (median EAUs/period): ", paste(BUDGETS, collapse = ", "),
        if (DRY_RUN) "   [DRY RUN]" else "",
        if (FORCE) "   [FORCE: rerun complete budgets]" else "")

for (b in BUDGETS) {
  cfg       <- config_for(b)                      # also sets PPR_BUDGET for the child
  res_file  <- file.path(cfg$DIR_OUT,  "model_results.csv")
  meta_file <- file.path(cfg$DIR_FIGS, "run_metadata.json")
  log_file  <- file.path(LOG_DIR, paste0(cfg$RUN_ID, ".log"))

  # 08 + 10 if results are missing (or --force); 10 only if just the assets are
  # missing; nothing if both exist.
  from <- if (FORCE || !file.exists(res_file)) "08"
          else if (!file.exists(meta_file))    "10"
          else NA_character_

  row <- data.frame(budget = b, run_id = cfg$RUN_ID, action = NA_character_,
                    minutes = NA_real_, log = log_file, stringsAsFactors = FALSE)

  if (is.na(from)) {
    message("\n[", cfg$RUN_ID, "] complete; skipping (use --force to rerun).")
    row$action <- "skipped (complete)"
  } else if (DRY_RUN) {
    message("\n[", cfg$RUN_ID, "] would run:  PPR_BUDGET=", format(b),
            " Rscript 00_run_all.R ", from, "\n    -> ", cfg$DIR_RUN,
            "\n    log: ", log_file)
    row$action <- paste("dry run: from", from)
  } else {
    message("\n[", cfg$RUN_ID, "] running from script ", from, "  (started ",
            format(Sys.time(), "%H:%M"), ")\n    log: ", log_file)
    t0 <- Sys.time()
    status <- system2(rscript, c("--no-save", "--no-restore", "00_run_all.R", from),
                      stdout = log_file, stderr = log_file)
    row$minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
    row$action  <- if (status == 0) paste("ran from", from)
                   else paste0("FAILED (exit ", status, ")")
    message("    ", row$action, " in ", row$minutes, " min")
    if (status != 0) message("    See the end of ", log_file, " for the error.")
  }
  log_rows[[length(log_rows) + 1]] <- row
}
Sys.unsetenv("PPR_BUDGET")

sweep_log <- do.call(rbind, log_rows)
message("\n══ Sweep summary (", round(as.numeric(difftime(Sys.time(), t_sweep,
        units = "mins")), 1), " min) ══")
print(sweep_log[, c("budget", "run_id", "action", "minutes")], row.names = FALSE)

if (DRY_RUN) quit(save = "no", status = 0)

# ── Policy ranking per run (climate scenarios) ────────────────────────────────
# Same check as 08 prints: rolling >= myopic >= greedy on total welfare J.
message("\n══ Policy ranking (rolling >= myopic >= greedy, climate scenarios) ══")
for (b in BUDGETS) {
  cfg <- config_for(b)
  f   <- file.path(cfg$DIR_OUT, "model_results.csv")
  if (!file.exists(f)) { message("  ", cfg$RUN_ID, ": no results"); next }
  r <- read.csv(f, stringsAsFactors = FALSE)
  r <- r[r$scenario != "stationary", ]
  w <- reshape(r[, c("scenario", "model", "J")], idvar = "scenario",
               timevar = "model", direction = "wide")
  ok <- w$J.rolling >= w$J.myopic - 1e-6 & w$J.myopic >= w$J.greedy - 1e-6
  message(sprintf("  %-14s %s", cfg$RUN_ID,
                  if (all(ok)) "holds in all climate scenarios"
                  else paste("FAILS in:", paste(w$scenario[!ok], collapse = ", "))))
}
Sys.unsetenv("PPR_BUDGET")

# ── Budget 5 vs canonical ─────────────────────────────────────────────────────
if (5 %in% BUDGETS) {
  b5         <- config_for(5); Sys.unsetenv("PPR_BUDGET")
  canon_file <- file.path(canonical$DIR_OUT, "model_results.csv")
  b5_file    <- file.path(b5$DIR_OUT,        "model_results.csv")
  if (file.exists(canon_file) && file.exists(b5_file)) {
    keep <- c("scenario", "model", "value_added", "n_acquired")
    m <- merge(read.csv(canon_file)[, keep], read.csv(b5_file)[, keep],
               by = c("scenario", "model"), suffixes = c("_canon", "_b05"))
    m$rel_diff_pct <- round(100 * (m$value_added_b05 - m$value_added_canon) /
                              m$value_added_canon, 3)
    m$flag <- ifelse(m$scenario != "stationary" &
                       abs(m$rel_diff_pct) > 100 * CANON_TOL, "CHECK", "")
    m <- data.frame(scenario = m$scenario, model = m$model,
                    va_canon = round(m$value_added_canon),
                    va_b05   = round(m$value_added_b05),
                    diff_pct = m$rel_diff_pct,
                    n_canon  = m$n_acquired_canon,
                    n_b05    = m$n_acquired_b05,
                    flag     = m$flag)
    message("\n══ ", b5$RUN_ID, " vs canonical ", canonical$RUN_ID,
            " (value_added; flag if |diff| > ", 100 * CANON_TOL,
            "% on climate scenarios) ══")
    print(m[order(m$scenario, m$model), ], row.names = FALSE)
    if (any(m$flag == "CHECK"))
      message("\nSome differences exceed the tolerance. Compare the solver ",
              "convergence tables (solver_convergence.csv) of the two runs.")
  } else {
    message("\nb05 vs canonical check skipped: results file(s) missing.")
  }
}

if (any(grepl("^FAILED", sweep_log$action))) quit(save = "no", status = 1)
