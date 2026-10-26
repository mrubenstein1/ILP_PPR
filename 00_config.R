# ══════════════════════════════════════════════════════════════════════════════
# 00_config.R — run configuration and folder layout
# Sourced by 00_setup.R (and so by every script), by 00_run_all.R, and by
# 11_results.qmd. Plain assignments only: no packages, nothing is created here.
#
# TO SWITCH LAND-COVER SCENARIO: edit LC_SCENARIO and LC_LABEL below.
# Each run writes to its own folder, runs/<RUN_ID>/, so runs never overwrite
# each other. The EAU grid (script 01) is shared by all runs in runs/common/.
# ══════════════════════════════════════════════════════════════════════════════

# ── Run settings ──────────────────────────────────────────────────────────────
# FOREsce scenario prefix. Input folders must be named
#   input_data/prairie_potholes_<LC_SCENARIO>_{rcp45,rcp85,meanclim}
# with files prairie_potholes_<LC_SCENARIO>_{rcp45,rcp85}_YYYY.tif and
# prairie_potholes_<LC_SCENARIO>_meanclim_2020.tif.
LC_SCENARIO <- "gcam_ref"

# Human-readable name, printed in the results document under Data Sources.
LC_LABEL    <- "GCAM reference scenario"

# ── Budget ────────────────────────────────────────────────────────────────────
# Per-period budget = BUDGET_EAUS_PER_PERIOD x median EAU cost in that period
# (make_budget() in 07_ilp_core.R). Canonical value: 5.
# For the budget sensitivity sweep, set the environment variable PPR_BUDGET when
# launching a script, e.g.   PPR_BUDGET=2 Rscript 00_run_all.R 08
# Leave PPR_BUDGET unset for the canonical run.
BUDGET_CANONICAL <- 5
.budget_env  <- Sys.getenv("PPR_BUDGET", unset = "")
BUDGET_SWEEP <- nzchar(.budget_env)          # TRUE when budget is set from outside
BUDGET_EAUS_PER_PERIOD <- if (BUDGET_SWEEP)
  suppressWarnings(as.numeric(.budget_env)) else BUDGET_CANONICAL
if (!is.finite(BUDGET_EAUS_PER_PERIOD) || BUDGET_EAUS_PER_PERIOD <= 0)
  stop("PPR_BUDGET must be a positive number; got '", .budget_env, "'.",
       call. = FALSE)

# ── Run identifiers ───────────────────────────────────────────────────────────
# DERIVED_RUN_ID: the run whose derived_data/ (scripts 02-06) is used. The budget
# does not enter scripts 02-06, so budget runs reuse the land-cover run's derived
# data instead of reprocessing land cover.
DERIVED_RUN_ID <- LC_SCENARIO

# RUN_ID: where model results (08) and assets (10) are written.
#   PPR_BUDGET unset -> runs/<LC_SCENARIO>/       canonical (unchanged)
#   PPR_BUDGET=2     -> runs/<LC_SCENARIO>_b02/   sweep run
# PPR_BUDGET=5 deliberately writes to _b05, not the canonical folder, so a sweep
# can never overwrite the canonical results.
.budget_tag <- if (BUDGET_EAUS_PER_PERIOD %% 1 == 0)
  sprintf("%02d", as.integer(BUDGET_EAUS_PER_PERIOD)) else
    sub(".", "p", format(BUDGET_EAUS_PER_PERIOD), fixed = TRUE)   # e.g. 2.5 -> 2p5
RUN_ID <- if (BUDGET_SWEEP) paste0(LC_SCENARIO, "_b", .budget_tag) else LC_SCENARIO

# ── Folder layout ─────────────────────────────────────────────────────────────
DIR_IN      <- "input_data"                                      # external, read-only
DIR_COMMON  <- file.path("runs", "common")                       # EAU grid + lookup (01), shared
DIR_RUN     <- file.path("runs", RUN_ID)
DIR_DERIVED <- file.path("runs", DERIVED_RUN_ID, "derived_data") # intermediates (02-06), may be shared
DIR_OUT     <- file.path(DIR_RUN, "output_data")                 # model results (08)
DIR_FIGS    <- file.path(DIR_RUN, "output_figs")                 # tables, figures, maps (10)

# Land-cover input folder for one variant: "rcp45", "rcp85" or "meanclim".
LC_PREFIX <- paste0("prairie_potholes_", LC_SCENARIO)
lc_dir    <- function(variant) file.path(DIR_IN, paste0(LC_PREFIX, "_", variant))