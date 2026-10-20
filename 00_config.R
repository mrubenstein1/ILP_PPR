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

# Output folder name. Defaults to the scenario prefix; change only if you want
# more than one run of the same land cover (e.g. a sensitivity analysis).
RUN_ID      <- LC_SCENARIO

# ── Folder layout ─────────────────────────────────────────────────────────────
DIR_IN      <- "input_data"                          # external, read-only
DIR_COMMON  <- file.path("runs", "common")           # EAU grid + lookup (01), shared
DIR_RUN     <- file.path("runs", RUN_ID)
DIR_DERIVED <- file.path(DIR_RUN, "derived_data")    # pipeline intermediates (02-06)
DIR_OUT     <- file.path(DIR_RUN, "output_data")     # model results (08)
DIR_FIGS    <- file.path(DIR_RUN, "output_figs")     # tables, figures, maps (10)

# Land-cover input folder for one variant: "rcp45", "rcp85" or "meanclim".
LC_PREFIX <- paste0("prairie_potholes_", LC_SCENARIO)
lc_dir    <- function(variant) file.path(DIR_IN, paste0(LC_PREFIX, "_", variant))
