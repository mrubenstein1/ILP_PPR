# ══════════════════════════════════════════════════════════════════════════════
# 15_run_hazard_sweep.R — conversion-hazard sensitivity sweep
# ══════════════════════════════════════════════════════════════════════════════
#
# For each target ANNUAL conversion rate in RATES below, builds the rescaled panel
# (14_build_hazard_panel.R) and runs the models and per-run assets (08 then 10) in
# its own folder, runs/<LC_SCENARIO>_haz<tag>/. Land-cover processing (01-06) is
# NOT rerun. The canonical run (k = 1) is the sweep's lowest point and is not
# rerun here.
#
# PREREQUISITE: the canonical run must exist (Rscript 00_run_all.R, with
# PPR_HAZARD and PPR_BUDGET unset).
#
# USAGE (from the project root):
#   Rscript 15_run_hazard_sweep.R               run the rates in RATES below
#   Rscript 15_run_hazard_sweep.R 0.09 0.63     run only these rates (% per year)
#   Rscript 15_run_hazard_sweep.R --dry-run     show what would run; run nothing
#   Rscript 15_run_hazard_sweep.R --force       rebuild and rerun even if complete
#
# On a Mac, stop the machine sleeping during the multi-hour sweep:
#   caffeinate -i Rscript 15_run_hazard_sweep.R
#
# Each rate runs exactly what you would type by hand:
#   PPR_HAZARD=<rate> Rscript 14_build_hazard_panel.R
#   PPR_HAZARD=<rate> Rscript 00_run_all.R 08
# Finished rates are skipped, so an interrupted sweep restarts with the same
# command. A rate whose scaled panel was built from an older canonical panel
# (source md5 differs) is rebuilt and rerun automatically.
#
# Console output goes to logs/hazard_sweep/<RUN_ID>_14.log (panel build) and
# logs/hazard_sweep/<RUN_ID>.log (models). Follow a run live with:
#   tail -f logs/hazard_sweep/<RUN_ID>.log
#
# Then, with PPR_HAZARD unset:  Rscript 10e_results_assets_hazard_sweep.R
# Base R only: no packages needed to run this driver.
# ══════════════════════════════════════════════════════════════════════════════

# ── Settings ──────────────────────────────────────────────────────────────────
# Target annual conversion rates (% per year), with the published estimates they
# stand for. The canonical run (~0.17%/yr, habitat-weighted 2020 rate) is the
# lowest point and needs no run. Edit freely; 10e reads the rates from the runs.
RATES <- c(
  0.28,   # Johnston 2013: wetland loss to row crops, Dakota PPR (NWI basis)
  0.40,   # Stephens et al. 2008: native grassland conversion, Missouri Coteau
  0.57,   # Kemink et al. 2023: upper bound of wetland drainage
  0.95,   # Rashford et al. 2011: grassland conversion, 2006-2011
  1.33    # Rashford et al. 2011: grassland conversion (high-hazard run)
)

LOG_DIR <- file.path("logs", "hazard_sweep")

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
  RATES <- suppressWarnings(as.numeric(nums))
  if (anyNA(RATES) || any(RATES <= 0 | RATES >= 100))
    stop("Rates are annual % in (0, 100), e.g. ",
         "'Rscript 15_run_hazard_sweep.R 0.4 0.95'.", call. = FALSE)
}

# ── Checks ────────────────────────────────────────────────────────────────────
if (!file.exists("00_run_all.R") || !file.exists("14_build_hazard_panel.R"))
  stop("00_run_all.R / 14_build_hazard_panel.R not found. Run this from the ",
       "project root.", call. = FALSE)

for (v in c("PPR_HAZARD", "PPR_BUDGET"))
  if (nzchar(Sys.getenv(v))) {
    message("Note: ", v, " was set in this shell (", Sys.getenv(v),
            "); ignoring it. Rates come from this script.")
    Sys.unsetenv(v)
  }

# Resolve the run configuration for a rate (NULL = canonical) by sourcing
# 00_config.R, so run IDs and paths come from one place only.
config_for <- function(rate = NULL) {
  Sys.unsetenv("PPR_BUDGET")
  if (is.null(rate)) Sys.unsetenv("PPR_HAZARD")
  else Sys.setenv(PPR_HAZARD = format(rate))
  e <- new.env()
  sys.source("00_config.R", envir = e)
  e
}

canonical  <- config_for(NULL)
panel_path <- file.path(canonical$DIR_DERIVED, "data_panel.rds")
if (!file.exists(panel_path))
  stop("Canonical derived data not found at ", panel_path, ".\n",
       "  Run the canonical pipeline first (PPR_HAZARD unset): Rscript 00_run_all.R",
       call. = FALSE)
canon_md5 <- unname(tools::md5sum(panel_path))

if (!DRY_RUN && !requireNamespace("gurobi", quietly = TRUE))
  stop("The 'gurobi' R package is not available; script 08 needs it.",
       call. = FALSE)

# md5 of the canonical panel a scaled panel was built from (NA if unknown).
built_from_md5 <- function(json) {
  if (!file.exists(json)) return(NA_character_)
  txt <- paste(readLines(json, warn = FALSE), collapse = "")
  m   <- regmatches(txt, regexpr('"source_panel_md5"\\s*:\\s*"[0-9a-f]+"', txt))
  if (!length(m)) NA_character_ else sub('.*"([0-9a-f]+)"$', "\\1", m)
}

# ── Sweep ─────────────────────────────────────────────────────────────────────
rscript <- file.path(R.home("bin"), "Rscript")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
t_sweep  <- Sys.time()
log_rows <- list()

message("Hazard sweep | land cover: ", canonical$LC_LABEL,
        " | canonical derived data: ", canonical$DIR_DERIVED,
        "\nRates (% per year): ", paste(RATES, collapse = ", "),
        if (DRY_RUN) "   [DRY RUN]" else "",
        if (FORCE) "   [FORCE: rebuild and rerun complete rates]" else "")

for (r in RATES) {
  cfg       <- config_for(r)                     # also sets PPR_HAZARD for the children
  json_file <- file.path(cfg$DIR_DERIVED, "hazard_scaling.json")
  res_file  <- file.path(cfg$DIR_OUT,  "model_results.csv")
  meta_file <- file.path(cfg$DIR_FIGS, "run_metadata.json")
  log14     <- file.path(LOG_DIR, paste0(cfg$RUN_ID, "_14.log"))
  log_file  <- file.path(LOG_DIR, paste0(cfg$RUN_ID, ".log"))

  # Rebuild the panel if forced, missing, or built from an older canonical panel;
  # a rebuilt panel always means rerunning the models.
  stale <- !is.na(built_from_md5(json_file)) && built_from_md5(json_file) != canon_md5
  build <- FORCE || !file.exists(json_file) || stale
  from  <- if (build || !file.exists(res_file)) "08"
           else if (!file.exists(meta_file))    "10"
           else NA_character_

  row <- data.frame(rate = r, run_id = cfg$RUN_ID, action = NA_character_,
                    minutes = NA_real_, log = log_file, stringsAsFactors = FALSE)

  if (is.na(from)) {
    message("\n[", cfg$RUN_ID, "] complete; skipping (use --force to rerun).")
    row$action <- "skipped (complete)"
  } else if (DRY_RUN) {
    message("\n[", cfg$RUN_ID, "] would run:",
            if (build) paste0("  PPR_HAZARD=", format(r), " Rscript 14_build_hazard_panel.R",
                              if (stale) "   (stale: canonical panel changed)" else "") else "",
            "\n    PPR_HAZARD=", format(r), " Rscript 00_run_all.R ", from,
            "\n    -> ", cfg$DIR_RUN)
    row$action <- paste0("dry run: ", if (build) "build + " else "", "from ", from)
  } else {
    message("\n[", cfg$RUN_ID, "] ", if (build) "building panel, then " else "",
            "running from script ", from, "  (started ", format(Sys.time(), "%H:%M"),
            ")\n    log: ", log_file)
    t0 <- Sys.time()
    status <- 0
    if (build)
      status <- system2(rscript, c("--no-save", "--no-restore", "14_build_hazard_panel.R"),
                        stdout = log14, stderr = log14)
    if (status == 0)
      status <- system2(rscript, c("--no-save", "--no-restore", "00_run_all.R", from),
                        stdout = log_file, stderr = log_file)
    row$minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
    row$action  <- if (status == 0) paste0(if (build) "built + " else "", "ran from ", from)
                   else paste0("FAILED (exit ", status, ")")
    message("    ", row$action, " in ", row$minutes, " min")
    if (status != 0) message("    See the end of ", log14, " and ", log_file, " for the error.")
  }
  log_rows[[length(log_rows) + 1]] <- row
}
Sys.unsetenv("PPR_HAZARD")

sweep_log <- do.call(rbind, log_rows)
message("\n══ Sweep summary (", round(as.numeric(difftime(Sys.time(), t_sweep,
        units = "mins")), 1), " min) ══")
print(sweep_log[, c("rate", "run_id", "action", "minutes")], row.names = FALSE)

if (DRY_RUN) quit(save = "no", status = 0)

# ── Policy ranking per run (climate scenarios) ────────────────────────────────
# Same check as 08 prints: rolling >= myopic >= greedy on total welfare J.
message("\n══ Policy ranking (rolling >= myopic >= greedy, climate scenarios) ══")
for (r in RATES) {
  cfg <- config_for(r)
  f   <- file.path(cfg$DIR_OUT, "model_results.csv")
  if (!file.exists(f)) { message("  ", cfg$RUN_ID, ": no results"); next }
  res <- read.csv(f, stringsAsFactors = FALSE)
  res <- res[res$scenario != "stationary", ]
  w <- reshape(res[, c("scenario", "model", "J")], idvar = "scenario",
               timevar = "model", direction = "wide")
  ok <- w$J.rolling >= w$J.myopic - 1e-6 & w$J.myopic >= w$J.greedy - 1e-6
  message(sprintf("  %-18s %s", cfg$RUN_ID,
                  if (all(ok)) "holds in all climate scenarios"
                  else paste("FAILS in:", paste(w$scenario[!ok], collapse = ", "))))
}
Sys.unsetenv("PPR_HAZARD")

message("\nNext (PPR_HAZARD unset):  Rscript 10e_results_assets_hazard_sweep.R")
if (any(grepl("^FAILED", sweep_log$action))) quit(save = "no", status = 1)
