# ══════════════════════════════════════════════════════════════════════════════
# 00_run_all.R — run the pipeline, each script in its own fresh R session
# ══════════════════════════════════════════════════════════════════════════════
#
# USAGE (from the project root):
#
#   Rscript 00_run_all.R           run the whole pipeline
#   Rscript 00_run_all.R 05        run 05 through the end
#   Rscript 00_run_all.R 05 08     run 05 through 08
#
# Each script is launched as a separate process, so nothing can accidentally
# depend on objects left in memory by an earlier script. That is what makes a
# successful run evidence of reproducibility.
#
# After this finishes, render 11_results.qmd to build the results document.
# ══════════════════════════════════════════════════════════════════════════════

RUN_VALIDATION <- FALSE   # TRUE to also run 09 (slow: re-solves every scenario)

all_scripts <- c(
  "01_create_EAUs.R",
  "02_EAU_prop_suitable.R",
  "03_create_data_panel.R",
  "04_benefit_data.R",
  "05_risk_of_loss.R",
  "06_cost_data.R",
  if (RUN_VALIDATION) "09_validation.R",
  "08_run_models.R",
  "10_results_assets.R"
)

# ── Optional range from the command line ──────────────────────────────────────
args    <- commandArgs(trailingOnly = TRUE)
scripts <- all_scripts

if (length(args) >= 1) {
  from <- grep(paste0("^", args[1]), all_scripts)
  to   <- if (length(args) >= 2) grep(paste0("^", args[2]), all_scripts)
          else length(all_scripts)
  if (!length(from) || !length(to))
    stop("Script number not recognised. Use e.g. 'Rscript 00_run_all.R 05 08'.",
         call. = FALSE)
  scripts <- all_scripts[from[1]:to[1]]
}

# ── Sanity check: are we in the project root? ─────────────────────────────────
if (!file.exists("00_setup.R"))
  stop("00_setup.R not found. Run this from the project root directory.",
       call. = FALSE)

# ── Show which run this is, so a wrong configuration is caught before a long run ─
source("00_config.R")
message("Run: ", RUN_ID, " | land cover: ", LC_LABEL, " (", LC_PREFIX, ")",
        " | budget: ", BUDGET_EAUS_PER_PERIOD, " median EAUs/period",
        "\nDerived data from: ", DIR_DERIVED,
        "\nOutputs go to: ", DIR_RUN, "\n")

# A budget run reuses another run's derived data; it must not rerun 01-06.
upstream_scripts <- c("01_create_EAUs.R", "02_EAU_prop_suitable.R",
                      "03_create_data_panel.R", "04_benefit_data.R",
                      "05_risk_of_loss.R", "06_cost_data.R")
if (RUN_ID != DERIVED_RUN_ID && any(scripts %in% upstream_scripts))
  stop("Budget run '", RUN_ID, "' reuses derived data from '", DERIVED_RUN_ID,
       "'. Run only from 08 onward, e.g.  PPR_BUDGET=",
       BUDGET_EAUS_PER_PERIOD, " Rscript 00_run_all.R 08", call. = FALSE)

# A high-hazard run has its own derived_data/, so the check above does not catch
# it. Its panel is a rescaled copy of the canonical one (14_build_hazard_panel.R);
# rerunning 01-06 here would reprocess land cover into the wrong folder.
if (HAZARD_RUN && any(scripts %in% upstream_scripts))
  stop("High-hazard run '", RUN_ID, "' uses a rescaled copy of the canonical panel ",
       "and must not rerun 01-06. Build the panel, then run from 08:\n",
       "  PPR_HAZARD=", .hazard_env, " Rscript 14_build_hazard_panel.R\n",
       "  PPR_HAZARD=", .hazard_env, " Rscript 00_run_all.R 08", call. = FALSE)

missing <- scripts[!file.exists(scripts)]
if (length(missing))
  stop("Missing script(s): ", paste(missing, collapse = ", "), call. = FALSE)

# ── Run ───────────────────────────────────────────────────────────────────────
rscript <- file.path(R.home("bin"), "Rscript")
t_start <- Sys.time()

message("Pipeline: ", paste(scripts, collapse = " -> "), "\n")

for (s in scripts) {
  t_script <- Sys.time()
  message("\n================================================================")
  message("RUNNING: ", s, "   (started ", format(Sys.time(), "%H:%M:%S"), ")")
  message("================================================================")

  status <- system2(rscript, c("--no-save", "--no-restore", shQuote(s)))

  if (status != 0)
    stop("FAILED at ", s, " (exit status ", status, "). Pipeline halted.",
         call. = FALSE)

  message(sprintf("--- %s finished in %.1f min ---", s,
                  as.numeric(difftime(Sys.time(), t_script, units = "mins"))))
}

message(sprintf("\nPipeline complete in %.1f min. Now render 11_results.qmd.",
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
