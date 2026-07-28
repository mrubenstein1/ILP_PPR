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
