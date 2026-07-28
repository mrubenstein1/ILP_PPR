# ══════════════════════════════════════════════════════════════════════════════
# 00_setup.R — packages, paths, and directory scaffolding
# Sourced by every other script. Contains no analysis.
# ══════════════════════════════════════════════════════════════════════════════

# ── Packages ──────────────────────────────────────────────────────────────────
.core <- c("dplyr", "tidyr", "readr", "stringr", "purrr", "ggplot2",
           "scales", "Matrix", "jsonlite")
.spatial <- c("terra", "sf")          # 01, 02, 06 require; 10 optional
.report  <- c("flextable", "knitr")   # 11 only

.check <- function(pkgs, label, hard = TRUE) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    msg <- sprintf("Missing %s package(s): %s\n  install.packages(c(%s))",
                   label, paste(miss, collapse = ", "),
                   paste(sprintf('"%s"', miss), collapse = ", "))
    if (hard) stop(msg, call. = FALSE) else message(msg)
    return(FALSE)
  }
  TRUE
}
.check(.core, "core")
HAS_SPATIAL <- .check(.spatial, "spatial", hard = FALSE)

.check(.report, "reporting", hard = FALSE)

suppressPackageStartupMessages(
  invisible(lapply(.core, library, character.only = TRUE))
)

if (HAS_SPATIAL) suppressPackageStartupMessages(
  invisible(lapply(.spatial, library, character.only = TRUE))
)

# ── Paths ─────────────────────────────────────────────────────────────────────
DIR_IN      <- "input_data"     # external, read-only — never written by scripts
DIR_DERIVED <- "derived_data"   # pipeline intermediates (01–06)
DIR_OUT     <- "output_data"    # model results (08)
DIR_FIGS    <- "output_figs"    # tables, figures, maps (10)

for (d in c(DIR_DERIVED, DIR_OUT, DIR_FIGS))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

.SETUP_DONE <- TRUE
message("Setup complete | R ", getRversion(),
        " | spatial stack: ", if (HAS_SPATIAL) "available" else "NOT available")