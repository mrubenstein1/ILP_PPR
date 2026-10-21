# ═══════════════════════════════════════════════════════════════════════════════
# 10b_results_assets_landcover.R — OPTIONAL: assets comparing two land-cover runs
# ═══════════════════════════════════════════════════════════════════════════════
# Compares this run (the BASE: RUN_ID in 00_config.R, normally "gcam_ref") with an
# ALTERNATIVE land-cover run (COMPARE_RUN below, normally "bau"). Differences are
# reported as alternative − base, with the same formulas as 12_compare_runs.qmd.
#
# PREREQUISITES — both complete pipeline runs must exist on disk:
#   runs/<RUN_ID>/output_data, runs/<RUN_ID>/derived_data       (base)
#   runs/<COMPARE_RUN>/output_data, runs/<COMPARE_RUN>/derived_data  (alternative)
# Typical order: run the pipeline for the alternative scenario, then for the base
# scenario, then this script (see README).
#
# Writes to the BASE run's output_figs (raw numbers; formatting lives in
# 11b_landcover_comparison.qmd):
#   cmp_value_added.csv   acquisition value by policy, both runs, all scenarios
#   cmp_headline.csv      % improvements for the climate scenarios + ranking check
#   cmp_risk_table.csv    per-period risk summary by RCP x decade, both runs
#   fig_cmp_risk_box      per-period risk, both runs side by side (PNG + PDF)
#   cmp_metadata.json     run ids, labels, short names, threshold, epsilons
# 11_results.qmd includes the comparison section only if cmp_metadata.json exists.
# ═══════════════════════════════════════════════════════════════════════════════

if (!exists(".SETUP_DONE")) source("00_setup.R")

# ── CONFIG ─────────────────────────────────────────────────────────────────────
COMPARE_RUN   <- "bau"                               # RUN_ID of the alternative run
CMP_SHORT     <- c(base = "GCAM ref", alt = "BAU")   # short names: table headers, legends
CMP_RUN_COL   <- c(base = "#4D4D4D", alt = "#C0392B")
CMP_THRESH_PP <- 2       # flag a headline gap that moves by more than this (pp)

OUT_DIR     <- DIR_FIGS                                   # the base run's output_figs
BASE_RES    <- file.path(DIR_OUT, "model_results.csv")
BASE_PANEL  <- file.path(DIR_DERIVED, "data_panel.rds")
ALT_RES     <- file.path("runs", COMPARE_RUN, "output_data", "model_results.csv")
ALT_PANEL   <- file.path("runs", COMPARE_RUN, "derived_data", "data_panel.rds")
ALT_META    <- file.path("runs", COMPARE_RUN, "output_figs", "run_metadata.json")

DPI         <- 300
BASE_SIZE   <- 12
BASE_FAMILY <- "serif"          # same as 10_results_assets.R
SAVE_PDF    <- TRUE
# ────────────────────────────────────────────────────────────────────────────────

if (COMPARE_RUN == RUN_ID)
  stop("COMPARE_RUN ('", COMPARE_RUN, "') is the same as this run. Set 00_config.R to ",
       "the base run (normally gcam_ref) before running this script.", call. = FALSE)
for (f in c(BASE_RES, ALT_RES))
  if (!file.exists(f))
    stop("Missing '", f, "'. Both pipeline runs must be complete before running 10b ",
         "(see the header of this script).", call. = FALSE)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

## packages ---------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(scales)
})

## labels, helpers (same definitions as 10_results_assets.R) --------------------
SCEN_LEVELS <- c("rcp45_wet", "rcp45_dry", "rcp85_wet", "rcp85_dry", "stationary")
CLIM <- SCEN_LEVELS[1:4]
pct_improve <- function(roll, base) ifelse(abs(base) < 1e-9, NA_real_, 100 * (roll - base) / base)

theme_ppr <- function(base_size = BASE_SIZE, base_family = BASE_FAMILY) {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.title    = element_text(colour = "grey25"),
      axis.title.y  = element_text(angle = 90, margin = margin(r = 6)),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "top", legend.justification = "right",
      legend.title = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0, size = rel(0.95)),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(10, 14, 8, 10)
    )
}
theme_set(theme_ppr())

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = DPI, bg = "white")
  if (SAVE_PDF) {
    ok <- tryCatch(
      { ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h,
               device = cairo_pdf, bg = "white"); TRUE },
      warning = function(w) FALSE,
      error   = function(e) FALSE)
    if (!ok)
      suppressWarnings(ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p,
                              width = w, height = h, bg = "white"))
  }
  invisible(p)
}

# ═══════════════════════════════════════════════════════════════════════════════
# MODEL RESULTS  —  improvement of rolling over greedy and myopic, both runs
# ═══════════════════════════════════════════════════════════════════════════════
summ_cmp <- function(r) {
  r %>% mutate(scenario = as.character(scenario)) %>%
    select(scenario, model, value_added) %>%
    pivot_wider(names_from = model, values_from = value_added) %>%
    transmute(scenario,
              rg_dpairs = rolling - greedy,
              rm_dpairs = rolling - myopic,
              # "<1 pair" differences (stationary) are not meaningful as %
              rg_pct = ifelse(abs(rg_dpairs) < 1, NA_real_, pct_improve(rolling, greedy)),
              rm_pct = ifelse(abs(rm_dpairs) < 1, NA_real_, pct_improve(rolling, myopic)),
              order_ok = rolling >= myopic - 1e-6 & myopic >= greedy - 1e-6)
}
cmp <- inner_join(summ_cmp(read_csv(BASE_RES, show_col_types = FALSE)),
                  summ_cmp(read_csv(ALT_RES,  show_col_types = FALSE)),
                  by = "scenario", suffix = c("_base", "_alt")) %>%
  mutate(scenario = factor(scenario, levels = SCEN_LEVELS)) %>%
  arrange(scenario)

tbl_cmp_improvement <- cmp %>%
  transmute(scenario,
            rg_dpairs_base, rg_dpairs_alt, rg_pct_base, rg_pct_alt,
            rg_pct_diff = rg_pct_alt - rg_pct_base,
            rm_dpairs_base, rm_dpairs_alt, rm_pct_base, rm_pct_alt,
            rm_pct_diff = rm_pct_alt - rm_pct_base)

## acquisition value (value added over the do-nothing baseline) by policy, both runs
va_of <- function(f) read_csv(f, show_col_types = FALSE) %>%
  mutate(scenario = as.character(scenario)) %>%
  select(scenario, model, value_added) %>%
  pivot_wider(names_from = model, values_from = value_added)
rel_chg <- function(a, b) ifelse(abs(b) < 1, NA_real_, 100 * (a - b) / b)  # % change, alt vs base
tbl_cmp_value <- inner_join(va_of(BASE_RES), va_of(ALT_RES), by = "scenario",
                            suffix = c("_base", "_alt")) %>%
  mutate(scenario = factor(scenario, levels = SCEN_LEVELS)) %>%
  arrange(scenario) %>%
  transmute(scenario,
            greedy_base,  greedy_alt,  greedy_pct_chg  = rel_chg(greedy_alt,  greedy_base),
            myopic_base,  myopic_alt,  myopic_pct_chg  = rel_chg(myopic_alt,  myopic_base),
            rolling_base, rolling_alt, rolling_pct_chg = rel_chg(rolling_alt, rolling_base))
write_csv(tbl_cmp_value, file.path(OUT_DIR, "cmp_value_added.csv"))

tbl_cmp_headline <- tbl_cmp_improvement %>%
  filter(scenario %in% CLIM) %>%
  left_join(cmp %>% select(scenario, order_ok_base, order_ok_alt), by = "scenario") %>%
  transmute(scenario, rg_pct_base, rg_pct_alt, rg_pct_diff,
            rm_pct_base, rm_pct_alt, rm_pct_diff,
            ranking_holds = order_ok_base & order_ok_alt)
write_csv(tbl_cmp_headline, file.path(OUT_DIR, "cmp_headline.csv"))

# ═══════════════════════════════════════════════════════════════════════════════
# CONVERSION RISK  —  as Fig. 3 of the main results (RCP 4.5/8.5, 2100 excluded)
# ═══════════════════════════════════════════════════════════════════════════════
risk_of <- function(pnl) {
  sv  <- pnl %>% filter(rcp == "stationary", year < 2100) %>% pull(trans_prob)
  eps <- min(sv[sv > 0])           # smallest POSITIVE value (structural zeros excluded)
  r <- pnl %>% distinct(eau_id, year, rcp, value = trans_prob) %>%
    filter(year < 2100, rcp %in% c("45", "85")) %>%
    mutate(rcp_lab = factor(recode(rcp, `45` = "RCP 4.5", `85` = "RCP 8.5"),
                            levels = c("RCP 4.5", "RCP 8.5")),
           at_floor = value <= eps * 1.000001,
           pct      = value * 100)
  list(r = r, eps = eps)
}

eps_base <- NA_real_; eps_alt <- NA_real_; drew_risk <- FALSE
if (file.exists(BASE_PANEL) && file.exists(ALT_PANEL)) {
  rb <- risk_of(readRDS(BASE_PANEL)); ra <- risk_of(readRDS(ALT_PANEL))
  eps_base <- rb$eps; eps_alt <- ra$eps
  rk <- bind_rows(rb$r %>% mutate(run = "base"), ra$r %>% mutate(run = "alt")) %>%
    mutate(run = factor(run, levels = c("base", "alt")))

  tbl_cmp_risk <- rk %>%
    group_by(rcp_lab, year, run) %>%
    summarise(mean = mean(pct), median = median(pct), max = max(pct),
              pct_floor = 100 * mean(at_floor), .groups = "drop") %>%
    pivot_wider(names_from = run, values_from = c(mean, median, max, pct_floor)) %>%
    arrange(rcp_lab, year) %>%
    mutate(rcp_lab = as.character(rcp_lab))
  write_csv(tbl_cmp_risk, file.path(OUT_DIR, "cmp_risk_table.csv"))

  ## figure: boxes dodged by run; diamond = mean, dash = maximum
  dg <- position_dodge(width = .75)
  fig_cmp_risk_box <- ggplot(rk, aes(factor(year), pct, fill = run, colour = run)) +
    geom_boxplot(aes(group = interaction(year, run)), outlier.shape = NA,
                 position = dg, width = .62, alpha = .35, linewidth = .3) +
    stat_summary(aes(group = run), fun = max, geom = "point", shape = 95, size = 6,
                 position = dg, show.legend = FALSE) +
    stat_summary(aes(group = run), fun = mean, geom = "point", shape = 23, size = 1.6,
                 fill = "white", position = dg, show.legend = FALSE) +
    facet_wrap(~ rcp_lab, ncol = 2) +
    scale_fill_manual(values = CMP_RUN_COL, labels = CMP_SHORT) +
    scale_colour_manual(values = CMP_RUN_COL, labels = CMP_SHORT) +
    scale_y_continuous("Per-period conversion risk",
                       labels = label_number(suffix = "%", accuracy = .1),
                       expand = expansion(mult = c(.02, .06))) +
    scale_x_discrete("year") +
    labs(caption = "Box = median and interquartile range; diamond = mean; dash = maximum.") +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = .5, size = rel(.8)),
          plot.caption = element_text(colour = "grey40", size = rel(.72), hjust = 0))
  save_fig(fig_cmp_risk_box, "fig_cmp_risk_box", 11, 5.0)
  drew_risk <- TRUE
} else {
  message("Risk figure/table skipped: inputs panel missing ('", BASE_PANEL,
          "' or '", ALT_PANEL, "'). Rerun that pipeline through script 06.")
}

# ═══════════════════════════════════════════════════════════════════════════════
# METADATA  —  11_results.qmd includes the comparison section only if this exists
# ═══════════════════════════════════════════════════════════════════════════════
alt_meta <- tryCatch(jsonlite::read_json(ALT_META, simplifyVector = TRUE),
                     error = function(e) NULL)
jsonlite::write_json(
  list(generated_date = format(Sys.Date(), "%Y-%m-%d"),
       base_run = RUN_ID, base_label = LC_LABEL, base_short = CMP_SHORT[["base"]],
       alt_run = COMPARE_RUN,
       alt_label = if (!is.null(alt_meta$lc_label)) alt_meta$lc_label else COMPARE_RUN,
       alt_short = CMP_SHORT[["alt"]],
       alt_generated_date = if (!is.null(alt_meta$generated_date)) alt_meta$generated_date else NA,
       thresh_pp = CMP_THRESH_PP, eps_base = eps_base, eps_alt = eps_alt,
       risk_drawn = drew_risk),
  file.path(OUT_DIR, "cmp_metadata.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")

message("\nDone. Comparison of '", RUN_ID, "' (base) with '", COMPARE_RUN, "' written to '",
        OUT_DIR, "/': 2 tables", if (drew_risk) " + risk table + risk figure" else "",
        " + cmp_metadata.json.")
