# ═══════════════════════════════════════════════════════════════════════════════
# 10d_results_assets_hazard.R — OPTIONAL: assets comparing the canonical run with
#                                the high-hazard run
# ═══════════════════════════════════════════════════════════════════════════════
# Compares this run (the BASE: RUN_ID in 00_config.R, normally "gcam_ref", run with
# PPR_HAZARD and PPR_BUDGET UNSET) with the high-hazard run named by COMPARE_RUN.
# The two runs differ only in trans_prob, rescaled by 14_build_hazard_panel.R.
#
# PREREQUISITES — both runs complete:
#   Rscript 00_run_all.R                                   (canonical)
#   PPR_HAZARD=1.33 Rscript 14_build_hazard_panel.R        (scaled panel)
#   PPR_HAZARD=1.33 Rscript 00_run_all.R 08                (high-hazard models + assets)
# then, with PPR_HAZARD unset:
#   Rscript 10d_results_assets_hazard.R
#
# Writes to the BASE run's output_figs (raw numbers; formatting lives in
# 11d_hazard_sensitivity.qmd):
#   haz_headline.csv      % improvement of rolling, both runs, climate scenarios + ranking
#   haz_value_added.csv   acquisition value by policy, both runs, ratio high/canonical
#   haz_overlap.csv       overlap (Jaccard) of the acquired EAU sets, by scenario x policy
#   fig_haz_risk_box      per-period risk, both runs side by side (PNG + PDF)
#   haz_metadata.json     run ids, target rate, k, staleness check of the source panel
# 11_results.qmd includes the hazard section only if haz_metadata.json exists.
# ═══════════════════════════════════════════════════════════════════════════════

if (!exists(".SETUP_DONE")) source("00_setup.R")

# ── CONFIG ─────────────────────────────────────────────────────────────────────
COMPARE_RUN   <- paste0(LC_SCENARIO, "_haz1p33")          # RUN_ID of the high-hazard run
HAZ_SHORT     <- c(base = "Canonical", alt = "High hazard") # table headers, legends
HAZ_RUN_COL   <- c(base = "#4D4D4D", alt = "#C0392B")
HAZ_THRESH_PP <- 2       # flag a headline gap that moves by more than this (pp)

OUT_DIR     <- DIR_FIGS                                     # the base run's output_figs
BASE_RES    <- file.path(DIR_OUT, "model_results.csv")
BASE_SCHED  <- file.path(DIR_OUT, "acquisition_schedule_spatial.csv")
BASE_PANEL  <- file.path(DIR_DERIVED, "data_panel.rds")
ALT_RES     <- file.path("runs", COMPARE_RUN, "output_data", "model_results.csv")
ALT_SCHED   <- file.path("runs", COMPARE_RUN, "output_data", "acquisition_schedule_spatial.csv")
ALT_PANEL   <- file.path("runs", COMPARE_RUN, "derived_data", "data_panel.rds")
ALT_SCALING <- file.path("runs", COMPARE_RUN, "derived_data", "hazard_scaling.json")

DPI         <- 300
BASE_SIZE   <- 12
BASE_FAMILY <- "serif"          # same as 10_results_assets.R
SAVE_PDF    <- TRUE
# ────────────────────────────────────────────────────────────────────────────────

if (HAZARD_RUN || BUDGET_SWEEP)
  stop("Run 10d from the canonical configuration: unset PPR_HAZARD and PPR_BUDGET.",
       call. = FALSE)
for (f in c(BASE_RES, ALT_RES, ALT_SCALING))
  if (!file.exists(f))
    stop("Missing '", f, "'. Both runs must be complete before running 10d ",
         "(see the header of this script).", call. = FALSE)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

## packages ---------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(scales)
})

## labels, helpers (same definitions as 10b_results_assets_landcover.R) ---------
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

## scaling record + staleness check ---------------------------------------------
## The scaled panel was built from the canonical panel as it stood then. If the
## canonical pipeline has been rerun since, the high-hazard run is stale.
scaling   <- jsonlite::read_json(ALT_SCALING, simplifyVector = TRUE)
md5_match <- if (file.exists(BASE_PANEL))
  unname(tools::md5sum(BASE_PANEL)) == scaling$source_panel_md5 else NA
if (isFALSE(md5_match))
  warning("The canonical data_panel.rds has changed since the high-hazard panel was ",
          "built. Rebuild it (PPR_HAZARD=", scaling$target_annual_pct,
          " Rscript 14_build_hazard_panel.R) and rerun its models before reporting.",
          call. = FALSE)

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

tbl_haz_headline <- cmp %>%
  filter(scenario %in% CLIM) %>%
  transmute(scenario,
            rg_pct_base, rg_pct_alt, rg_pct_diff = rg_pct_alt - rg_pct_base,
            rm_pct_base, rm_pct_alt, rm_pct_diff = rm_pct_alt - rm_pct_base,
            ranking_holds = order_ok_base & order_ok_alt)
write_csv(tbl_haz_headline, file.path(OUT_DIR, "haz_headline.csv"))

## acquisition value (value added over the do-nothing baseline) by policy, both runs
va_of <- function(f) read_csv(f, show_col_types = FALSE) %>%
  mutate(scenario = as.character(scenario)) %>%
  select(scenario, model, value_added) %>%
  pivot_wider(names_from = model, values_from = value_added)
ratio <- function(a, b) ifelse(abs(b) < 1, NA_real_, a / b)   # high hazard / canonical
tbl_haz_value <- inner_join(va_of(BASE_RES), va_of(ALT_RES), by = "scenario",
                            suffix = c("_base", "_alt")) %>%
  mutate(scenario = factor(scenario, levels = SCEN_LEVELS)) %>%
  arrange(scenario) %>%
  transmute(scenario,
            greedy_base,  greedy_alt,  greedy_ratio  = ratio(greedy_alt,  greedy_base),
            myopic_base,  myopic_alt,  myopic_ratio  = ratio(myopic_alt,  myopic_base),
            rolling_base, rolling_alt, rolling_ratio = ratio(rolling_alt, rolling_base))
write_csv(tbl_haz_value, file.path(OUT_DIR, "haz_value_added.csv"))

# ═══════════════════════════════════════════════════════════════════════════════
# SCHEDULE OVERLAP  —  same EAUs bought under both hazard levels?
# Greedy ignores the hazard, so its overlap must be 1 (a built-in check).
# ═══════════════════════════════════════════════════════════════════════════════
drew_overlap <- FALSE
if (file.exists(BASE_SCHED) && file.exists(ALT_SCHED)) {
  acq_of <- function(f) read_csv(f, show_col_types = FALSE) %>%
    mutate(scenario = as.character(scenario)) %>%
    filter(!is.na(acquired_period)) %>%
    distinct(scenario, model, eau_id)
  ab <- acq_of(BASE_SCHED) %>% mutate(base = TRUE)
  aa <- acq_of(ALT_SCHED)  %>% mutate(alt  = TRUE)
  tbl_haz_overlap <- full_join(ab, aa, by = c("scenario", "model", "eau_id")) %>%
    mutate(base = coalesce(base, FALSE), alt = coalesce(alt, FALSE)) %>%
    group_by(scenario, model) %>%
    summarise(n_base = sum(base), n_alt = sum(alt), n_both = sum(base & alt),
              jaccard = n_both / sum(base | alt), .groups = "drop") %>%
    mutate(scenario = factor(scenario, levels = SCEN_LEVELS),
           model    = factor(model, levels = c("greedy", "myopic", "rolling"))) %>%
    arrange(scenario, model)
  write_csv(tbl_haz_overlap, file.path(OUT_DIR, "haz_overlap.csv"))
  drew_overlap <- TRUE
  g_off <- tbl_haz_overlap %>% filter(model == "greedy", jaccard < 1)
  if (nrow(g_off))
    warning("Greedy schedules differ between the runs (", nrow(g_off), " scenario(s)). ",
            "Greedy does not use the hazard, so the runs differ in something else.",
            call. = FALSE)
} else {
  message("Overlap table skipped: schedule missing ('", BASE_SCHED, "' or '",
          ALT_SCHED, "'). Rerun 08 for that run.")
}

# ═══════════════════════════════════════════════════════════════════════════════
# CONVERSION RISK  —  as Fig. 3 of the main results (RCP 4.5/8.5, 2100 excluded)
# ═══════════════════════════════════════════════════════════════════════════════
risk_of <- function(pnl) {
  pnl %>% distinct(eau_id, year, rcp, value = trans_prob) %>%
    filter(year < 2100, rcp %in% c("45", "85")) %>%
    mutate(rcp_lab = factor(recode(rcp, `45` = "RCP 4.5", `85` = "RCP 8.5"),
                            levels = c("RCP 4.5", "RCP 8.5")),
           pct = value * 100)
}

drew_risk <- FALSE
if (file.exists(BASE_PANEL) && file.exists(ALT_PANEL)) {
  rk <- bind_rows(risk_of(readRDS(BASE_PANEL)) %>% mutate(run = "base"),
                  risk_of(readRDS(ALT_PANEL))  %>% mutate(run = "alt")) %>%
    mutate(run = factor(run, levels = c("base", "alt")))

  ## figure: boxes dodged by run; diamond = mean, dash = maximum
  dg <- position_dodge(width = .75)
  fig_haz_risk_box <- ggplot(rk, aes(factor(year), pct, fill = run, colour = run)) +
    geom_boxplot(aes(group = interaction(year, run)), outlier.shape = NA,
                 position = dg, width = .62, alpha = .35, linewidth = .3) +
    stat_summary(aes(group = run), fun = max, geom = "point", shape = 95, size = 6,
                 position = dg, show.legend = FALSE) +
    stat_summary(aes(group = run), fun = mean, geom = "point", shape = 23, size = 1.6,
                 fill = "white", position = dg, show.legend = FALSE) +
    facet_wrap(~ rcp_lab, ncol = 2) +
    scale_fill_manual(values = HAZ_RUN_COL, labels = HAZ_SHORT) +
    scale_colour_manual(values = HAZ_RUN_COL, labels = HAZ_SHORT) +
    scale_y_continuous("Per-period conversion risk",
                       labels = label_number(suffix = "%", accuracy = .1),
                       expand = expansion(mult = c(.02, .06))) +
    scale_x_discrete("year") +
    labs(caption = "Box = median and interquartile range; diamond = mean; dash = maximum.") +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = .5, size = rel(.8)),
          plot.caption = element_text(colour = "grey40", size = rel(.72), hjust = 0))
  save_fig(fig_haz_risk_box, "fig_haz_risk_box", 11, 5.0)
  drew_risk <- TRUE
} else {
  message("Risk figure skipped: inputs panel missing ('", BASE_PANEL,
          "' or '", ALT_PANEL, "').")
}

# ═══════════════════════════════════════════════════════════════════════════════
# METADATA  —  11_results.qmd includes the hazard section only if this exists
# ═══════════════════════════════════════════════════════════════════════════════
jsonlite::write_json(
  list(generated_date    = format(Sys.Date(), "%Y-%m-%d"),
       base_run          = RUN_ID, base_short = HAZ_SHORT[["base"]],
       alt_run           = COMPARE_RUN, alt_short = HAZ_SHORT[["alt"]],
       target_annual_pct = scaling$target_annual_pct,
       target_decadal    = scaling$target_decadal,
       k                 = scaling$k,
       k_unweighted      = scaling$k_unweighted,
       canonical_mean_2020_weighted = scaling$canonical_mean_2020_weighted,
       scaled_mean_2020_weighted    = scaling$scaled_mean_2020_weighted,
       source_panel_md5_match       = md5_match,
       thresh_pp         = HAZ_THRESH_PP,
       overlap_drawn     = drew_overlap,
       risk_drawn        = drew_risk),
  file.path(OUT_DIR, "haz_metadata.json"), auto_unbox = TRUE, pretty = TRUE,
  na = "null", digits = NA)

message("\nDone. Comparison of '", RUN_ID, "' (canonical) with '", COMPARE_RUN,
        "' written to '", OUT_DIR, "/': 2 tables",
        if (drew_overlap) " + overlap table" else "",
        if (drew_risk) " + risk figure" else "", " + haz_metadata.json.")
