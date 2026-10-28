# ═══════════════════════════════════════════════════════════════════════════════
# 10e_results_assets_hazard_sweep.R — OPTIONAL: assets across the hazard sweep
# ═══════════════════════════════════════════════════════════════════════════════
# Collects the canonical run (k = 1) and every complete high-hazard run,
# runs/<LC_SCENARIO>_haz<tag>/ (made by 15_run_hazard_sweep.R or by hand), and
# traces how policy performance changes with the level of conversion pressure.
# Each run's rate, k and source-panel md5 are read from its hazard_scaling.json.
#
# The gap between rolling and greedy is split into its two parts:
#   myopic vs greedy  = value of using conversion risk at all (risk frozen at today's level)
#   rolling vs myopic = value of foresight (knowing how risk and benefit change)
#   rolling vs greedy = the two together
#
# USAGE (PPR_HAZARD and PPR_BUDGET unset):
#   Rscript 10e_results_assets_hazard_sweep.R
#
# Writes to the canonical run's output_figs (raw numbers):
#   haz_sweep.csv            one row per run x scenario: rate, k, value by policy,
#                            % improvement and Δpairs for the three comparisons
#   fig_haz_sweep_pct        % improvement vs annual rate, three panels (PNG + PDF)
#   fig_haz_sweep_pairs      Δpairs vs annual rate, three panels (PNG + PDF)
#   fig_haz_sweep_risk_box   per-period risk by decade, one box per run, by RCP (PNG + PDF)
#   haz_sweep_metadata.json  runs included, canonical rate, staleness check
# ═══════════════════════════════════════════════════════════════════════════════

if (!exists(".SETUP_DONE")) source("00_setup.R")

# ── CONFIG ─────────────────────────────────────────────────────────────────────
SWEEP_PATTERN <- paste0("^", LC_SCENARIO, "_haz[0-9]+(p[0-9]+)?$")   # hazard run folders
# Short labels for the top axis, keyed by the rate as typed in PPR_HAZARD.
# Rates not listed here are labelled by their value.
BENCH_LAB <- c(`0.09` = "Dahl 2014",
               `0.28` = "Johnston 2013",
               `0.4`  = "Stephens 2008",
               `0.57` = "Kemink 2023",
               `0.95` = "Rashford 2011 (2006\u201311)",
               `1.33` = "Rashford 2011")
CANON_LAB <- "Canonical"

OUT_DIR   <- DIR_FIGS                        # the canonical run's output_figs
CANON_RES <- file.path(DIR_OUT, "model_results.csv")
CANON_PNL <- file.path(DIR_DERIVED, "data_panel.rds")

DPI         <- 300
BASE_SIZE   <- 12
BASE_FAMILY <- "serif"                       # same as 10_results_assets.R
SAVE_PDF    <- TRUE
# ────────────────────────────────────────────────────────────────────────────────

if (HAZARD_RUN || BUDGET_SWEEP)
  stop("Run 10e from the canonical configuration: unset PPR_HAZARD and PPR_BUDGET.",
       call. = FALSE)
if (!file.exists(CANON_RES))
  stop("Canonical results not found at ", CANON_RES, ". Run the canonical pipeline first.",
       call. = FALSE)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(scales)
})

SCEN_LEVELS <- c("rcp45_wet", "rcp45_dry", "rcp85_wet", "rcp85_dry", "stationary")
CLIM        <- SCEN_LEVELS[1:4]
SCEN_LAB1   <- c(rcp45_wet = "RCP 4.5 wet", rcp45_dry = "RCP 4.5 dry",
                 rcp85_wet = "RCP 8.5 wet", rcp85_dry = "RCP 8.5 dry")
SCEN_COL    <- c(rcp45_wet = "#3B7EA1", rcp45_dry = "#9DBBCB",
                 rcp85_wet = "#9E2B3A", rcp85_dry = "#E0A090")
CMP_LAB     <- c(mg = "Myopic vs. greedy\n(value of using risk)",
                 rm = "Rolling vs. myopic\n(value of foresight)",
                 rg = "Rolling vs. greedy\n(total)")

theme_ppr <- function(base_size = BASE_SIZE, base_family = BASE_FAMILY) {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.title    = element_text(colour = "grey25"),
      axis.title.y  = element_text(angle = 90, margin = margin(r = 6)),
      panel.grid.minor = element_blank(),
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
# RUNS  —  canonical + every complete hazard run
# ═══════════════════════════════════════════════════════════════════════════════
cand <- grep(SWEEP_PATTERN, list.dirs("runs", recursive = FALSE, full.names = FALSE),
             value = TRUE)
runs <- bind_rows(lapply(cand, function(id) {
  js  <- file.path("runs", id, "derived_data", "hazard_scaling.json")
  res <- file.path("runs", id, "output_data", "model_results.csv")
  if (!file.exists(js) || !file.exists(res)) return(NULL)
  s <- jsonlite::read_json(js, simplifyVector = TRUE)
  tibble(run_id = id, res = res, rate_typed = s$target_annual_pct, k = s$k,
         decadal = s$target_decadal, canon_decadal = s$canonical_mean_2020_weighted,
         source_md5 = s$source_panel_md5)
}))
if (!nrow(runs))
  stop("No complete hazard runs matching runs/", LC_SCENARIO, "_haz*. ",
       "Run 15_run_hazard_sweep.R first.", call. = FALSE)

## staleness: every hazard run must come from the canonical panel as it is now
canon_md5 <- if (file.exists(CANON_PNL)) unname(tools::md5sum(CANON_PNL)) else NA_character_
runs <- runs %>% mutate(stale = !is.na(canon_md5) & source_md5 != canon_md5)
if (any(runs$stale))
  warning("Built from an older canonical panel (rerun them with 15_run_hazard_sweep.R): ",
          paste(runs$run_id[runs$stale], collapse = ", "), call. = FALSE)

## canonical point: its weighted 2020 rate, as recorded by 14 (identical across runs
## built from the same canonical panel)
canon_decadal <- runs$canon_decadal[!runs$stale][1]
if (is.na(canon_decadal)) canon_decadal <- runs$canon_decadal[1]
to_annual <- function(d) 100 * (1 - (1 - d)^(1 / 10))     # decadal prob -> annual %

pts <- bind_rows(
  tibble(run_id = RUN_ID, res = CANON_RES, rate = to_annual(canon_decadal), k = 1,
         decadal = canon_decadal, label = CANON_LAB, stale = FALSE),
  runs %>% transmute(run_id, res, rate = rate_typed, k, decadal,
                     # format() each rate on its own: on a vector it pads to a common width
                     label = coalesce(unname(BENCH_LAB[vapply(rate_typed, format, "")]),
                                      paste0(vapply(rate_typed, format, ""), "%")),
                     stale)
) %>% arrange(rate)
message("Runs: ", paste0(pts$run_id, " (", signif(pts$rate, 3), "%/yr, k = ",
                         signif(pts$k, 3), ")", collapse = ", "))

# ═══════════════════════════════════════════════════════════════════════════════
# TABLE  —  value by policy and the three comparisons, per run x scenario
# ═══════════════════════════════════════════════════════════════════════════════
pct_improve <- function(a, b, d) ifelse(abs(d) < 1 | abs(b) < 1e-9, NA_real_, 100 * (a - b) / b)

sweep <- bind_rows(lapply(seq_len(nrow(pts)), function(i) {
  read_csv(pts$res[i], show_col_types = FALSE) %>%
    mutate(scenario = as.character(scenario)) %>%
    select(scenario, model, value_added) %>%
    pivot_wider(names_from = model, values_from = value_added) %>%
    mutate(run_id = pts$run_id[i], label = pts$label[i], annual_rate_pct = pts$rate[i],
           decadal_rate = pts$decadal[i], k = pts$k[i], stale = pts$stale[i])
})) %>%
  transmute(run_id, label, annual_rate_pct, decadal_rate, k, stale,
            scenario = factor(scenario, levels = SCEN_LEVELS),
            greedy, myopic, rolling,
            mg_pairs = myopic - greedy, rm_pairs = rolling - myopic, rg_pairs = rolling - greedy,
            mg_pct = pct_improve(myopic,  greedy, mg_pairs),
            rm_pct = pct_improve(rolling, myopic, rm_pairs),
            rg_pct = pct_improve(rolling, greedy, rg_pairs),
            ranking_holds = rolling >= myopic - 1e-6 & myopic >= greedy - 1e-6) %>%
  arrange(annual_rate_pct, scenario)
write_csv(sweep, file.path(OUT_DIR, "haz_sweep.csv"))

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURES  —  the three comparisons against the annual rate (climate scenarios)
# x-axis labels give the annual rate and the benchmark each run stands for.
# ═══════════════════════════════════════════════════════════════════════════════
long_of <- function(kind) {
  sweep %>% filter(scenario %in% CLIM) %>%
    select(annual_rate_pct, scenario, ends_with(paste0("_", kind))) %>%
    pivot_longer(-c(annual_rate_pct, scenario), names_to = "cmp", values_to = "value") %>%
    mutate(cmp = factor(sub(paste0("_", kind), "", cmp), levels = names(CMP_LAB)))
}
x_scale <- scale_x_continuous("Annual conversion rate (%)", breaks = pts$rate,
                              labels = paste0(sprintf("%.2f", pts$rate), "  ", pts$label))
sweep_plot <- function(d, y_lab, y_labels) {
  ggplot(d, aes(annual_rate_pct, value, colour = scenario)) +
    geom_vline(xintercept = pts$rate, colour = "grey88", linewidth = .3) +
    geom_line(linewidth = .7) + geom_point(size = 1.6) +
    facet_wrap(~ cmp, ncol = 3, scales = "free_y", labeller = as_labeller(CMP_LAB)) +
    scale_colour_manual(values = SCEN_COL, labels = SCEN_LAB1) +
    x_scale +
    scale_y_continuous(y_lab, labels = y_labels, limits = c(0, NA),
                       expand = expansion(mult = c(0, .06))) +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = rel(.75)))
}
save_fig(sweep_plot(long_of("pct"), "% improvement", label_number(suffix = "%")),
         "fig_haz_sweep_pct", 13, 5.2)
save_fig(sweep_plot(long_of("pairs"), "\u0394pairs (discounted)", label_comma()),
         "fig_haz_sweep_pairs", 13, 5.2)

# ═══════════════════════════════════════════════════════════════════════════════
# CONVERSION RISK  —  per-period risk by decade, one box per run, panels by RCP.
# trans_prob varies by RCP only (not wet/dry); 2100 has no next period, so excluded.
# ═══════════════════════════════════════════════════════════════════════════════
pnl_of <- function(id) if (id == RUN_ID) CANON_PNL else
  file.path("runs", id, "derived_data", "data_panel.rds")
risk_pnl <- vapply(pts$run_id, pnl_of, "")
if (!all(file.exists(risk_pnl))) {
  message("Risk figure skipped: panel missing for ",
          paste(pts$run_id[!file.exists(risk_pnl)], collapse = ", "), ".")
} else {
  rate_lab <- paste0(sprintf("%.2f", pts$rate), "%/yr")
  rk <- bind_rows(lapply(seq_len(nrow(pts)), function(i) {
    readRDS(risk_pnl[i]) %>%
      distinct(eau_id, year, rcp, value = trans_prob) %>%
      filter(year < 2100, rcp %in% c("45", "85")) %>%
      mutate(run = rate_lab[i])
  })) %>%
    mutate(rcp_lab = factor(recode(rcp, `45` = "RCP 4.5", `85` = "RCP 8.5"),
                            levels = c("RCP 4.5", "RCP 8.5")),
           run = factor(run, levels = rate_lab),          # pts is sorted by rate
           pct = value * 100)

  dg <- position_dodge(width = .85)
  fig_haz_sweep_risk_box <- ggplot(rk, aes(factor(year), pct, fill = run, colour = run)) +
    geom_boxplot(aes(group = interaction(year, run)), outlier.shape = NA,
                 position = dg, width = .78, alpha = .35, linewidth = .3) +
    stat_summary(aes(group = run), fun = max, geom = "point", shape = 95, size = 4,
                 position = dg, show.legend = FALSE) +
    stat_summary(aes(group = run), fun = mean, geom = "point", shape = 23, size = 1.2,
                 fill = "white", position = dg, show.legend = FALSE) +
    facet_wrap(~ rcp_lab, ncol = 1) +
    scale_fill_viridis_d(option = "rocket", begin = .85, end = .15) +
    scale_colour_viridis_d(option = "rocket", begin = .85, end = .15) +
    scale_y_continuous("Per-period conversion risk",
                       labels = label_number(suffix = "%", accuracy = 1),
                       expand = expansion(mult = c(.02, .06))) +
    scale_x_discrete("year") +
    labs(caption = "Box = median and interquartile range; diamond = mean; dash = maximum.") +
    theme(panel.grid.major.x = element_blank(),
          plot.caption = element_text(colour = "grey40", size = rel(.72), hjust = 0))
  save_fig(fig_haz_sweep_risk_box, "fig_haz_sweep_risk_box", 13, 8)
}

# ═══════════════════════════════════════════════════════════════════════════════
# METADATA
# ═══════════════════════════════════════════════════════════════════════════════
jsonlite::write_json(
  list(generated_date       = format(Sys.Date(), "%Y-%m-%d"),
       canonical_run        = RUN_ID,
       canonical_annual_pct = to_annual(canon_decadal),
       runs                 = pts %>% select(run_id, annual_rate_pct = rate, k, label, stale),
       any_stale            = any(pts$stale)),
  file.path(OUT_DIR, "haz_sweep_metadata.json"), auto_unbox = TRUE, pretty = TRUE,
  digits = NA)

message("\nDone. ", nrow(pts), " runs (canonical + ", nrow(runs), " hazard) written to '",
        OUT_DIR, "/': haz_sweep.csv + figures + haz_sweep_metadata.json.")
