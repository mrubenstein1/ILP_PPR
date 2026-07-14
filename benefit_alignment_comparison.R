# ══════════════════════════════════════════════════════════════════════════════
# benefit_alignment_comparison.R  —  Layer-1 (pre-solver) benefit comparison
# ══════════════════════════════════════════════════════════════════════════════
#
# Compares the benefit trajectory under the two temporal alignments produced by
# 04__benefit_data.R, BEFORE any model is run. It touches only the benefit columns
# (scaled_abundance / abs_abundance); no Gurobi, no 05/06, no survival or cost.
#
# Prerequisite: run 04 under BOTH alignments so the tagged benefit panels exist:
#     set ALIGNMENT <- "current"  in 04, run it  -> data_panel_benefit_current.rds
#     set ALIGNMENT <- "midpoint" in 04, run it  -> data_panel_benefit_midpoint.rds
#
# Produces three views (the three things model behaviour actually depends on):
#   (1) LEVEL   — distribution of benefit values per year x scenario (mean/median/
#                 spread). Moves J_baseline and value_added scale, but not
#                 necessarily the model's decisions.
#   (2) SHAPE   — the benefit TRAJECTORY over 2020-2100, landscape-total and
#                 per-WMD. This is what myopic freezes and mis-anticipates, so it
#                 is the strongest pre-solver predictor of whether the foresight
#                 gap will move. Front-loading (midpoint pulls anchors earlier and
#                 flattens the post-2085 tail) shows up here.
#   (3) TIMING  — the DISCOUNTED benefit stream Sum_t delta^t * b (no survival, no
#                 protection: a pure benefit-side proxy). Earlier-arriving benefit
#                 carries more discounted weight, so this is the pre-solver proxy
#                 for "does the alignment change what foresight is worth".
#
# Outputs: CSV tables + (if ggplot2 is available) figures, under OUT_DIR.
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble); library(stringr)
})

# ── Config ────────────────────────────────────────────────────────────────────
IN_DIR  <- "input_data"
OUT_DIR <- "output_figs/benefit_alignment"
DELTA   <- 0.95                     # per-decade discount (match 07's DELTA)
DRAW_PLOTS      <- TRUE             # set FALSE for tables only
PERWMD_SCENARIO <- "rcp85_dry"      # which scenario the per-WMD trajectory panel uses

# Scenario coordinates (rcp, gcm), matching 07's SCENARIOS. Each scenario's
# trajectory is the shared 2020 baseline row (rcp = "baseline") stitched onto the
# scenario's own 2030-2100 rows — exactly how build_scenario_matrices assembles it.
SCEN <- list(
  rcp45_wet  = c(rcp = "45",         gcm = "wet"),
  rcp45_dry  = c(rcp = "45",         gcm = "dry"),
  rcp85_wet  = c(rcp = "85",         gcm = "wet"),
  rcp85_dry  = c(rcp = "85",         gcm = "dry"),
  stationary = c(rcp = "stationary", gcm = "stationary")
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ── 0. Load the two tagged benefit panels ─────────────────────────────────────
load_alignment <- function(tag) {
  f <- file.path(IN_DIR, sprintf("data_panel_benefit_%s.rds", tag))
  if (!file.exists(f))
    stop("Missing '", f, "'. Run 04__benefit_data.R with ALIGNMENT <- \"", tag,
         "\" first (both 'current' and 'midpoint' are required).")
  readRDS(f) %>%
    select(eau_id, wmd_id, year, rcp, gcm, prop_suitable,
           abs_abundance, scaled_abundance) %>%
    mutate(alignment = tag)
}
panel <- bind_rows(load_alignment("current"), load_alignment("midpoint"))
panel$alignment <- factor(panel$alignment, levels = c("current", "midpoint"))
cat(sprintf("Loaded benefit panels: %d rows total (%d EAUs, %d WMDs)\n",
            nrow(panel), n_distinct(panel$eau_id), n_distinct(panel$wmd_id)))


# ── Helper: assemble one scenario's rows (2020 baseline + scenario future) ─────
assemble_scenario <- function(dp, sc_name) {
  rcp_s <- SCEN[[sc_name]][["rcp"]]; gcm_s <- SCEN[[sc_name]][["gcm"]]
  base <- dp %>% filter(year == 2020, rcp == "baseline")
  fut  <- dp %>% filter(year  > 2020, rcp == rcp_s, gcm == gcm_s)
  bind_rows(base, fut) %>% mutate(scenario = sc_name)
}
# Long panel of every scenario x alignment, with the baseline stitched in.
scen_panel <- bind_rows(lapply(names(SCEN), function(sn)
  bind_rows(lapply(c("current", "midpoint"), function(al)
    assemble_scenario(panel %>% filter(alignment == al), sn) %>%
      mutate(alignment = al)))))
scen_panel$alignment <- factor(scen_panel$alignment, levels = c("current", "midpoint"))
scen_panel$scenario  <- factor(scen_panel$scenario,  levels = names(SCEN))


# ══ (1) LEVEL — distribution of benefit values per year x scenario ════════════
# EAU-level scaled_abundance distribution (the coefficient the model uses).
level_stats <- scen_panel %>%
  group_by(alignment, scenario, year) %>%
  summarise(
    n_eau      = dplyr::n(),
    mean_b     = mean(scaled_abundance),
    median_b   = median(scaled_abundance),
    sd_b       = sd(scaled_abundance),
    p10_b      = quantile(scaled_abundance, 0.10),
    p90_b      = quantile(scaled_abundance, 0.90),
    total_b    = sum(scaled_abundance),      # landscape total abundance
    .groups = "drop"
  )
write_csv(level_stats, file.path(OUT_DIR, "level_stats_by_year_scenario.csv"))

# Compact side-by-side of the landscape total (the headline level number).
level_wide <- level_stats %>%
  select(alignment, scenario, year, total_b) %>%
  pivot_wider(names_from = alignment, values_from = total_b) %>%
  mutate(diff_mid_minus_cur = midpoint - current,
         pct_change = 100 * (midpoint - current) / current)
write_csv(level_wide, file.path(OUT_DIR, "landscape_total_by_year.csv"))
cat("\n== Landscape-total abundance: current vs midpoint (by scenario x year) ==\n")
print(as.data.frame(level_wide %>% mutate(across(where(is.numeric), ~round(.x, 2)))),
      row.names = FALSE)


# ══ (2) SHAPE — trajectories ══════════════════════════════════════════════════
# Landscape-total trajectory per scenario x alignment (already in level_stats$total_b).
landscape_traj <- level_stats %>% select(alignment, scenario, year, total_b)

# Per-WMD trajectory (sum of scaled_abundance within WMD = WMD-level abundance).
perwmd_traj <- scen_panel %>%
  group_by(alignment, scenario, wmd_id, year) %>%
  summarise(wmd_total = sum(scaled_abundance), .groups = "drop")
write_csv(perwmd_traj, file.path(OUT_DIR, "perwmd_trajectory.csv"))

# A simple "shape divergence" metric per scenario: how much the normalised
# trajectory (each alignment scaled to its own 2020 value) differs. Isolates SHAPE
# change from LEVEL change — if this is ~0 the alignment only rescaled benefit.
shape_div <- landscape_traj %>%
  group_by(alignment, scenario) %>%
  arrange(year) %>%
  mutate(norm_b = total_b / total_b[year == 2020]) %>%
  ungroup() %>%
  select(alignment, scenario, year, norm_b) %>%
  pivot_wider(names_from = alignment, values_from = norm_b) %>%
  mutate(shape_gap = midpoint - current)
shape_summary <- shape_div %>%
  group_by(scenario) %>%
  summarise(max_abs_shape_gap = max(abs(shape_gap)),
            mean_abs_shape_gap = mean(abs(shape_gap)), .groups = "drop")
write_csv(shape_div, file.path(OUT_DIR, "shape_divergence_normalised.csv"))
cat("\n== Shape divergence (normalised-to-2020 trajectory, midpoint - current) ==\n")
print(as.data.frame(shape_summary %>% mutate(across(where(is.numeric), ~round(.x, 4)))),
      row.names = FALSE)


# ══ (3) TIMING — discounted benefit stream ════════════════════════════════════
# Pure benefit-side proxy: Sum_t delta^t * (landscape-total b). NO survival, NO
# protection weighting (those enter at the solver layer). t = (year - 2020)/10.
disc <- landscape_traj %>%
  mutate(t = (year - 2020) / 10,
         discounted_b = total_b * DELTA^t)

# Cumulative discounted benefit over time (for the timing plot).
disc_cum <- disc %>%
  group_by(alignment, scenario) %>%
  arrange(year) %>%
  mutate(cum_discounted_b = cumsum(discounted_b)) %>%
  ungroup()
write_csv(disc_cum %>% select(alignment, scenario, year, discounted_b, cum_discounted_b),
          file.path(OUT_DIR, "discounted_benefit_stream.csv"))

# Total discounted benefit per scenario x alignment (the headline timing number).
disc_total <- disc %>%
  group_by(alignment, scenario) %>%
  summarise(discounted_total = sum(discounted_b), .groups = "drop") %>%
  pivot_wider(names_from = alignment, values_from = discounted_total) %>%
  mutate(diff_mid_minus_cur = midpoint - current,
         pct_change = 100 * (midpoint - current) / current)
write_csv(disc_total, file.path(OUT_DIR, "discounted_total_by_scenario.csv"))
cat("\n== Discounted benefit stream total (delta =", DELTA, "): current vs midpoint ==\n")
print(as.data.frame(disc_total %>% mutate(across(where(is.numeric), ~round(.x, 2)))),
      row.names = FALSE)


# ══ Figures (skipped gracefully if ggplot2 is absent) ═════════════════════════
if (DRAW_PLOTS && requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  al_cols <- c(current = "#7570b3", midpoint = "#d95f02")
  base_thm <- theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey92", color = NA),
          panel.grid.minor = element_blank(), legend.position = "bottom")

  # (2a) Landscape-total trajectory, faceted by scenario.
  p1 <- ggplot(landscape_traj, aes(year, total_b, color = alignment)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.2) +
    facet_wrap(~ scenario, scales = "free_y") +
    scale_color_manual(values = al_cols) +
    labs(title = "Landscape-total benefit trajectory: current vs midpoint alignment",
         subtitle = "Sum of scaled_abundance across EAUs (2020 baseline stitched onto each scenario)",
         x = NULL, y = "Landscape-total abundance", color = "Alignment") +
    base_thm
  ggsave(file.path(OUT_DIR, "fig1_landscape_trajectory.png"), p1,
         width = 9, height = 6, dpi = 150)

  # (2b) Per-WMD trajectory for one scenario.
  p2 <- perwmd_traj %>% filter(scenario == PERWMD_SCENARIO) %>%
    ggplot(aes(year, wmd_total, color = alignment)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ wmd_id, scales = "free_y") +
    scale_color_manual(values = al_cols) +
    labs(title = paste0("Per-WMD benefit trajectory (", PERWMD_SCENARIO, "): current vs midpoint"),
         x = NULL, y = "WMD-total abundance", color = "Alignment") +
    base_thm
  ggsave(file.path(OUT_DIR, "fig2_perwmd_trajectory.png"), p2,
         width = 10, height = 7, dpi = 150)

  # (3) Cumulative discounted benefit over time, faceted by scenario.
  p3 <- ggplot(disc_cum, aes(year, cum_discounted_b, color = alignment)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ scenario, scales = "free_y") +
    scale_color_manual(values = al_cols) +
    labs(title = "Cumulative discounted benefit stream: current vs midpoint",
         subtitle = paste0("Sum_t delta^t * landscape-total b (delta = ", DELTA,
                           "); benefit-side only, no survival/protection"),
         x = NULL, y = "Cumulative discounted abundance", color = "Alignment") +
    base_thm
  ggsave(file.path(OUT_DIR, "fig3_cumulative_discounted.png"), p3,
         width = 9, height = 6, dpi = 150)

  cat("\nFigures written to", OUT_DIR, "\n")
} else if (DRAW_PLOTS) {
  cat("\n[note] ggplot2 not available — tables written, figures skipped.\n")
}

cat("\nDone. CSV tables + figures in:", OUT_DIR, "\n")
