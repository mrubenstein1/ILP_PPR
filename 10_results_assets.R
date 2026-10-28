# ═══════════════════════════════════════════════════════════════════════════════
# 10_results_assets.R  —  Asset generator (figs, tables, maps) for the PPR acquisition results document
# ═══════════════════════════════════════════════════════════════════════════════
# Writes to OUT_DIR:
#
# Tables are written as tidy CSVs of RAW NUMBERS; all formatting and the grouped-header
# styling happen in results.qmd. Same approach with Figs/Maps
#
# After this run's own assets, three COMPARISON SECTIONS (formerly the separate
# scripts 10b, 10c and 10d) run at the end of this file:
#   land-cover comparison   this run vs the BAU land-cover run      (cmp_* files)
#   budget sensitivity      the budget runs from 13_run_budget_sweep.R  (bud_* files)
#   hazard sweep            the hazard runs from 15_run_hazard_sweep.R  (haz_sweep* files)
# Each runs only in the canonical configuration and only when its input runs exist;
# otherwise it is skipped with a message. See each section's own header.

# ═══════════════════════════════════════════════════════════════════════════════

if (!exists(".SETUP_DONE")) source("00_setup.R")

# ── CONFIG ─────────────────────────────────────────────────────────────────────
IN_DIR   <- DIR_OUT      # runs/<RUN_ID>/output_data
OUT_DIR  <- DIR_FIGS     # runs/<RUN_ID>/output_figs
RESULTS  <- file.path(IN_DIR, "model_results.csv")
TRAJ     <- file.path(IN_DIR, "model_trajectories.csv")
SCHEDULE <- file.path(IN_DIR, "acquisition_schedule_spatial.csv")

DPI         <- 300
BASE_SIZE   <- 12
# Figure/map text (axes, legends, strip titles, colour-bar, on-plot labels) is drawn
# in this family. "serif" is a portable serif (Times-like) that matches the PDF's
# roman body font across every ggsave backend (grDevices, cairo, ragg) with no extra
# packages. For a *pixel-exact* match to the PDF body (Latin Modern Roman), see the
# note in results.qmd and swap in a registered "Latin Modern Roman" family instead.
BASE_FAMILY <- "serif"
SAVE_PDF    <- TRUE           # also write vector PDF alongside each PNG

# palette (models) -------------------------------------------------------------
COL_MODEL <- c(greedy = "#B0413E", myopic = "#E2A23B", rolling = "#2E7D5B")

# maps -------------------------------------------------------------------------
DRAW_MAPS    <- TRUE                       # set FALSE to skip maps (tables + figs only)
# EAU tile size in metres = res() of derived_data/wmd_raster_equal_area.tif
# (18,513.54 x 18,485.21 m, ~342.2 km^2; cells are very slightly non-square)
CELL_W       <- 18513.54                  # EAU width  (x resolution), metres
CELL_H       <- 18485.21                  # EAU height (y resolution), metres
MAP_BG       <- "grey91"                  # unacquired-EAU landscape colour
WMD_OVERLAY  <- TRUE                       # draw WMD outlines if raster + terra/sf available
WMD_RASTER   <- file.path(DIR_COMMON, "wmd_raster_equal_area.tif")
WMD_LINE_COL <- "grey25"
WMD_LINE_W   <- 0.25
PAIR_COL <- c("Acquired by both"               = "#5E4FA2",
              "Rolling only (foresight adds)"  = "#2C7FB8",
              "Myopic only (foresight avoids)" = "#E6892B",
              "Neither"                        = MAP_BG)
# per-parcel prevented-loss difference (ΔV = rolling − myopic), diverging
DIFF_POS   <- "#2C7FB8"   # blue: rolling banks MORE prevented loss here
DIFF_NEG   <- "#B0413E"   # red:  myopic banks more (budget landed here, not on a blue cell)
DIFF_CAP_Q <- 0.98        # symmetric colour cap at this quantile of |ΔV| (squish beyond)

# descriptive-summary assets (benefit / cost / risk inputs) --------------------
PANEL_RDS     <- file.path(DIR_DERIVED, "data_panel.rds")      # inputs panel from script 06
EAU_LOOKUP    <- file.path(DIR_COMMON, "eau_wmd_lookup.rds")  # eau_id -> x_coord/y_coord
MAP_YEAR      <- 2020        # decision period shown by the cost map
N_COST_BINS   <- 6           # equal-count classes for the cost map
COST_RAMP     <- "rocket"    # viridis option, warm -> cost
BEN_RAMP      <- "mako"      # viridis option, cool -> benefit
BEN_CAP_Q     <- 1           # 1 = no cap (colour strictly proportional to pairs);
                             # set e.g. .98 to squish the top 2% if the ramp reads flat
BEN_MAP_SCEN  <- "rcp85_dry" # scenario shown in the 2100 panel of the benefit map
BEN_MAP_YEARS <- c(2020, 2100)
MAP_NA_COL    <- "grey88"    # EAUs with no estimate
EPS_LAB       <- "background floor"   # on-plot wording for the epsilon floor
# ────────────────────────────────────────────────────────────────────────────────

## packages ---------------------------------------------------------------------
.need <- c("readr", "dplyr", "tidyr", "ggplot2", "scales", "jsonlite")
.miss <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
if (length(.miss))
  stop("Install missing package(s): ", paste(.miss, collapse = ", "),
       "\n  install.packages(c(", paste(sprintf('\"%s\"', .miss), collapse = ", "), "))")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(scales)
})

if (!file.exists(RESULTS) || !file.exists(TRAJ))
  stop("Could not find inputs in IN_DIR = '", IN_DIR, "'.\n  Expected: ",
       RESULTS, " and ", TRAJ, "\n  Edit IN_DIR at the top if they live elsewhere.")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

DELTA <- "\u0394"  # Δ, written portably

## labels / order ---------------------------------------------------------------
SCEN_LEVELS <- c("rcp45_wet", "rcp45_dry", "rcp85_wet", "rcp85_dry", "stationary")
SCEN_LAB1   <- c(rcp45_wet = "RCP 4.5 wet", rcp45_dry = "RCP 4.5 dry",
                 rcp85_wet = "RCP 8.5 wet", rcp85_dry = "RCP 8.5 dry",
                 stationary = "Stationary")
CLIM <- SCEN_LEVELS[1:4]

## formatter used for on-plot labels (table formatting lives in results.qmd) -----
fmt_pairs <- function(x) ifelse(abs(x) < 1, "<1", formatC(round(x), format = "d", big.mark = ","))

## theme ------------------------------------------------------------------------
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

## on-plot text (e.g. the end-of-line value labels) doesn't read the theme's
## base_family, so set it explicitly here to keep every glyph in the same font.
ggplot2::update_geom_defaults("text",  list(family = BASE_FAMILY))
ggplot2::update_geom_defaults("label", list(family = BASE_FAMILY))

# Saves PNG (300 dpi) + optional vector PDF. Prefers cairo_pdf for font embedding and
# falls back to the base 'pdf' device (catching the Windows "no cairo" warning) quietly.
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

## data -------------------------------------------------------------------------
res <- read_csv(RESULTS, show_col_types = FALSE) %>%
  mutate(scenario = factor(scenario, levels = SCEN_LEVELS))
trj <- read_csv(TRAJ, show_col_types = FALSE) %>%
  mutate(scenario = factor(scenario, levels = SCEN_LEVELS))

rollVA <- res %>% filter(model == "rolling") %>% select(scenario, rollVA = value_added)
gaps <- res %>%
  filter(model %in% c("greedy", "myopic")) %>%
  left_join(rollVA, by = "scenario") %>%
  mutate(dpairs = rollVA - value_added,
         model  = factor(model, levels = c("greedy", "myopic")))

# ═══════════════════════════════════════════════════════════════════════════════
# TABLES  —  tidy CSVs of RAW NUMBERS (formatting + grouped-header styling: results.qmd)
# ═══════════════════════════════════════════════════════════════════════════════


# Table 1 — improvement of rolling foresight over greedy and myopic
# Value added (VA) = J − J_baseline: cumulative discounted duck pairs secured relative
# to the do-nothing baseline, 2020–2100. J_baseline is common to all policies within a
# scenario, so rolling − model differences are identical whether taken on VA or on J.
#   *_va           : each policy's value added (pairs)
#   rg_/rm_dpairs  : rolling VA − greedy/myopic VA (pairs)
#   rg_/rm_pct     : % improvement of rolling over greedy/myopic, i.e. dpairs divided
#                    by the COMPARISON policy's VA (not rolling's). NA if that VA is ~0.
pct_improve <- function(roll, base) ifelse(abs(base) < 1e-9, NA_real_, 100 * (roll - base) / base)
pct_short <- function(roll, base) ifelse(abs(roll) < 1e-9, NA_real_, 100 * (roll - base) / roll)

tbl_rolling_improvement <- res %>%
  select(scenario, model, value_added) %>%
  pivot_wider(names_from = model, values_from = value_added) %>%
  arrange(scenario) %>%
  transmute(scenario,
            greedy_va  = greedy, myopic_va = myopic, rolling_va = rolling,
            rg_dpairs  = rolling - greedy, rg_pct = pct_improve(rolling, greedy),
            rm_dpairs  = rolling - myopic, rm_pct = pct_improve(rolling, myopic),
            rg_short_pct = pct_short(rolling, greedy),   # % of ROLLING's VA
            rm_short_pct = pct_short(rolling, myopic))
write_csv(tbl_rolling_improvement, file.path(OUT_DIR, "table_rolling_improvement.csv"))

# Table 2 — value created over the do-nothing baseline (all three models)
tbl_over_baseline <- res %>%
  transmute(scenario, model,
            pct    = value_added / J_baseline * 100,   # % improvement over baseline
            dpairs = value_added) %>%                  # Δpairs over baseline (= value_added)
  pivot_wider(names_from = model, values_from = c(pct, dpairs)) %>%
  arrange(scenario) %>%
  transmute(scenario,
            greedy_pct  = pct_greedy,   greedy_dpairs  = dpairs_greedy,
            myopic_pct  = pct_myopic,   myopic_dpairs  = dpairs_myopic,
            rolling_pct = pct_rolling,  rolling_dpairs = dpairs_rolling)
write_csv(tbl_over_baseline, file.path(OUT_DIR, "table_over_baseline.csv"))

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE (→ Figure 1)  —  cumulative performance relative to rolling (climate scenarios)
# Rolling is the benchmark and sits at zero by construction; greedy and myopic are
# plotted as (model − rolling), so they fall below zero as their shortfall accumulates.
# ═══════════════════════════════════════════════════════════════════════════════
roll_disc <- trj %>% filter(model == "rolling") %>% select(scenario, year, roll = discounted)
cumg <- trj %>%
  filter(model %in% c("greedy", "myopic", "rolling"), scenario %in% CLIM) %>%
  left_join(roll_disc, by = c("scenario", "year")) %>%
  arrange(scenario, model, year) %>%
  group_by(scenario, model) %>%
  mutate(cum_gap = cumsum(discounted - roll)) %>%     # rolling: exactly 0
  ungroup() %>%
  mutate(model = factor(model, levels = c("rolling", "myopic", "greedy")),
         scenario = factor(scenario, levels = CLIM))
ends <- cumg %>% filter(model != "rolling") %>%
  group_by(scenario, model) %>% filter(year == max(year)) %>% ungroup()

fig_overtime <- ggplot(cumg, aes(year, cum_gap, colour = model)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  geom_text(data = ends, aes(label = fmt_pairs(cum_gap)), hjust = -0.12, size = 2.8,
            show.legend = FALSE) +
  facet_wrap(~scenario, ncol = 2, labeller = as_labeller(SCEN_LAB1)) +
  scale_colour_manual(values = COL_MODEL, breaks = c("rolling", "myopic", "greedy"),
                      labels = c("Rolling (benchmark)", "Myopic", "Greedy")) +
  scale_x_continuous("year", expand = expansion(mult = c(0.02, 0.14))) +
  scale_y_continuous("Cumulative difference vs. rolling\n(discounted duck-pairs)",
                     labels = label_comma()) +
  theme(panel.grid.major.x = element_line(colour = "grey92"))
save_fig(fig_overtime, "overtime_gap", 10.5, 7.4)

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE (→ Figure 2)  —  landscape abundance over time (rolling)
# ═══════════════════════════════════════════════════════════════════════════════
SCEN_COL <- c(rcp45_wet = "#3F7FA6", rcp45_dry = "#9BB7C4",
              rcp85_wet = "#9C3848", rcp85_dry = "#D79A86", stationary = "#777777")
dec <- trj %>% filter(model == "rolling") %>%
  mutate(scenario = factor(scenario, levels = SCEN_LEVELS))

fig_decline <- ggplot(dec, aes(year, ducks / 1e6, colour = scenario)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
  scale_colour_manual(values = SCEN_COL, labels = SCEN_LAB1) +
  scale_y_continuous("Landscape abundance\n(million duck-pairs)") +
  scale_x_continuous("year", breaks = seq(2020, 2100, 20)) +
  guides(colour = guide_legend(nrow = 2)) +
  theme(panel.grid.major.x = element_line(colour = "grey92"))
save_fig(fig_decline, "decline", 9, 5.2)

# ═══════════════════════════════════════════════════════════════════════════════
# SHARED MAP MACHINERY  —  used by the descriptive maps below AND by Figures 3–4
# Hoisted to top level (it previously lived inside the Figures 3–4 block) so both
# map families share one WMD overlay, one theme, and one coordinate system. The
# descriptive maps need this machinery but not the acquisition schedule, so the
# "active WMD" filter now reads EAU centroids from the script-01 lookup rather than
# from the schedule. Both cover the full landscape, so the outlines are unchanged.
# ═══════════════════════════════════════════════════════════════════════════════
eau_xy <- if (file.exists(EAU_LOOKUP))
  readRDS(EAU_LOOKUP)[, c("eau_id", "x_coord", "y_coord")] else NULL

wmd_sf <- NULL; map_crs <- NULL
if (WMD_OVERLAY && file.exists(WMD_RASTER) &&
    requireNamespace("terra", quietly = TRUE) && requireNamespace("sf", quietly = TRUE)) {
  ok <- tryCatch({
    r     <- terra::rast(WMD_RASTER)
    polys <- terra::as.polygons(r, dissolve = TRUE)
    if (!is.null(eau_xy)) {          # keep only WMDs that actually contain EAUs
      active <- terra::extract(r, as.matrix(eau_xy[, c("x_coord", "y_coord")]))[, 1]
      active <- sort(unique(active[!is.na(active)]))
      polys  <- polys[polys[[names(polys)[1]]] %in% active, ]
    }
    wmd_sf <<- sf::st_as_sf(polys); map_crs <<- sf::st_crs(wmd_sf); TRUE
  }, error = function(e) { message("WMD overlay skipped: ", conditionMessage(e)); FALSE })
  if (!ok) { wmd_sf <- NULL; map_crs <- NULL }
} else if (WMD_OVERLAY) {
  message("WMD overlay: raster or terra/sf unavailable — drawing maps without outlines.")
}

theme_map <- function() {
  theme_void(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
    theme(strip.text      = element_text(face = "bold", size = rel(0.95)),
          legend.position = "bottom", panel.spacing = unit(4, "pt"),
          plot.background = element_rect(fill = "white", colour = NA))
}
# coord_sf (raster CRS) when outlines drawn so tiles + polygons share one space.
add_geo <- function(p) {
  if (!is.null(wmd_sf)) {
    p + geom_sf(data = wmd_sf, fill = NA, colour = WMD_LINE_COL,
                linewidth = WMD_LINE_W, inherit.aes = FALSE) +   # ggplot2 < 3.4: size=
      coord_sf(crs = map_crs, datum = NA, expand = FALSE)
  } else p + coord_equal()
}

# ═══════════════════════════════════════════════════════════════════════════════
# DESCRIPTIVE SUMMARY  —  the three per-EAU model inputs, before any optimisation
#   fig_benefit_ribbon        benefit over time, by scenario (median / IQR / 10–90)
#   map_benefit_2020_2100     benefit in space, 2020 vs 2100, one shared colour scale
#   map_cost_binned           cost in space, 2020 $/ha, equal-count classes
#   fig_risk_box              risk distribution, linear axis, by RCP
#
# Reads the inputs panel (script 06 output), not the model results, so this section
# is independent of 07/08 and renders even before any model has been run.
#
# STRUCTURAL FACTS THAT SHAPE THESE DISPLAYS:
#  * benefit varies by eau × year × rcp × gcm. In the ribbon, quantiles are taken
#    ACROSS EAUs within a scenario-year, so the bands show how unevenly benefit is
#    spread over the landscape, NOT projection or GCM uncertainty; the legend says so.
#  * cost varies by eau × year only (scenario-invariant), and its time dimension is a
#    deterministic 1.02^(t−2017) rescale of one 2017 base, so space is the only axis
#    along which it meaningfully varies. Per script 06, cost = mean_fmv_per_ha ×
#    eau_area_ha with eau_area_ha a SCALAR (EAUs are equal-area), so $/ha is an exact
#    constant rescale of total cost: the map is identical either way, and $/ha is used
#    because it is the directly quotable unit.
#  * risk varies by eau × year × rcp (NOT wet/dry), spans ~5 orders of magnitude, is
#    exactly 0 at 2100 by construction (excluded), and has a large mass at the ε floor.
# ═══════════════════════════════════════════════════════════════════════════════
drew_desc <- FALSE

if (file.exists(PANEL_RDS)) {

  panel <- readRDS(PANEL_RDS)

  ## palettes / labels used only by this section --------------------------------
  RCP_COL <- c(`RCP 4.5` = "#3F7FA6", `RCP 8.5` = "#9C3848")
  LAB_MED <- "Median across EAUs"
  LAB_IQR <- "Interquartile range (25th\u201375th)"
  LAB_B90 <- "10th\u201390th percentile"
  # EAU area derives from CELL_W x CELL_H so there is a single source of truth for cell size
  EAU_AREA_HA <- CELL_W * CELL_H / 1e4
  # captions + legend titles are off by default in theme_map(); the descriptive maps
  # need both, added here rather than in theme_map() so Figures 3–4 are untouched
  DESC_MAP_EXTRA <- theme(
    legend.title = element_text(colour = "grey25", size = rel(0.9)),
    legend.text  = element_text(size = rel(0.8)),
    plot.caption = element_text(colour = "grey40", size = rel(0.72), hjust = 0))

  scen_of <- function(rcp, gcm) dplyr::case_when(
    rcp == "stationary" ~ "stationary",
    rcp == "baseline"   ~ "baseline",
    TRUE                ~ paste0("rcp", rcp, "_", gcm))

  ## benefit — the 2020 row is a SHARED cross-scenario baseline (rcp = "baseline"),
  ## exactly as script 07 treats period 0. Replicated onto all five scenarios so every
  ## trajectory starts from the common anchor; all five facets are therefore identical
  ## at 2020 BY CONSTRUCTION, not by coincidence.
  benefit <- bind_rows(
    panel %>% filter(year == 2020, rcp == "baseline") %>%
      transmute(eau_id, year, value = scaled_abundance) %>%
      tidyr::crossing(scenario = SCEN_LEVELS),
    panel %>% filter(year > 2020) %>%
      transmute(eau_id, year, scenario = scen_of(rcp, gcm), value = scaled_abundance)
  ) %>% mutate(scenario = factor(scenario, levels = SCEN_LEVELS))

  ## cost — keys are eau × year only
  cost <- panel %>% distinct(eau_id, year, value = cost)

  ## risk — keys are eau × year × rcp. Script 05 assigns trans_prob = epsilon under
  ## stationary EXCEPT where prop_suitable == 0, which gets exactly 0; a plain min()
  ## would therefore return 0 rather than epsilon and silently break the floor
  ## accounting, so take the smallest POSITIVE value instead.
  stat_vals <- panel %>% filter(rcp == "stationary", year < 2100) %>% pull(trans_prob)
  EPS <- min(stat_vals[stat_vals > 0])
  risk <- panel %>%
    distinct(eau_id, year, rcp, value = trans_prob) %>%
    filter(year < 2100, rcp %in% c("45", "85")) %>%    # 2100 is structurally zero
    mutate(rcp_lab = factor(recode(rcp, `45` = "RCP 4.5", `85` = "RCP 8.5"),
                            levels = c("RCP 4.5", "RCP 8.5")),
           at_floor = value <= EPS * 1.000001,          # includes structural zeros
           pct      = value * 100)

  # ─────────────────────────────────────────────────────────────────────────────
  # FIGURE — benefit over time: median + IQR + 10–90 ribbons, faceted by scenario
  # The statistical definitions live in a LEGEND, not the y-axis title. colour/fill
  # are already spoken for by `scenario` (and suppressed, since the facet strips carry
  # it), so the bands map to a dummy `alpha` aesthetic and the central line to a dummy
  # `linetype`; override.aes then draws the keys in neutral grey so they read as
  # statistics rather than as a sixth series. Key alphas match the panels, so the
  # darker key is visibly the IQR and the lighter one the 10–90 band.
  # ─────────────────────────────────────────────────────────────────────────────
  ben_band <- benefit %>% group_by(scenario, year) %>%
    summarise(p10 = quantile(value, .10), p25 = quantile(value, .25),
              med = median(value),
              p75 = quantile(value, .75), p90 = quantile(value, .90), .groups = "drop")

  fig_benefit_ribbon <- ggplot(ben_band, aes(year, med, colour = scenario, fill = scenario)) +
    geom_ribbon(aes(ymin = p10, ymax = p90, alpha = LAB_B90), colour = NA) +
    geom_ribbon(aes(ymin = p25, ymax = p75, alpha = LAB_IQR), colour = NA) +
    geom_line(aes(linetype = LAB_MED), linewidth = .8) +
    facet_wrap(~ scenario, ncol = 5, labeller = as_labeller(SCEN_LAB1)) +
    scale_colour_manual(values = SCEN_COL) +
    scale_fill_manual(values = SCEN_COL) +
    scale_alpha_manual(name = NULL, values = setNames(c(.30, .16), c(LAB_IQR, LAB_B90)),
                       breaks = c(LAB_IQR, LAB_B90)) +
    scale_linetype_manual(name = NULL, values = setNames("solid", LAB_MED)) +
    scale_x_continuous("year", breaks = seq(2020, 2100, 40)) +
    scale_y_continuous("Duck pairs per EAU", labels = label_comma()) +
    guides(colour = "none", fill = "none",
           linetype = guide_legend(order = 1,
                                   override.aes = list(colour = "grey25", fill = NA)),
           alpha    = guide_legend(order = 2,
                                   override.aes = list(fill = "grey25", colour = NA))) +
    theme(legend.position = "top", legend.box = "horizontal",
          legend.key.height = unit(10, "pt"), legend.text = element_text(size = rel(.8)))
  save_fig(fig_benefit_ribbon, "fig_benefit_ribbon", 13.5, 4.4)

  # ─────────────────────────────────────────────────────────────────────────────
  # FIGURE — conversion risk on a LINEAR axis, faceted by RCP
  # A log axis makes a right-skewed distribution look symmetric, which is what a log
  # axis is for: it shows the bulk at the cost of hiding the tail's severity. On a
  # linear axis the box collapses toward zero and the tail sprays upward — that IS the
  # distribution, shown rather than transformed away. Two annotations make the
  # collapsed box interpretable:
  #  * the number above each box is the % of EAUs at the background floor. Without it a
  #    box flattened against the axis cannot be told from "low risk everywhere". WATCH:
  #    if that share exceeds 75% the IQR is exactly zero, the box degenerates to a line,
  #    and whiskers have zero length (1.5 × 0), so every non-floor EAU is an outlier.
  #  * the horizontal tick is the MAXIMUM. Every EAU is already drawn by the jitter
  #    layer (outlier.shape = NA only suppresses the boxplot's duplicate rendering of
  #    the outliers, it does not drop them), but the alpha needed to keep ~1,000
  #    overlapping floor points legible renders a lone extreme point nearly invisible.
  #    The MINIMUM is deliberately unmarked: it is pinned at the floor every decade, so
  #    it is a structural constant with hundreds of ties rather than an observation.
  # ─────────────────────────────────────────────────────────────────────────────
  LAB_Y <- max(risk$pct, na.rm = TRUE) * 1.10
  floor_lab <- risk %>% group_by(rcp_lab, year) %>%
    summarise(pct_floor = 100 * mean(at_floor), .groups = "drop") %>%
    mutate(lab = sprintf("%.0f", pct_floor), lab_y = LAB_Y)

  fig_risk_box <- ggplot(risk, aes(factor(year), pct, fill = rcp_lab, colour = rcp_lab)) +
    geom_jitter(width = .18, size = .3, alpha = .10, stroke = 0) +
    geom_boxplot(outlier.shape = NA, alpha = .45, linewidth = .32, width = .62) +
    # explicit maximum: a zero-height errorbar sits on the exact value, unlike a glyph
    # shape whose baseline offset would place it slightly off
    stat_summary(fun.min = max, fun.max = max, geom = "errorbar",
                 width = .34, linewidth = .45, colour = "grey15") +
    stat_summary(aes(shape = "Mean"), fun = mean, geom = "point",
                 size = 1.7, fill = "white", colour = "grey15", stroke = .4) +
    geom_text(data = floor_lab, aes(factor(year), lab_y, label = lab),
              size = 2.6, vjust = 0, show.legend = FALSE) +
    facet_wrap(~ rcp_lab, ncol = 2) +
    scale_fill_manual(values = RCP_COL, guide = "none") +
    scale_colour_manual(values = RCP_COL, guide = "none") +
    scale_shape_manual(name = NULL, values = c(Mean = 23)) +
    scale_y_continuous("Per-period conversion risk",
                       labels = label_number(suffix = "%", accuracy = .1),
                       expand = expansion(mult = c(.02, .12))) +
    scale_x_discrete("year") +
    labs(caption = paste0(
      "Box = median and interquartile range; whiskers = 1.5 \u00d7 IQR; horizontal tick = ",
      "maximum; points = individual EAUs. Numbers above each box give the % of EAUs at ",
      "the ", EPS_LAB, ".")) +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = .5, size = rel(.8)),
          plot.caption = element_text(colour = "grey40", size = rel(.72), hjust = 0))
  save_fig(fig_risk_box, "fig_risk_box", 11, 5.0)

  # ── maps need EAU centroids from the script-01 lookup ────────────────────────
  if (!is.null(eau_xy)) {

    # shared plate skeleton: `d` needs x_coord, y_coord, fill_val. Callers add the fill
    # scale and caption, so benefit and cost stay identical in every other respect.
    eau_map <- function(d) add_geo(
      ggplot(d, aes(x_coord, y_coord, fill = fill_val)) +
        geom_tile(width = CELL_W, height = CELL_H))
    MAP_CAP <- paste0("One tile = one EAU (", format(round(EAU_AREA_HA), big.mark = ","),
                      " ha); outlines are Wetland Management Districts.")

    # ───────────────────────────────────────────────────────────────────────────
    # MAP — benefit in space: 2020 baseline vs 2100, ONE SHARED colour scale
    # The shared scale is the point: separately-normalised panels would each show
    # their own internal gradient and the decline between them would vanish. Panel
    # labels are asymmetric on purpose — 2020 is the shared cross-scenario anchor
    # (the surface under all five facets of fig_benefit_ribbon), 2100 is one
    # scenario's projection; bare years would imply a matched pair.
    # ───────────────────────────────────────────────────────────────────────────
    BEN_MAP_LABS <- c(paste0(BEN_MAP_YEARS[1], " (shared baseline)"),
                      paste0(BEN_MAP_YEARS[2], " (", SCEN_LAB1[[BEN_MAP_SCEN]], ")"))
    ben_map_dat <- benefit %>%
      filter(scenario == BEN_MAP_SCEN, year %in% BEN_MAP_YEARS) %>%
      distinct(eau_id, year, fill_val = value) %>%
      inner_join(eau_xy, by = "eau_id") %>%
      mutate(panel = factor(year, levels = BEN_MAP_YEARS, labels = BEN_MAP_LABS))

    ben_lim <- if (BEN_CAP_Q >= 1) NULL else
      c(0, as.numeric(quantile(ben_map_dat$fill_val, BEN_CAP_Q, na.rm = TRUE)))

    map_benefit <- eau_map(ben_map_dat) +
      facet_wrap(~ panel, ncol = 2) +
      scale_fill_viridis_c(option = BEN_RAMP, direction = -1, begin = .08, end = .96,
                           labels = label_comma(), na.value = MAP_NA_COL,
                           limits = ben_lim, oob = scales::squish,
                           name = "Duck pairs per EAU",
                           guide = guide_colourbar(title.position = "top",
                                                   barwidth = 16, barheight = .6)) +
      theme_map() + DESC_MAP_EXTRA +
      labs(caption = paste0(
        MAP_CAP, " Both panels share one colour scale, so tiles are directly comparable",
        " across years. Colour is strictly proportional to duck pairs (no transform).",
        if (BEN_CAP_Q < 1)
          sprintf("\nRamp capped at the %.0fth percentile; EAUs above it share the top colour.",
                  BEN_CAP_Q * 100) else ""))
    save_fig(map_benefit, "map_benefit_2020_2100", 13, 5.4)

    # ───────────────────────────────────────────────────────────────────────────
    # MAP — cost in space, 2020 $/ha, equal-count classes
    # Quantile classes rather than a proportional ramp: cost is heavily right-skewed,
    # so a continuous scale puts almost every EAU at the low end and leaves the
    # interior gradient unreadable. Each class holds ~1/N_COST_BINS of the EAUs, so
    # contrast is spread across the map; break labels are raw 2020 dollars, so the
    # values stay legible. What is given up is proportionality of colour distance.
    # ───────────────────────────────────────────────────────────────────────────
    cost_map_dat <- cost %>% filter(year == MAP_YEAR) %>%
      inner_join(eau_xy, by = "eau_id") %>%
      mutate(fill_val = value / EAU_AREA_HA)        # $/ha (constant rescale of total)

    fmt_ha <- label_dollar(accuracy = 1)
    qb <- unique(quantile(cost_map_dat$fill_val,
                          probs = seq(0, 1, length.out = N_COST_BINS + 1), na.rm = TRUE))
    cost_map_dat$fill_val <- cut(cost_map_dat$fill_val, breaks = qb, include.lowest = TRUE)
    levels(cost_map_dat$fill_val) <- sprintf("%s \u2013 %s",
                                             fmt_ha(head(qb, -1)), fmt_ha(qb[-1]))

    map_cost <- eau_map(cost_map_dat) +
      scale_fill_viridis_d(option = COST_RAMP, direction = -1, begin = .08, end = .96,
                           na.value = MAP_NA_COL,
                           name = paste0("Acquisition cost per hectare, ", MAP_YEAR,
                                         " USD  (", nlevels(cost_map_dat$fill_val),
                                         " equal-count classes)"),
                           guide = guide_legend(title.position = "top", nrow = 1,
                                                override.aes = list(colour = NA),
                                                keywidth = unit(20, "pt"),
                                                keyheight = unit(9, "pt"))) +
      theme_map() + DESC_MAP_EXTRA +
      labs(caption = paste0(
        MAP_CAP, " Raw 2020 dollars, no log transform. EAUs are equal-area, so",
        " per-hectare cost is a constant rescale of total acquisition cost \u2014 the",
        " spatial pattern is identical either way."))
    save_fig(map_cost, "map_cost_binned", 10, 6.4)

  } else {
    message("EAU lookup '", EAU_LOOKUP, "' not found \u2014 skipping descriptive maps. ",
            "Run 01_create_EAUs.R, then re-run.")
  }

  drew_desc <- TRUE

} else {
  message("Inputs panel '", PANEL_RDS, "' not found \u2014 skipping descriptive assets. ",
          "Run 06_cost_data.R, then re-run.")
}

# ═══════════════════════════════════════════════════════════════════════════════
# MAPS (→ Figures 3–4)  —  acquisition footprints from 08's persisted schedule
# Reads acquisition_schedule_spatial.csv (scenario, model, eau_id, wmd_id, x/y,
# acquired_period/year, prevented_loss). Cells are ~342 km^2 EAUs in USGS Conus Albers
# (ESRI:102039); the optional WMD overlay dissolves wmd_raster_equal_area.tif to outlines.
# ═══════════════════════════════════════════════════════════════════════════════
if (DRAW_MAPS && file.exists(SCHEDULE)) {

  sched <- read_csv(SCHEDULE, show_col_types = FALSE)
  land  <- sched %>% distinct(eau_id, x_coord, y_coord)   # landscape backdrop
  acq   <- function(df) df %>% filter(!is.na(acquired_year))

  ## MAP (→ Figure 3) — foresight difference (rolling vs myopic), per climate scenario
  classify_pair <- function(sc) {
    r <- acq(sched) %>% filter(scenario == sc, model == "rolling") %>% pull(eau_id)
    m <- acq(sched) %>% filter(scenario == sc, model == "myopic")  %>% pull(eau_id)
    land %>% mutate(scenario = sc,
      cat = factor(case_when(
        eau_id %in% r & eau_id %in% m ~ "Acquired by both",
        eau_id %in% r                 ~ "Rolling only (foresight adds)",
        eau_id %in% m                 ~ "Myopic only (foresight avoids)",
        TRUE                          ~ "Neither"), levels = names(PAIR_COL)))
  }
  md_fore <- bind_rows(lapply(CLIM, classify_pair)) %>%
    mutate(scenario = factor(SCEN_LAB1[scenario], levels = SCEN_LAB1[CLIM])) %>%
    arrange(scenario, cat %in% c("Rolling only (foresight adds)",
                                 "Myopic only (foresight avoids)"))
  mp_fore <- ggplot(md_fore, aes(x_coord, y_coord, fill = cat)) +
    geom_tile(width = CELL_W, height = CELL_H) +
    scale_fill_manual(values = PAIR_COL, name = NULL) +
    facet_wrap(~ scenario, ncol = 2) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(colour = NA)))
  mp_fore <- add_geo(mp_fore) + theme_map()
  save_fig(mp_fore, "foresight_difference", 13, 11)

  ## MAP (→ Figure 4) — per-parcel prevented-loss difference (ΔV = rolling − myopic)
  # Each cell's signed value is the true-future prevented loss rolling banks there minus
  # myopic's (0 if that policy never buys the parcel). Blue = rolling banks more; red =
  # myopic banks more (allocative, not "worthless"); grey = no difference. Signed cells
  # sum to the value_added gap. Heavy-tailed, so colour is capped symmetrically at the
  # DIFF_CAP_Q quantile of |ΔV| (cells beyond squished to full saturation).
  md_dv <- sched %>%
    filter(scenario %in% CLIM, model %in% c("rolling", "myopic")) %>%
    select(scenario, eau_id, x_coord, y_coord, model, prevented_loss) %>%
    pivot_wider(names_from = model, values_from = prevented_loss) %>%
    mutate(rolling = coalesce(rolling, 0), myopic = coalesce(myopic, 0),
           dV = rolling - myopic,
           scenario = factor(SCEN_LAB1[scenario], levels = SCEN_LAB1[CLIM]))
  nz   <- md_dv$dV[md_dv$dV != 0]
  dlim <- if (length(nz)) as.numeric(quantile(abs(nz), DIFF_CAP_Q, na.rm = TRUE)) else 1
  if (!is.finite(dlim) || dlim <= 0) dlim <- max(abs(md_dv$dV), 1)
  mp_dv <- ggplot() +
    geom_tile(data = land, aes(x_coord, y_coord), width = CELL_W, height = CELL_H, fill = MAP_BG) +
    geom_tile(data = dplyr::filter(md_dv, dV != 0),
              aes(x_coord, y_coord, fill = dV), width = CELL_W, height = CELL_H) +
    scale_fill_gradient2(low = DIFF_NEG, mid = "white", high = DIFF_POS, midpoint = 0,
                         limits = c(-dlim, dlim), oob = scales::squish,
                         labels = scales::label_comma(),
                         name = "\u0394 prevented loss\n(rolling \u2212 myopic, pairs)",
                         guide = guide_colourbar(barwidth = 14)) +
    facet_wrap(~ scenario, ncol = 2)
  mp_dv <- add_geo(mp_dv) + theme_map()
  save_fig(mp_dv, "prevented_loss_diff", 13, 11)

} else if (DRAW_MAPS) {
  message("Schedule '", SCHEDULE, "' not found \u2014 skipping maps. ",
          "Run an 08 that persists the schedule, then re-run. (Tables + figures were written.)")
}

## run metadata (results.qmd reads this to stamp the document) -------------------
drew_maps <- DRAW_MAPS && file.exists(SCHEDULE)
HAZARD_K <- if (HAZARD_RUN)   # high-hazard run: k solved by 14_build_hazard_panel.R
  jsonlite::read_json(file.path(DIR_DERIVED, "hazard_scaling.json"))$k else 1
jsonlite::write_json(
  list(generated      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
       generated_date = format(Sys.Date(), "%Y-%m-%d"),
       maps_drawn     = drew_maps,
       run_id         = RUN_ID,
       derived_run_id = DERIVED_RUN_ID,
       budget_eaus_per_period = BUDGET_EAUS_PER_PERIOD,
       hazard_target_annual_pct = HAZARD_TARGET_ANNUAL_PCT,
       hazard_k       = HAZARD_K,
       lc_scenario    = LC_SCENARIO,
       lc_label       = LC_LABEL,
       r_version      = R.version.string),
  file.path(OUT_DIR, "run_metadata.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")

message("\nDone. Wrote 2 tables + 2 figures",
        if (drew_maps) " + 2 maps" else "",
        if (drew_desc) " + 4 descriptive assets" else "",
        " + run_metadata.json to '", OUT_DIR, "/'.")

# ═══════════════════════════════════════════════════════════════════════════════
# COMPARISON SECTIONS (formerly 10b, 10c, 10d)
# ═══════════════════════════════════════════════════════════════════════════════
# Each section below is the former stand-alone script, unchanged except that a
# missing input or a non-canonical configuration skips the section (with a message)
# instead of stopping. Real errors still stop. Each section runs in its own
# environment: several of its settings and helpers share names with the ones above
# (OUT_DIR, theme_ppr, CLIM, ...) but differ, so they must not overwrite each other.

skip_section <- function(...)
  stop(structure(class = c("section_skip", "error", "condition"),
                 list(message = paste0(...), call = NULL)))

run_section <- function(name, code) {
  message("\n── ", name, " ", strrep("─", max(3, 66 - nchar(name))))
  tryCatch(eval(substitute(code), envir = new.env(parent = globalenv())),
           section_skip = function(e) message("Skipped ", name, ": ", conditionMessage(e)))
  invisible(NULL)
}

run_section("Land-cover comparison", {
# ═══════════════════════════════════════════════════════════════════════════════
# SECTION: land-cover comparison (formerly 10b) — assets comparing two land-cover runs
# ═══════════════════════════════════════════════════════════════════════════════
# Compares this run (the BASE: RUN_ID in 00_config.R, normally "gcam_ref") with an
# ALTERNATIVE land-cover run (COMPARE_RUN below, normally "bau"). Differences are
# reported as alternative − base, with the same formulas as 12_compare_runs.qmd.
#
# PREREQUISITES — both complete pipeline runs must exist on disk:
#   runs/<RUN_ID>/output_data, runs/<RUN_ID>/derived_data       (base)
#   runs/<COMPARE_RUN>/output_data, runs/<COMPARE_RUN>/derived_data  (alternative)
# Typical order: run the pipeline for the alternative scenario, then for the base
# scenario; this section then runs at the end of the base run's script 10 (see README).
#
# Writes to the BASE run's output_figs (raw numbers; formatting lives in
# 11_results.qmd):
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

if (HAZARD_RUN || BUDGET_SWEEP)
  skip_section("runs only in the canonical configuration (PPR_HAZARD and PPR_BUDGET unset).")
if (COMPARE_RUN == RUN_ID)
  skip_section("this run is the comparison run (COMPARE_RUN = '", COMPARE_RUN, "').")
for (f in c(BASE_RES, ALT_RES))
  if (!file.exists(f))
    skip_section("missing '", f, "'. Both land-cover runs must be complete ",
                 "(see the header of this section).")
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
})

run_section("Budget sensitivity", {
# ═══════════════════════════════════════════════════════════════════════════════
# SECTION: budget sensitivity (formerly 10c) — assets comparing budget-sensitivity runs
# ═══════════════════════════════════════════════════════════════════════════════
# Compares the budget runs made by 13_run_budget_sweep.R: every complete folder
# runs/<LC_SCENARIO>_bXX/ (e.g. gcam_ref_b02, gcam_ref_b05, gcam_ref_b10). The budget
# of each run is read from its own model_results.rds (params$budget_eaus_per_period),
# not from the folder name. Formulas are the same as in 10_results_assets.R.
#
# PREREQUISITES:
#   - the canonical run, runs/<LC_SCENARIO>/ (Rscript 00_run_all.R, PPR_BUDGET unset)
#   - at least two complete budget runs (Rscript 13_run_budget_sweep.R)
# Runs as part of script 10 in the canonical configuration (PPR_BUDGET unset):
#   Rscript 10_results_assets.R
#
# Writes to the CANONICAL run's output_figs (raw numbers; formatting lives in the
# budget-comparison .qmd), all prefixed bud_:
#   bud_value_added.csv     value added, % over baseline, EAUs and hectares bought,
#                           by budget x scenario x policy
#   bud_headline.csv        Tables 1-2 by budget: rolling vs greedy/myopic in pairs,
#                           % improvement and % shortfall, + ranking check
#   bud_convergence.csv     solver convergence by budget x scenario x policy
#   fig_bud_improvement     % improvement of rolling vs budget (PNG + PDF)
#   fig_bud_value_added     acquisition value vs budget, by policy (PNG + PDF)
#   fig_bud_overtime        Figure 5 by budget: cumulative gap vs rolling (PNG + PDF)
#   map_bud_foresight       Figure 6 by budget: rolling vs myopic acquisition sets
#   bud_metadata.json       runs, budgets, check of the canonical-budget run
# No dollar amounts are written (total_spend is deliberately not carried over).
# ═══════════════════════════════════════════════════════════════════════════════

if (!exists(".SETUP_DONE")) source("00_setup.R")
source("00_config.R")   # always re-read: a long-lived R session can hold an older config
if (HAZARD_RUN || BUDGET_SWEEP)
  skip_section("runs only in the canonical configuration (PPR_HAZARD and PPR_BUDGET unset).")

# ── CONFIG ─────────────────────────────────────────────────────────────────────
SWEEP_PATTERN <- paste0("^", LC_SCENARIO, "_b[0-9]+(p[0-9]+)?$")   # budget run folders
OUT_DIR       <- DIR_FIGS                                         # canonical output_figs
CANON_RES     <- file.path(DIR_OUT, "model_results.csv")          # canonical results
EAU_LOOKUP    <- file.path(DIR_COMMON, "eau_wmd_lookup.rds")

DPI         <- 300
BASE_SIZE   <- 12
BASE_FAMILY <- "serif"          # same as 10_results_assets.R
SAVE_PDF    <- TRUE
DRAW_MAPS   <- TRUE             # FALSE skips map_bud_foresight

# palettes and map settings: same as 10_results_assets.R
COL_MODEL <- c(greedy = "#B0413E", myopic = "#E2A23B", rolling = "#2E7D5B")
SCEN_COL  <- c(rcp45_wet = "#3F7FA6", rcp45_dry = "#9BB7C4",
               rcp85_wet = "#9C3848", rcp85_dry = "#D79A86", stationary = "#777777")
CELL_W       <- 18513.54        # EAU width  (x resolution), metres
CELL_H       <- 18485.21        # EAU height (y resolution), metres
EAU_AREA_HA  <- CELL_W * CELL_H / 1e4
MAP_BG       <- "grey91"
WMD_OVERLAY  <- TRUE
WMD_RASTER   <- file.path(DIR_COMMON, "wmd_raster_equal_area.tif")
WMD_LINE_COL <- "grey25"
WMD_LINE_W   <- 0.25
PAIR_COL <- c("Acquired by both"               = "#5E4FA2",
              "Rolling only (foresight adds)"  = "#2C7FB8",
              "Myopic only (foresight avoids)" = "#E6892B",
              "Neither"                        = MAP_BG)
# ────────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(scales)
})
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

## labels, helpers (same definitions as 10_results_assets.R) --------------------
SCEN_LEVELS <- c("rcp45_wet", "rcp45_dry", "rcp85_wet", "rcp85_dry", "stationary")
SCEN_LAB1   <- c(rcp45_wet = "RCP 4.5 wet", rcp45_dry = "RCP 4.5 dry",
                 rcp85_wet = "RCP 8.5 wet", rcp85_dry = "RCP 8.5 dry",
                 stationary = "Stationary")
CLIM <- SCEN_LEVELS[1:4]
pct_improve <- function(roll, base) ifelse(abs(base) < 1e-9, NA_real_, 100 * (roll - base) / base)
pct_short   <- function(roll, base) ifelse(abs(roll) < 1e-9, NA_real_, 100 * (roll - base) / roll)
fmt_pairs   <- function(x) ifelse(abs(x) < 1, "<1", formatC(round(x), format = "d", big.mark = ","))
bud_lab     <- function(b) paste0(format(b), " EAUs / period")

theme_ppr <- function(base_size = BASE_SIZE, base_family = BASE_FAMILY) {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.title    = element_text(colour = "grey25"),
      axis.title.y  = element_text(angle = 90, margin = margin(r = 6)),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "top", legend.justification = "right",
      legend.title = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0, size = rel(0.95),
                                margin = margin(2, 0, 4, 0)),   # stops descenders clipping
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(10, 14, 8, 10)
    )
}
theme_set(theme_ppr())
ggplot2::update_geom_defaults("text", list(family = BASE_FAMILY))

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
# FIND THE BUDGET RUNS
# ═══════════════════════════════════════════════════════════════════════════════
cand <- grep(SWEEP_PATTERN, list.dirs("runs", recursive = FALSE, full.names = FALSE),
             value = TRUE)
runs <- bind_rows(lapply(cand, function(id) {
  od  <- file.path("runs", id, "output_data")
  rds <- file.path(od, "model_results.rds")
  if (!file.exists(rds) || !file.exists(file.path(od, "model_results.csv"))) {
    message("Skipping ", id, ": no complete model results.")
    return(NULL)
  }
  b <- readRDS(rds)$params$budget_eaus_per_period
  if (is.null(b)) { message("Skipping ", id, ": budget not recorded in params."); return(NULL) }
  data.frame(run_id = id, budget = as.numeric(b), dir = od, stringsAsFactors = FALSE)
}))
if (nrow(runs) < 2)
  skip_section("found ", nrow(runs), " complete budget run(s) matching runs/", LC_SCENARIO,
               "_bXX. Run 13_run_budget_sweep.R first.")
if (anyDuplicated(runs$budget))
  stop("Two runs share a budget: ", paste(runs$run_id[duplicated(runs$budget) |
       duplicated(runs$budget, fromLast = TRUE)], collapse = ", "), call. = FALSE)
runs <- arrange(runs, budget)
BUD_LEVELS <- bud_lab(runs$budget)
message("Budget runs: ", paste0(runs$run_id, " (", runs$budget, ")", collapse = ", "))

read_runs <- function(file) bind_rows(lapply(seq_len(nrow(runs)), function(k) {
  p <- file.path(runs$dir[k], file)
  if (!file.exists(p)) { message("  (", runs$run_id[k], ": no ", file, ")"); return(NULL) }
  read_csv(p, show_col_types = FALSE) %>%
    mutate(budget = runs$budget[k], run_id = runs$run_id[k], .before = 1)
}))
res  <- read_runs("model_results.csv")
trj  <- read_runs("model_trajectories.csv")
conv <- read_runs("solver_convergence.csv")

# ═══════════════════════════════════════════════════════════════════════════════
# TABLES — tidy CSVs of RAW NUMBERS
# ═══════════════════════════════════════════════════════════════════════════════

# value added, % over baseline, acquisitions -- long, one row per budget x scenario x policy
tbl_value <- res %>%
  transmute(budget, scenario = factor(scenario, levels = SCEN_LEVELS), model,
            value_added,
            pct_over_baseline = 100 * value_added / J_baseline,
            va_per_budget_unit = value_added / budget,   # pairs per median-EAU of budget
            n_acquired,
            hectares = n_acquired * EAU_AREA_HA) %>%
  arrange(budget, scenario, model)
write_csv(tbl_value, file.path(OUT_DIR, "bud_value_added.csv"))

# Tables 1-2 by budget ("<1 pair" differences are not meaningful as %)
tbl_headline <- res %>%
  select(budget, scenario, model, value_added) %>%
  pivot_wider(names_from = model, values_from = value_added) %>%
  transmute(budget, scenario = factor(scenario, levels = SCEN_LEVELS),
            greedy_va = greedy, myopic_va = myopic, rolling_va = rolling,
            rg_dpairs = rolling - greedy, rm_dpairs = rolling - myopic,
            rg_pct       = ifelse(abs(rg_dpairs) < 1, NA_real_, pct_improve(rolling, greedy)),
            rm_pct       = ifelse(abs(rm_dpairs) < 1, NA_real_, pct_improve(rolling, myopic)),
            rg_short_pct = ifelse(abs(rg_dpairs) < 1, NA_real_, pct_short(rolling, greedy)),
            rm_short_pct = ifelse(abs(rm_dpairs) < 1, NA_real_, pct_short(rolling, myopic)),
            ranking_holds = rolling >= myopic - 1e-6 & myopic >= greedy - 1e-6) %>%
  arrange(budget, scenario)
write_csv(tbl_headline, file.path(OUT_DIR, "bud_headline.csv"))

# solver convergence (tag = "<scenario>/<model>"; greedy uses no solver)
if (!is.null(conv) && nrow(conv)) {
  tbl_conv <- conv %>%
    separate(tag, into = c("scenario", "model"), sep = "/", remove = TRUE) %>%
    select(budget, scenario, model, n_solves, n_hit_timelimit, worst_gap_pct) %>%
    mutate(scenario = factor(scenario, levels = SCEN_LEVELS)) %>%
    arrange(budget, scenario, model)
  write_csv(tbl_conv, file.path(OUT_DIR, "bud_convergence.csv"))
}

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE — % improvement of rolling vs budget (climate scenarios)
# ═══════════════════════════════════════════════════════════════════════════════
d_imp <- tbl_headline %>%
  filter(scenario %in% CLIM) %>%
  select(budget, scenario, rg_pct, rm_pct) %>%
  pivot_longer(c(rg_pct, rm_pct), names_to = "cmp", values_to = "pct") %>%
  mutate(cmp = factor(cmp, levels = c("rg_pct", "rm_pct"),
                      labels = c("Rolling vs. greedy", "Rolling vs. myopic")),
         scenario = factor(scenario, levels = CLIM))
fig_imp <- ggplot(d_imp, aes(budget, pct, colour = scenario)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  facet_wrap(~ cmp, nrow = 1) +
  scale_colour_manual(values = SCEN_COL, labels = SCEN_LAB1) +
  scale_x_continuous("Budget (median-priced EAUs per period)", breaks = runs$budget) +
  scale_y_continuous("% improvement of rolling", labels = function(x) paste0(x, "%"),
                     limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  theme(panel.grid.major.x = element_line(colour = "grey92"))
save_fig(fig_imp, "fig_bud_improvement", 10, 4.8)

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE — acquisition value vs budget, by policy (climate scenarios)
# ═══════════════════════════════════════════════════════════════════════════════
d_va <- tbl_value %>%
  filter(scenario %in% CLIM) %>%
  mutate(model = factor(model, levels = c("rolling", "myopic", "greedy")),
         scenario = factor(scenario, levels = CLIM))
fig_va <- ggplot(d_va, aes(budget, value_added, colour = model)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  facet_wrap(~ scenario, ncol = 2, labeller = as_labeller(SCEN_LAB1)) +
  scale_colour_manual(values = COL_MODEL, breaks = c("rolling", "myopic", "greedy"),
                      labels = c("Rolling", "Myopic", "Greedy")) +
  scale_x_continuous("Budget (median-priced EAUs per period)", breaks = runs$budget) +
  scale_y_continuous("Acquisition value\n(discounted duck-pairs over do-nothing)",
                     labels = label_comma(), limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  theme(panel.grid.major.x = element_line(colour = "grey92"))
save_fig(fig_va, "fig_bud_value_added", 10, 7.4)

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE (Figure 5 by budget) — cumulative performance relative to rolling
# Rows = climate scenario, columns = budget. Same calculation as overtime_gap in 10.
# ═══════════════════════════════════════════════════════════════════════════════
roll_disc <- trj %>% filter(model == "rolling") %>%
  select(budget, scenario, year, roll = discounted)
cumg <- trj %>%
  filter(model %in% c("greedy", "myopic", "rolling"), scenario %in% CLIM) %>%
  left_join(roll_disc, by = c("budget", "scenario", "year")) %>%
  arrange(budget, scenario, model, year) %>%
  group_by(budget, scenario, model) %>%
  mutate(cum_gap = cumsum(discounted - roll)) %>%     # rolling: exactly 0
  ungroup() %>%
  mutate(model    = factor(model, levels = c("rolling", "myopic", "greedy")),
         scenario = factor(scenario, levels = CLIM),
         bud      = factor(bud_lab(budget), levels = BUD_LEVELS))
ends <- cumg %>% filter(model != "rolling") %>%
  group_by(budget, scenario, model) %>% filter(year == max(year)) %>% ungroup()

fig_overtime <- ggplot(cumg, aes(year, cum_gap, colour = model)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
  geom_text(data = ends, aes(label = fmt_pairs(cum_gap)), hjust = -0.12, size = 2.5,
            show.legend = FALSE) +
  facet_grid(scenario ~ bud, labeller = labeller(scenario = SCEN_LAB1)) +
  scale_colour_manual(values = COL_MODEL, breaks = c("rolling", "myopic", "greedy"),
                      labels = c("Rolling (benchmark)", "Myopic", "Greedy")) +
  scale_x_continuous("year", breaks = seq(2020, 2100, 40),
                     expand = expansion(mult = c(0.02, 0.30))) +
  scale_y_continuous("Cumulative difference vs. rolling\n(discounted duck-pairs)",
                     labels = label_comma()) +
  theme(panel.grid.major.x = element_line(colour = "grey92"),
        strip.text.y = element_text(face = "bold", angle = -90, size = rel(0.9)))
save_fig(fig_overtime, "fig_bud_overtime", 11, 11)

# ═══════════════════════════════════════════════════════════════════════════════
# MAP (Figure 6 by budget) — rolling vs myopic acquisition sets
# Rows = climate scenario, columns = budget. Same classification as 10.
# ═══════════════════════════════════════════════════════════════════════════════
drew_maps <- FALSE
if (DRAW_MAPS) {
  sched <- read_runs("acquisition_schedule_spatial.csv")
  if (is.null(sched) || !nrow(sched) || all(is.na(sched$x_coord))) {
    message("No mappable schedules found; skipping map_bud_foresight.")
  } else {
    ## WMD outlines (same approach as 10_results_assets.R) ----------------------
    eau_xy <- if (file.exists(EAU_LOOKUP))
      readRDS(EAU_LOOKUP)[, c("eau_id", "x_coord", "y_coord")] else NULL
    wmd_sf <- NULL; map_crs <- NULL
    if (WMD_OVERLAY && file.exists(WMD_RASTER) &&
        requireNamespace("terra", quietly = TRUE) && requireNamespace("sf", quietly = TRUE)) {
      tryCatch({
        r     <- terra::rast(WMD_RASTER)
        polys <- terra::as.polygons(r, dissolve = TRUE)
        if (!is.null(eau_xy)) {
          active <- terra::extract(r, as.matrix(eau_xy[, c("x_coord", "y_coord")]))[, 1]
          active <- sort(unique(active[!is.na(active)]))
          # as.data.frame(...)[[1]] gives a plain vector in every terra version;
          # polys[[name]] returns a one-column data.frame in recent terra, which
          # makes %in% FALSE for every polygon and silently drops all outlines.
          polys  <- polys[as.data.frame(polys)[[1]] %in% active, ]
        }
        wmd_sf  <- sf::st_as_sf(polys)
        map_crs <- sf::st_crs(wmd_sf)
      }, error = function(e) {
        message("WMD overlay skipped: ", conditionMessage(e))
        wmd_sf <<- NULL; map_crs <<- NULL
      })
    }

    ## classify each EAU, per budget x scenario ---------------------------------
    land <- sched %>% distinct(eau_id, x_coord, y_coord)
    acq  <- sched %>% filter(!is.na(acquired_year), model %in% c("rolling", "myopic"),
                             scenario %in% CLIM) %>%
      distinct(budget, scenario, model, eau_id) %>%
      mutate(flag = TRUE) %>%
      pivot_wider(names_from = model, values_from = flag)
    if (!"rolling" %in% names(acq)) acq$rolling <- NA
    if (!"myopic"  %in% names(acq)) acq$myopic  <- NA
    md <- tidyr::expand_grid(budget = runs$budget, scenario = CLIM, land) %>%
      left_join(acq, by = c("budget", "scenario", "eau_id")) %>%
      mutate(r = coalesce(rolling, FALSE), m = coalesce(myopic, FALSE),
             cat = factor(case_when(r & m ~ "Acquired by both",
                                    r     ~ "Rolling only (foresight adds)",
                                    m     ~ "Myopic only (foresight avoids)",
                                    TRUE  ~ "Neither"), levels = names(PAIR_COL)),
             scenario = factor(scenario, levels = CLIM),
             bud = factor(bud_lab(budget), levels = BUD_LEVELS)) %>%
      arrange(budget, scenario, cat != "Neither")     # coloured cells drawn on top

    mp <- ggplot(md, aes(x_coord, y_coord, fill = cat)) +
      geom_tile(width = CELL_W, height = CELL_H) +
      scale_fill_manual(values = PAIR_COL, name = NULL, drop = FALSE) +
      facet_grid(scenario ~ bud, labeller = labeller(scenario = SCEN_LAB1)) +
      guides(fill = guide_legend(nrow = 1, override.aes = list(colour = NA)))
    mp <- if (!is.null(wmd_sf)) {
      mp + geom_sf(data = wmd_sf, fill = NA, colour = WMD_LINE_COL,
                   linewidth = WMD_LINE_W, inherit.aes = FALSE) +
        coord_sf(crs = map_crs, datum = NA, expand = FALSE)
    } else mp + coord_equal()
    mp <- mp + theme_void(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
      theme(strip.text      = element_text(face = "bold", size = rel(0.95),
                                           margin = margin(2, 2, 4, 2)),
            strip.text.y    = element_text(face = "bold", angle = -90, size = rel(0.95)),
            legend.position = "bottom", panel.spacing = unit(4, "pt"),
            plot.background = element_rect(fill = "white", colour = NA))
    save_fig(mp, "map_bud_foresight", 10, 12.5)
    drew_maps <- TRUE
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# CHECK — the run at the canonical budget should reproduce the canonical run
# ═══════════════════════════════════════════════════════════════════════════════
canon_check <- NULL
if (BUDGET_CANONICAL %in% runs$budget && file.exists(CANON_RES)) {
  cid <- runs$run_id[runs$budget == BUDGET_CANONICAL]
  cc <- inner_join(
    read_csv(CANON_RES, show_col_types = FALSE) %>% select(scenario, model, va_c = value_added),
    res %>% filter(budget == BUDGET_CANONICAL) %>% select(scenario, model, va_b = value_added),
    by = c("scenario", "model")) %>%
    filter(scenario %in% CLIM) %>%
    mutate(diff_pct = 100 * (va_b - va_c) / va_c)
  canon_check <- list(run_id = cid, canonical_run_id = RUN_ID,
                      max_abs_diff_pct = max(abs(cc$diff_pct)),
                      identical = all(cc$va_b == cc$va_c))
  message("Canonical-budget check: ", cid, " vs ", RUN_ID, " | max |diff| = ",
          signif(canon_check$max_abs_diff_pct, 3), "% of value added",
          if (canon_check$identical) " (identical)" else "")
}

## metadata (the budget-comparison .qmd reads this) ------------------------------
jsonlite::write_json(
  list(generated_date    = format(Sys.Date(), "%Y-%m-%d"),
       base_run_id       = RUN_ID,
       lc_label          = LC_LABEL,
       canonical_budget  = BUDGET_CANONICAL,
       budgets           = runs$budget,
       run_ids           = runs$run_id,
       eau_area_ha       = EAU_AREA_HA,
       maps_drawn        = drew_maps,
       canon_check       = canon_check,
       r_version         = R.version.string),
  file.path(OUT_DIR, "bud_metadata.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)

message("\nDone. Wrote 3 tables + 3 figures", if (drew_maps) " + 1 map" else "",
        " + bud_metadata.json to '", OUT_DIR, "/'.")
})

run_section("Hazard sweep", {
# ═══════════════════════════════════════════════════════════════════════════════
# SECTION: hazard sweep (formerly 10d/10e) — assets across the hazard sweep
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
#   Rscript 10_results_assets.R
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
  skip_section("runs only in the canonical configuration (PPR_HAZARD and PPR_BUDGET unset).")
if (!file.exists(CANON_RES))
  skip_section("canonical results not found at ", CANON_RES, ". Run the canonical pipeline first.")
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
  skip_section("no complete hazard runs matching runs/", LC_SCENARIO, "_haz*. ",
               "Run 15_run_hazard_sweep.R first.")

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
})
