# ═══════════════════════════════════════════════════════════════════════════════
# 10_results_assets.R  —  Asset generator (figs, tables, maps) for the PPR acquisition results document
# ═══════════════════════════════════════════════════════════════════════════════
# Writes to OUT_DIR:
#
# Tables are written as tidy CSVs of RAW NUMBERS; all formatting and the grouped-header
# styling happen in results.qmd. Same approach with Figs/Maps

# ═══════════════════════════════════════════════════════════════════════════════

if (!exists(".SETUP_DONE")) source("00_setup.R")

# ── CONFIG ─────────────────────────────────────────────────────────────────────
IN_DIR   <- DIR_OUT      # "output_data"
OUT_DIR  <- DIR_FIGS     # "output_figs"
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
CELL_M       <- sqrt(282e6)               # EAU side length in metres (tile width/height)
MAP_BG       <- "grey91"                  # unacquired-EAU landscape colour
WMD_OVERLAY  <- TRUE                       # draw WMD outlines if raster + terra/sf available
WMD_RASTER   <- "derived_data/wmd_raster_equal_area.tif"
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
PANEL_RDS     <- "derived_data/data_panel.rds"      # inputs panel from script 06
EAU_LOOKUP    <- "derived_data/eau_wmd_lookup.rds"  # eau_id -> x_coord/y_coord
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

# Table 1 — value created over the do-nothing baseline (all three models)
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

# Table 2 — value-added gap vs. rolling foresight (greedy, myopic)
tbl_gap_vs_rolling <- gaps %>%
  transmute(scenario, model,
            gap_pct = gap_value_added_pct, dpairs = dpairs) %>%
  pivot_wider(names_from = model, values_from = c(gap_pct, dpairs)) %>%
  arrange(scenario) %>%
  transmute(scenario,
            greedy_gap_pct = gap_pct_greedy, greedy_dpairs = dpairs_greedy,
            myopic_gap_pct = gap_pct_myopic, myopic_dpairs = dpairs_myopic)
write_csv(tbl_gap_vs_rolling, file.path(OUT_DIR, "table_gap_vs_rolling.csv"))

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE (→ Figure 1)  —  cumulative gap vs. rolling over the horizon (climate scenarios)
# ═══════════════════════════════════════════════════════════════════════════════
roll_disc <- trj %>% filter(model == "rolling") %>% select(scenario, year, roll = discounted)
cumg <- trj %>%
  filter(model %in% c("greedy", "myopic"), scenario %in% CLIM) %>%
  left_join(roll_disc, by = c("scenario", "year")) %>%
  arrange(scenario, model, year) %>%
  group_by(scenario, model) %>%
  mutate(cum_gap = cumsum(roll - discounted)) %>%
  ungroup() %>%
  mutate(model = factor(model, levels = c("greedy", "myopic")),
         scenario = factor(scenario, levels = CLIM))
ends <- cumg %>% group_by(scenario, model) %>% filter(year == max(year)) %>% ungroup()

fig_overtime <- ggplot(cumg, aes(year, cum_gap, colour = model)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  geom_text(data = ends, aes(label = fmt_pairs(cum_gap)), hjust = -0.12, size = 2.8,
            show.legend = FALSE) +
  facet_wrap(~scenario, ncol = 2, labeller = as_labeller(SCEN_LAB1)) +
  scale_colour_manual(values = COL_MODEL, labels = c("Greedy", "Myopic")) +
  scale_x_continuous("year", expand = expansion(mult = c(0.02, 0.14))) +
  scale_y_continuous("Cumulative gap vs. rolling\n(discounted duck-pairs)", labels = label_comma()) +
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
  # EAU area derives from CELL_M so there is a single source of truth for cell size
  EAU_AREA_HA <- CELL_M^2 / 1e4
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
        geom_tile(width = CELL_M, height = CELL_M))
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
# acquired_period/year, prevented_loss). Cells are 282 km^2 EAUs in USGS Conus Albers
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
    geom_tile(width = CELL_M, height = CELL_M) +
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
    geom_tile(data = land, aes(x_coord, y_coord), width = CELL_M, height = CELL_M, fill = MAP_BG) +
    geom_tile(data = dplyr::filter(md_dv, dV != 0),
              aes(x_coord, y_coord, fill = dV), width = CELL_M, height = CELL_M) +
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
jsonlite::write_json(
  list(generated      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
       generated_date = format(Sys.Date(), "%Y-%m-%d"),
       maps_drawn     = drew_maps,
       r_version      = R.version.string),
  file.path(OUT_DIR, "run_metadata.json"), auto_unbox = TRUE, pretty = TRUE)

message("\nDone. Wrote 2 tables + 2 figures",
        if (drew_maps) " + 2 maps" else "",
        if (drew_desc) " + 4 descriptive assets" else "",
        " + run_metadata.json to '", OUT_DIR, "/'.")
