# ═══════════════════════════════════════════════════════════════════════════════
# 10_results_assets.R  —  Asset generator (figs, tables, maps) for the PPR acquisition results document
# ═══════════════════════════════════════════════════════════════════════════════
# Writes to OUT_DIR:
#
# Tables are written as tidy CSVs of RAW NUMBERS; all formatting and the grouped-header
# styling happen in results.qmd. Same approach with Figs/Maps

# ═══════════════════════════════════════════════════════════════════════════════

if (!isTRUE(.SETUP_DONE)) source("00_setup.R")

# ── CONFIG ─────────────────────────────────────────────────────────────────────
RATE_LABEL <- "data_derived"   # or: wetland_high, grass_low, grass_high
IN_DIR     <- paste0("output_data_", RATE_LABEL)
OUT_DIR    <- paste0("output_figs/_", RATE_LABEL)
RESULTS  <- file.path(IN_DIR, "model_results.csv")
TRAJ     <- file.path(IN_DIR, "model_trajectories.csv")
SCHEDULE <- file.path(IN_DIR, "acquisition_schedule_spatial.csv")  # 08's persisted schedule (maps)

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
WMD_RASTER   <- "input_data/wmd_raster_equal_area.tif"
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
# MAPS (→ Figures 3–4)  —  acquisition footprints from 08's persisted schedule
# Reads acquisition_schedule_spatial.csv (scenario, model, eau_id, wmd_id, x/y,
# acquired_period/year, prevented_loss). Cells are 282 km^2 EAUs in USGS Conus Albers
# (ESRI:102039); the optional WMD overlay dissolves wmd_raster_equal_area.tif to outlines.
# ═══════════════════════════════════════════════════════════════════════════════
if (DRAW_MAPS && file.exists(SCHEDULE)) {

  sched <- read_csv(SCHEDULE, show_col_types = FALSE)
  land  <- sched %>% distinct(eau_id, x_coord, y_coord)   # landscape backdrop
  acq   <- function(df) df %>% filter(!is.na(acquired_year))

  ## optional WMD outline layer (dissolve labelled EAU raster -> active WMD polygons)
  wmd_sf <- NULL; map_crs <- NULL
  if (WMD_OVERLAY && file.exists(WMD_RASTER) &&
      requireNamespace("terra", quietly = TRUE) && requireNamespace("sf", quietly = TRUE)) {
    ok <- tryCatch({
      r      <- terra::rast(WMD_RASTER)
      active <- terra::extract(r, as.matrix(land[, c("x_coord", "y_coord")]))[, 1]
      active <- sort(unique(active[!is.na(active)]))
      polys  <- terra::as.polygons(r, dissolve = TRUE)
      polys  <- polys[polys[[names(polys)[1]]] %in% active, ]
      wmd_sf  <<- sf::st_as_sf(polys); map_crs <<- sf::st_crs(wmd_sf); TRUE
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
       rate_label     = RATE_LABEL,
       maps_drawn     = drew_maps,
       r_version      = R.version.string),
  file.path(OUT_DIR, "run_metadata.json"), auto_unbox = TRUE, pretty = TRUE)

message("\nDone. Wrote 2 tables + 2 figures",
        if (drew_maps) " + 2 maps" else "",
        " + run_metadata.json to '", OUT_DIR, "/'.")
