# ═══════════════════════════════════════════════════════════════════════════════
# 10c_results_assets_budget.R — OPTIONAL: assets comparing budget-sensitivity runs
# ═══════════════════════════════════════════════════════════════════════════════
# Compares the budget runs made by 13_run_budget_sweep.R: every complete folder
# runs/<LC_SCENARIO>_bXX/ (e.g. gcam_ref_b02, gcam_ref_b05, gcam_ref_b10). The budget
# of each run is read from its own model_results.rds (params$budget_eaus_per_period),
# not from the folder name. Formulas are the same as in 10_results_assets.R.
#
# PREREQUISITES:
#   - the canonical run, runs/<LC_SCENARIO>/ (Rscript 00_run_all.R, PPR_BUDGET unset)
#   - at least two complete budget runs (Rscript 13_run_budget_sweep.R)
# Run with PPR_BUDGET UNSET, from the project root:
#   Rscript 10c_results_assets_budget.R
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

# Checked before setup, which would otherwise create empty folders for a budget run.
if (nzchar(Sys.getenv("PPR_BUDGET")))
  stop("Run 10c with PPR_BUDGET unset: it writes into the canonical run's ",
       "output_figs. In this shell: unset PPR_BUDGET", call. = FALSE)

if (!exists(".SETUP_DONE")) source("00_setup.R")
source("00_config.R")   # always re-read: a long-lived R session can hold an older config

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
  stop("Found ", nrow(runs), " complete budget run(s) matching runs/", LC_SCENARIO,
       "_bXX. Run 13_run_budget_sweep.R first.", call. = FALSE)
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
