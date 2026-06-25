### DIAGNOSTIC: Gross vs. Net Pixel Transitions (Vote-Counting Check) ####
#
# PURPOSE
#   Decide whether to replace the current transition-probability proxy
#   (net proportional habitat loss, Script 05) with a "vote-counting" proxy
#   (fraction of suitable 30m pixels that flip suitable -> unsuitable).
#
#   The two proxies are arithmetically IDENTICAL unless there is simultaneous
#   gain and loss of suitable pixels within an EAU across a year pair. The
#   current method sees only NET change; vote-counting sees GROSS suitable->
#   unsuitable loss. This script quantifies the gap between them.
#
# WHAT IT COMPUTES, per EAU x consecutive-year-pair x RCP:
#   - gross_loss   : # pixels suitable at t that become unsuitable at t+1
#   - gross_gain   : # pixels unsuitable at t that become suitable at t+1
#   - suitable_t   : # suitable pixels at t (denominator)
#   - net_loss_rate    = (gross_loss - gross_gain) / suitable_t   [current proxy analogue]
#   - gross_loss_rate  =  gross_loss               / suitable_t   [vote-counting proxy]
#   - churn_ratio      =  gross_loss_rate / net_loss_rate (where net > 0)
#
# DECISION RULE (interpret after running):
#   If gross_loss_rate ~= net_loss_rate for most EAU-years (churn small),
#     -> vote-counting changes nothing; keep current proxy, document clearly.
#   If gross_loss_rate >> net_loss_rate (substantial churn),
#     -> vote-counting reveals hidden conversion risk; worth re-engineering
#        the transition data, and stochastic availability (Option B) becomes viable.
#
# INPUTS  : FOREsce rasters (RCP 4.5 and 8.5), wmd_raster_equal_area.tif,
#           eau_wmd_lookup.csv. Same paths/conventions as Script 02.
# OUTPUTS : input_data/diagnostic_pixel_transitions.csv  (full per-EAU-year table)
#           console summary tables + a diagnostic plot
#
# NOTE ON RESOLUTION: Unlike Script 02 (which resamples 30m -> EAU grid via
# "sum" and so only ever sees per-EAU aggregates), this script must compare
# pixels to pixels ACROSS years at native 30m resolution, then aggregate the
# transition counts up to the EAU. That is why it is a separate script and not
# a tweak to Script 02.
#
# RUNTIME: heavier than Script 02 because it crosses years at native resolution.
# Expect several minutes per RCP. Parallelized over year-pairs.
# ──────────────────────────────────────────────────────────────────────────────

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(parallel)
library(ggplot2)
library(scales)


# ══ PARAMETERS ════════════════════════════════════════════════════════════════
# Match Script 02 exactly so results are comparable.

scenarios <- list(
  `45` = "input_data/prairie_potholes_gcam_ref_rcp45",
  `85` = "input_data/prairie_potholes_gcam_ref_rcp85"
)

# Suitable land cover classes (identical to Script 02):
# 19 = Perennial Grass, 20 = Open Water, 26 = Grassland,
# 28 = Woody Wetland, 29 = Herbaceous Wetland
suitable_classes <- c(19, 20, 26, 28, 29)

# EAU zone raster (WMD names as cell values) and the FOREsce-shared mask.
eau_zone_path   <- "input_data/wmd_raster_equal_area.tif"
shared_mask_path <- "input_data/foresce_eau_shared_mask.tif"
lookup_path     <- "input_data/eau_wmd_lookup.csv"

# Restrict the diagnostic to the decision-relevant horizon. The transition that
# feeds trans_prob in the ILP is the forward difference between decadal decision
# years (2030->2040, ..., 2090->2100), plus the 2020->2030 boundary step.
# FOREsce is typically published in 10-year steps; if yours is annual, this still
# works (it just compares whatever consecutive files exist). Set to NULL to use
# ALL consecutive year-pairs found on disk.
restrict_to_years <- c(2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100)
# (2014 is intentionally excluded: it predates the decision horizon, plays no
#  role in Scripts 03-05, and its 2014->2020 pair spans 6 years, not a decade,
#  which would put its loss rate on a different footing from the decadal pairs.)

# Parallel cores.
n_cores <- max(1, parallel::detectCores() - 1)

# Output path for the full per-EAU-year transition table.
out_csv <- "input_data/diagnostic_pixel_transitions.csv"
# ══════════════════════════════════════════════════════════════════════════════


# ── 1. Load the EAU zone raster at native FOREsce resolution ──────────────────
# The challenge: wmd_raster_equal_area.tif is at the COARSE EAU resolution
# (one cell per EAU). To attribute 30m pixel transitions to EAUs, we need an
# eau_id label on the FINE grid. We build that by:
#   (a) loading the coarse EAU zone raster (WMD names),
#   (b) reconstructing sequential eau_id on the coarse grid in the EXACT order
#       Script 01 used (arrange by wmd_id_num, then original cell number),
#   (c) later: resampling that fine-labeled grid to each FOREsce raster via
#       nearest-neighbour so zone IDs are never interpolated.

cat("Loading spatial inputs...\n")

eau_zone <- rast(eau_zone_path)      # coarse: WMD names per cell
shared_mask <- rast(shared_mask_path) # coarse: EAU x FOREsce intersection
lookup <- read_csv(lookup_path, show_col_types = FALSE)

cat("  EAU zone CRS: ", crs(eau_zone, proj = TRUE), "\n")
cat("  EAU zone res: ", res(eau_zone)[1], "m\n")
cat("  EAUs in lookup: ", nrow(lookup), "\n")

# Reconstruct an eau_id raster on the COARSE grid, matching Script 06's approach:
# place each lookup eau_id at the cell containing its recorded centroid.
eau_cells <- cellFromXY(eau_zone, cbind(lookup$x_coord, lookup$y_coord))
if (any(is.na(eau_cells)))
  stop("cellFromXY returned NA for ", sum(is.na(eau_cells)),
       " centroid(s). CRS mismatch between lookup coords and eau_zone?")
if (any(duplicated(eau_cells)))
  stop("cellFromXY returned duplicate cells for ", sum(duplicated(eau_cells)),
       " EAU(s).")

eau_id_coarse <- rast(ext(eau_zone), res = res(eau_zone), crs = crs(eau_zone))
values(eau_id_coarse) <- NA_integer_
eau_id_coarse[eau_cells] <- lookup$eau_id

# Apply the FOREsce-shared mask so we only ever attribute pixels for EAUs that
# actually have land-cover data (drops the 42 out-of-extent EAUs, consistent
# with Script 02).
eau_id_coarse <- mask(eau_id_coarse, shared_mask)

cat("  Non-NA EAU cells after shared mask: ",
    sum(!is.na(values(eau_id_coarse))), "\n")


# ── 2. Helper: list FOREsce files for a scenario, sorted by year ──────────────
list_foresce_files <- function(LU_in, suffix) {
  fnames <- list.files(
    LU_in,
    pattern = paste0("prairie_potholes_gcam_ref_rcp", suffix, "_[0-9]{4}\\.tif$"),
    full.names = FALSE
  )
  if (length(fnames) == 0)
    stop("Scenario ", suffix, " has no matching rasters in ", LU_in)
  
  years <- as.integer(substr(fnames, nchar(fnames) - 7, nchar(fnames) - 4))
  o <- order(years)
  tibble(fname = fnames[o], year = years[o], path = file.path(LU_in, fnames[o]))
}


# ── 3. Helper: classify a FOREsce raster to a binary suitable mask ────────────
# 1 = suitable, 0 = unsuitable, NA = outside coverage. Cropped to the EAU extent
# and snapped to the fine FOREsce grid (NOT resampled to the coarse EAU grid —
# that is the whole point: we keep 30m pixels).
to_suitable_mask <- function(path, template_extent) {
  r <- rast(path)
  r <- crop(r, template_extent)
  # suitable -> 1, everything else -> 0, preserving NA outside data
  m <- classify(r,
                rcl = cbind(suitable_classes, 1),
                others = 0)
  m
}


# ── 4. Core: transition counts for one consecutive year-pair ──────────────────
# Returns a data.frame: eau_id, year (= start year), rcp, gross_loss,
# gross_gain, suitable_t, suitable_tp1, total_pixels.
#
# Logic per pixel:
#   suitable_t = 1 & suitable_tp1 = 0  -> LOSS  (suitable -> unsuitable)
#   suitable_t = 0 & suitable_tp1 = 1  -> GAIN  (unsuitable -> suitable)
# We zonal-sum each indicator over the fine eau_id grid.
transition_for_pair <- function(path_t, path_tp1, year_t, suffix,
                                eau_id_coarse, eau_extent) {
  
  s_t   <- to_suitable_mask(path_t,   eau_extent)
  s_tp1 <- to_suitable_mask(path_tp1, eau_extent)
  
  # Align the two years to a common fine grid (they should already match; this
  # guards against off-by-one extent/origin differences between files).
  s_tp1 <- resample(s_tp1, s_t, method = "near")
  
  # Build the fine-resolution eau_id label by resampling the coarse eau_id grid
  # onto the FOREsce grid with NEAREST NEIGHBOUR (never interpolate zone IDs).
  eau_id_fine <- resample(eau_id_coarse, s_t, method = "near")
  
  # Per-pixel transition indicators (NA-safe: NA in either year -> NA here).
  loss_ind <- (s_t == 1) & (s_tp1 == 0)   # suitable -> unsuitable
  gain_ind <- (s_t == 0) & (s_tp1 == 1)   # unsuitable -> suitable
  
  # Zonal sums over EAUs. zonal() with a logical/0-1 raster + fun="sum" counts.
  z_loss <- zonal(loss_ind,        eau_id_fine, fun = "sum", na.rm = TRUE)
  z_gain <- zonal(gain_ind,        eau_id_fine, fun = "sum", na.rm = TRUE)
  z_st   <- zonal((s_t   == 1),    eau_id_fine, fun = "sum", na.rm = TRUE)
  z_stp1 <- zonal((s_tp1 == 1),    eau_id_fine, fun = "sum", na.rm = TRUE)
  z_tot  <- zonal(!is.na(s_t),     eau_id_fine, fun = "sum", na.rm = TRUE)
  
  names(z_loss) <- c("eau_id", "gross_loss")
  names(z_gain) <- c("eau_id", "gross_gain")
  names(z_st)   <- c("eau_id", "suitable_t")
  names(z_stp1) <- c("eau_id", "suitable_tp1")
  names(z_tot)  <- c("eau_id", "total_pixels")
  
  out <- reduce(
    list(z_loss, z_gain, z_st, z_stp1, z_tot),
    left_join, by = "eau_id"
  )
  out$year <- year_t
  out$rcp  <- suffix
  out
}


# ── 5. Run all year-pairs for both scenarios ──────────────────────────────────
eau_extent <- ext(eau_zone)

run_scenario <- function(suffix) {
  files <- list_foresce_files(scenarios[[suffix]], suffix)
  
  if (!is.null(restrict_to_years)) {
    files <- files %>% filter(year %in% restrict_to_years)
  }
  
  if (nrow(files) < 2)
    stop("Scenario ", suffix, " has fewer than 2 usable rasters; ",
         "cannot form a year-pair.")
  
  # Consecutive pairs: (row 1, row 2), (row 2, row 3), ...
  pair_idx <- seq_len(nrow(files) - 1)
  
  cat(sprintf("\nScenario RCP %s: %d rasters -> %d consecutive pairs\n",
              suffix, nrow(files), length(pair_idx)))
  print(files %>% select(year, fname))
  
  t0 <- Sys.time()
  
  res_list <- mclapply(pair_idx, function(k) {
    transition_for_pair(
      path_t    = files$path[k],
      path_tp1  = files$path[k + 1],
      year_t    = files$year[k],
      suffix    = suffix,
      eau_id_coarse = eau_id_coarse,
      eau_extent    = eau_extent
    )
  }, mc.cores = n_cores)
  
  cat(sprintf("  RCP %s done in %.1f min\n", suffix,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  
  bind_rows(res_list)
}

all_transitions <- map_dfr(names(scenarios), run_scenario)


# ── 6. Derive the two competing rates and the churn ratio ─────────────────────
# Keep only real EAUs (in the lookup). Guard against divide-by-zero.
trans_tbl <- all_transitions %>%
  filter(eau_id %in% lookup$eau_id) %>%
  mutate(
    net_loss        = gross_loss - gross_gain,
    net_loss_rate   = if_else(suitable_t > 0, net_loss   / suitable_t, NA_real_),
    gross_loss_rate = if_else(suitable_t > 0, gross_loss / suitable_t, NA_real_),
    # churn ratio only meaningful where there is genuine net loss (> 0)
    churn_ratio     = if_else(net_loss_rate > 0,
                              gross_loss_rate / net_loss_rate, NA_real_),
    # absolute gap in rate space (gross minus net); always >= 0 by construction
    rate_gap        = gross_loss_rate - net_loss_rate
  )

write_csv(trans_tbl, out_csv)
cat(sprintf("\nFull per-EAU-year transition table written to: %s\n", out_csv))


# ── 7. SUMMARY TABLES ─────────────────────────────────────────────────────────

cat("\n=======================================================\n")
cat("  VOTE-COUNTING DIAGNOSTIC: gross vs net suitable loss\n")
cat("=======================================================\n\n")

# How many EAU-years show simultaneous gain AND loss (i.e., any churn at all)?
churn_presence <- trans_tbl %>%
  summarise(
    n_eau_years          = n(),
    n_with_any_loss      = sum(gross_loss > 0, na.rm = TRUE),
    n_with_any_gain      = sum(gross_gain > 0, na.rm = TRUE),
    n_with_both          = sum(gross_loss > 0 & gross_gain > 0, na.rm = TRUE),
    pct_with_both        = round(100 * n_with_both / n_eau_years, 1)
  )
cat("-- Presence of simultaneous gain & loss --\n")
print(as.data.frame(churn_presence), row.names = FALSE)

# Rate comparison by RCP x year. The headline numbers: if mean gross_loss_rate
# is materially above mean net_loss_rate, vote-counting matters.
cat("\n-- Mean rates by RCP x year (the headline comparison) --\n")
rate_by_year <- trans_tbl %>%
  group_by(rcp, year) %>%
  summarise(
    n               = n(),
    mean_net_rate   = round(mean(net_loss_rate,   na.rm = TRUE), 6),
    mean_gross_rate = round(mean(gross_loss_rate, na.rm = TRUE), 6),
    median_net      = round(median(net_loss_rate,   na.rm = TRUE), 6),
    median_gross    = round(median(gross_loss_rate, na.rm = TRUE), 6),
    .groups = "drop"
  ) %>%
  mutate(gross_over_net = round(mean_gross_rate / mean_net_rate, 2))
print(as.data.frame(rate_by_year), row.names = FALSE)

# Distribution of the churn ratio among EAU-years with genuine net loss.
cat("\n-- Churn ratio (gross_loss_rate / net_loss_rate), net-losing EAU-years only --\n")
cat("   (1.0 = no hidden churn; >>1 = lots of offsetting gain hidden under net) \n\n")
churn_dist <- trans_tbl %>%
  filter(!is.na(churn_ratio), is.finite(churn_ratio)) %>%
  summarise(
    n        = n(),
    min      = round(min(churn_ratio), 2),
    p25      = round(quantile(churn_ratio, 0.25), 2),
    median   = round(median(churn_ratio), 2),
    mean     = round(mean(churn_ratio), 2),
    p75      = round(quantile(churn_ratio, 0.75), 2),
    p90      = round(quantile(churn_ratio, 0.90), 2),
    max      = round(max(churn_ratio), 2)
  )
print(as.data.frame(churn_dist), row.names = FALSE)

# How many EAUs would move OFF the epsilon floor under vote-counting?
# Compare against the Script 05 epsilon (data-derived, ~1e-7 to 1e-8). We can't
# read Script 05's epsilon here, so we report the share of EAU-years whose
# gross_loss_rate exceeds several candidate thresholds.
cat("\n-- Share of EAU-years with gross_loss_rate above candidate floors --\n")
thresholds <- c(1e-6, 1e-4, 1e-3, 1e-2, 5e-2)
floor_tbl <- tibble(threshold = thresholds) %>%
  rowwise() %>%
  mutate(
    pct_net_above   = round(100 * mean(trans_tbl$net_loss_rate   > threshold, na.rm = TRUE), 1),
    pct_gross_above = round(100 * mean(trans_tbl$gross_loss_rate > threshold, na.rm = TRUE), 1)
  ) %>%
  ungroup()
print(as.data.frame(floor_tbl), row.names = FALSE)


# ── 8. DIAGNOSTIC PLOT ────────────────────────────────────────────────────────
# Paired histogram of net vs gross loss rates (log-ish x via signed sqrt) so you
# can SEE whether the gross distribution shifts right of the net distribution.

plot_df <- trans_tbl %>%
  select(rcp, year, net_loss_rate, gross_loss_rate) %>%
  pivot_longer(c(net_loss_rate, gross_loss_rate),
               names_to = "proxy", values_to = "rate") %>%
  mutate(
    proxy = recode(proxy,
                   net_loss_rate   = "Net (current proxy)",
                   gross_loss_rate = "Gross (vote-counting)"),
    rcp_label = paste0("RCP ", ifelse(rcp == "45", "4.5", "8.5"))
  )

theme_thesis <- theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey92", color = NA),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom")

p_diag <- plot_df %>%
  filter(!is.na(rate)) %>%
  ggplot(aes(x = pmax(0, rate), fill = proxy)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 60,
                 color = NA) +
  scale_fill_manual(values = c("Net (current proxy)"   = "#2171b5",
                               "Gross (vote-counting)" = "#cb181d"),
                    name = NULL) +
  scale_x_sqrt(labels = percent_format(accuracy = 0.1)) +
  facet_grid(rcp_label ~ ., scales = "free_y") +
  labs(
    title    = "Net vs. Gross Suitable-Habitat Loss Rates per EAU-Year",
    subtitle = "If red (gross) sits noticeably right of blue (net), hidden churn is real and vote-counting matters.",
    x        = "Per-period loss rate (sqrt scale)",
    y        = "Count of EAU-years"
  ) +
  theme_thesis

print(p_diag)

ggsave("input_data/diagnostic_vote_counting_plot.png", p_diag,
       width = 9, height = 6, dpi = 150)

cat("\n=======================================================\n")
cat("  INTERPRETATION GUIDE\n")
cat("=======================================================\n")
cat("  * Look at 'gross_over_net' in the rate-by-year table and the\n")
cat("    median/mean churn ratio. Values near 1.0 => keep current proxy.\n")
cat("  * Values of, say, 1.5x-3x+ => vote-counting reveals real conversion\n")
cat("    risk; worth re-engineering trans_prob and reconsidering Option B.\n")
cat("  * Also check 'pct_gross_above' vs 'pct_net_above': if vote-counting\n")
cat("    pushes many more EAU-years above a ~1e-3/1e-2 floor, the risk\n")
cat("    landscape changes materially.\n")
cat("  * Plot saved to input_data/diagnostic_vote_counting_plot.png\n")
cat("=======================================================\n\n")