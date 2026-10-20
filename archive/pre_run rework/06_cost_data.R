### Cost Data (Script 06) ####
#
# Loads the PLACES FMV: vacant raster (Nolte 2020)
#
# Operationalization:
#   cost[i, t] = mean_fmv_per_ha[i] × area_ha × (1 + r)^(t − base_year)
#
# where:
#   mean_fmv_per_ha[i]  = arithmetic mean of PLACES FMV ($/ha) across all
#                         fine-resolution pixels within EAU i
#   area_ha             = EAU area in hectares (equal for all EAUs by construction)
#   r                   = 0.02 (annual inflation rate)
#   base_year           = 2017 (dollar denomination year for PLACES FMV)
#
# Cost is scenario-invariant: it does not vary by rcp or gcm.
# Cost varies by eau_id and year only.
#
# Data source:
#   Nolte, C. (2020). High-resolution land value maps reveal underestimation
#   of conservation costs in the United States. PNAS, 117(47), 29577–29583.
#   Data: https://datadryad.org/dataset/doi:10.5061/dryad.np5hqbzq9
#
#   File: places_fmv_vacant.tif
#   Values: natural log of USD per hectare (2017 USD), from 2010 market conditions.
#   use the vacant-land variant (trained on sales without structures), which is
#   the appropriate specification for conservation land acquisition cost estimation.
#   Pre-download to input_data/ before running this script.

#
# Edge EAU design note:
#   cost is computed here only over PLACES pixels that fall
#   within WMD boundaries (mask_to_wmd = TRUE). Edge EAUs are thereby costed only
#   on the portion of their footprint for which benefit data also exist. The
#   pixel-count diagnostic in Step 7 records how many pixels contribute to each
#   EAU's mean, making it straightforward to identify which EAUs are edge cases
#   and how much of their area falls outside WMD boundaries.
#
# INPUTS:  data_panel (from script 05), eau_wmd_lookup.csv,
#          wmd_raster_equal_area.tif, places_fmv_vacant.tif,
#          derived_data/Wetland_Management_Districts/FSMS_WMD.shp
# OUTPUTS: data_panel with cost column populated; saved to input_data/

#------ Setup ________
if (!exists(".SETUP_DONE")) source("00_setup.R")

# ── 0. Load data ──────────────────────────────────────####
data_panel <- readRDS("derived_data/panel_05_risk.rds")
eau_wmd    <- read_csv("derived_data/eau_wmd_lookup.csv", show_col_types = FALSE)


# ══ PARAMETERS ════════════════════════════════════════════════════════════════####
inflation_rate <- 0.02    # annual inflation rate for cost projection
base_year      <- 2017    # dollar denomination year for PLACES FMV: vacant
places_path    <- "input_data/places_fmv_pnas_dryad/1 estimates/places_fmv_vacant.tif"
wmd_shp_path   <- "input_data/Wetland_Management_Districts/FSMS_WMD.shp"

# Edge EAU masking (see design note in header):
#   TRUE  — cost averaged only over PLACES pixels inside WMD boundaries. This is
#           the primary approach: it ensures cost and benefit are computed over the
#           same spatial domain, eliminating the cost-benefit asymmetry for edge
#           EAUs whose rectangular footprints extend outside WMD polygons.
#           Requires wmd_shp_path above.
#   FALSE — cost averaged over full rectangular EAU footprint. Retained as a
#           toggle for sensitivity comparison; not the primary analysis.
mask_to_wmd <- TRUE
# ══════════════════════════════════════════════════════════════════════════════


# ── 1. Check that PLACES raster is present ────────────────────────────────────####

if (!file.exists(places_path)) {
  stop(
    "PLACES FMV raster not found at: ", places_path, "\n\n",
    "Download places_fmv_vacant.tif from:\n",
    "  https://datadryad.org/dataset/doi:10.5061/dryad.np5hqbzq9\n",
    "and save it to: input_data/places_fmv_vacant.tif\n\n",
    "Values in the file are ln($2017/ha); this script handles the ",
    "back-transformation."
  )
}


# ── 2. Load spatial inputs ─────────────────────────────────────────────────────####

cat("Loading spatial inputs...\n")

eau_r   <- rast(file.path(DIR_DERIVED, "wmd_raster_equal_area.tif"))  # EAU zone raster (WMD names)
fmv_raw <- rast(places_path)                             # PLACES FMV: ln($/ha), 2017 USD

cat("  EAU raster CRS:     ", crs(eau_r,   proj = TRUE), "\n")
cat("  EAU raster res:     ", res(eau_r)[1], "m ×", res(eau_r)[2], "m\n")
cat("  PLACES raster CRS:  ", crs(fmv_raw, proj = TRUE), "\n")
cat("  PLACES raster res:  ", res(fmv_raw)[1], "m\n")
cat(sprintf("  PLACES value range: %.2f to %.2f [ln($/ha)]\n",
            global(fmv_raw, "min", na.rm = TRUE)[[1]],
            global(fmv_raw, "max", na.rm = TRUE)[[1]]))
cat("  EAUs in lookup:     ", nrow(eau_wmd), "\n")


# ── 3. Crop PLACES raster to PPR extent ───────────────────────────────────────####
# Project the EAU extent into the PLACES CRS and crop before reprojection.
# This avoids reprojecting the full continental US raster, which is slow.

cat("\nCropping PLACES raster to PPR extent...\n")

ppr_ext_in_fmv_crs <- project(ext(eau_r), from = crs(eau_r), to = crs(fmv_raw))
fmv_crop <- crop(fmv_raw, ppr_ext_in_fmv_crs)

cat("  Cropped extent: ", as.character(ext(fmv_crop)), "\n")
cat("  Cells retained: ", ncell(fmv_crop), "\n")


# ── 4. Back-transform: ln($/ha) → $/ha ────────────────────────────────────────####
# PLACES stores values as natural log of USD per hectare.
# We exponentiate BEFORE aggregating so that the zonal mean is the arithmetic
# mean of $/ha values across pixels, not the exponential of the mean log (geometric
# mean). The arithmetic mean is the correct basis for estimating per-hectare
# acquisition cost across an EAU.

cat("\nBack-transforming: exp(ln($/ha)) → $/ha...\n")

fmv_per_ha_crop <- exp(fmv_crop)

cat(sprintf("  $/ha range: $%s  to  $%s\n",
            format(round(global(fmv_per_ha_crop, "min", na.rm = TRUE)[[1]]),
                   big.mark = ","),
            format(round(global(fmv_per_ha_crop, "max", na.rm = TRUE)[[1]]),
                   big.mark = ",")))


# ── 5. Reproject PLACES $/ha raster to EAU CRS ────────────────────────────────####
# Reproject at the native ~480m PLACES resolution. We preserve fine-grained
# values here and aggregate up to EAU scale via zonal statistics in step 7,
# rather than snapping to the coarse EAU grid (which would discard within-EAU
# heterogeneity used to compute a better mean).

cat("\nReprojecting PLACES raster to EAU CRS (bilinear interpolation)...\n")

fmv_reproj <- project(fmv_per_ha_crop, crs(eau_r), method = "bilinear")

cat("  Reprojected resolution: ", res(fmv_reproj)[1], "m\n")
cat("  Reprojected extent:     ", as.character(ext(fmv_reproj)), "\n")


# ── 6. Build eau_id raster ────────────────────────────────────────────────────####
# The EAU raster (wmd_raster_equal_area.tif) stores WMD names as cell values,
# not sequential eau_ids. We reconstruct an eau_id raster by:
#   (a) Using cellFromXY() to find the cell number in eau_r for each EAU centroid
#       recorded in the lookup table (x_coord / y_coord are those centroids).
#   (b) Creating a new raster with the same geometry as eau_r, then assigning
#       each cell the corresponding sequential eau_id.
#
# This exactly reproduces the Script 01 ordering (arrange by wmd_id_num, then
# original cell number → row_number() = eau_id).

cat("\nBuilding eau_id raster...\n")

eau_cells <- cellFromXY(eau_r, cbind(eau_wmd$x_coord, eau_wmd$y_coord))

# Internal consistency checks before proceeding
if (any(is.na(eau_cells))) {
  stop("cellFromXY() returned NA for ", sum(is.na(eau_cells)),
       " EAU centroid(s). CRS mismatch between eau_wmd coords and eau_r?")
}
if (any(duplicated(eau_cells))) {
  stop("cellFromXY() returned duplicate cell numbers for ",
       sum(duplicated(eau_cells)), " EAU(s). ",
       "Two EAUs may map to the same raster cell.")
}
if (length(eau_cells) != nrow(eau_wmd)) {
  stop("cellFromXY() returned ", length(eau_cells), " values but eau_wmd has ",
       nrow(eau_wmd), " rows.")
}

# Create blank raster template (same extent / res / CRS as eau_r)
eau_id_r <- rast(ext(eau_r), res = res(eau_r), crs = crs(eau_r))
values(eau_id_r) <- NA_integer_
eau_id_r[eau_cells] <- eau_wmd$eau_id

cat("  Non-NA cells in eau_id raster: ", sum(!is.na(values(eau_id_r))),
    "(expected:", nrow(eau_wmd), ")\n")


# ── 7. Zonal statistics: mean $/ha per EAU ────────────────────────────────────####
#
# Primary approach (mask_to_wmd = TRUE): resample the eau_id raster to PLACES
# resolution (nearest neighbor, so zone IDs are never interpolated), mask the
# PLACES raster to WMD boundaries so pixels outside any WMD are set to NA,
# then compute the zonal mean over only the within-WMD pixels for each EAU.
# This ensures cost and benefit are computed over the same spatial domain.
# For interior EAUs (fully within a WMD) the mask has no effect. For edge EAUs
# (footprints that clip a WMD boundary) the mask excludes outside pixels, so cost
# reflects only the portion of the EAU for which benefit data also exist.
# The sensitivity alternative (mask_to_wmd = FALSE) uses the full rectangular
# footprint and is retained for comparison.
#
# n_pixels per EAU is recorded in both cases to quantify edge EAU coverage.

cat("\nResampling eau_id raster to PLACES resolution (nearest neighbor)...\n")
eau_id_resampled <- resample(eau_id_r, fmv_reproj, method = "near")

# Mask PLACES raster to WMD boundaries (primary approach) or use full footprint
fmv_for_zonal <- fmv_reproj    # fallback; overwritten below if mask_to_wmd = TRUE

if (mask_to_wmd) {
  if (!file.exists(wmd_shp_path)) {
    stop("mask_to_wmd = TRUE but WMD shapefile not found at: ", wmd_shp_path)
  }
  cat("Applying WMD boundary mask to PLACES raster (mask_to_wmd = TRUE)...\n")
  wmd_sf        <- st_read(wmd_shp_path, quiet = TRUE)
  wmd_vect      <- project(vect(wmd_sf), crs(fmv_reproj))
  fmv_for_zonal <- mask(fmv_reproj, wmd_vect)
  n_inside  <- sum(!is.na(values(fmv_for_zonal)))
  n_outside <- sum(!is.na(values(fmv_reproj))) - n_inside
  cat(sprintf("  Pixels inside WMD boundaries:   %d\n", n_inside))
  cat(sprintf("  Pixels masked to NA (outside):  %d\n", n_outside))
} else {
  cat("Using full rectangular EAU footprints (mask_to_wmd = FALSE).\n")
}

# Compute zonal mean AND pixel count (count = non-NA PLACES pixels per zone).
# Pixel count is the key diagnostic for edge EAUs: an EAU with substantially
# fewer pixels than typical has a large fraction of its footprint outside either
# the WMD boundary (if mask_to_wmd = TRUE) or the PLACES coverage extent.
cat("Computing zonal mean FMV ($/ha) and pixel count per EAU...\n")

fmv_mean  <- zonal(fmv_for_zonal, eau_id_resampled, fun = "mean",  na.rm = TRUE)
fmv_count <- zonal(!is.na(fmv_for_zonal), eau_id_resampled, fun = "sum",   na.rm = TRUE)
names(fmv_mean)  <- c("eau_id", "fmv_per_ha_2017")
names(fmv_count) <- c("eau_id", "n_pixels")

fmv_by_eau <- left_join(fmv_mean, fmv_count, by = "eau_id") %>%
  filter(eau_id %in% eau_wmd$eau_id)

# Identify edge EAUs: those with < 50% of the median pixel count.
# Interior EAUs fill their full rectangle; edge EAUs that are clipped
# (by the WMD mask or study extent) will have fewer pixels.
typical_n_pixels <- median(fmv_by_eau$n_pixels, na.rm = TRUE)
low_coverage     <- fmv_by_eau %>% filter(n_pixels < 0.5 * typical_n_pixels)

cat(sprintf("  Typical pixels per EAU (median): %d\n", round(typical_n_pixels)))
cat(sprintf("  EAUs with < 50%% pixel coverage:  %d\n", nrow(low_coverage)))

if (nrow(low_coverage) > 0) {
  edge_by_wmd <- low_coverage %>%
    left_join(eau_wmd %>% select(eau_id, wmd_id), by = "eau_id") %>%
    count(wmd_id, name = "n_edge_eaus")
  cat("  Edge EAU count by WMD:\n")
  print(edge_by_wmd)
}

cat("  EAUs with cost data: ", sum(!is.na(fmv_by_eau$fmv_per_ha_2017)), "\n")
cat("  EAUs missing cost:   ", sum(is.na(fmv_by_eau$fmv_per_ha_2017)), "\n")
cat(sprintf("  FMV range ($/ha):    $%s  to  $%s\n",
            format(round(min(fmv_by_eau$fmv_per_ha_2017, na.rm = TRUE)),
                   big.mark = ","),
            format(round(max(fmv_by_eau$fmv_per_ha_2017, na.rm = TRUE)),
                   big.mark = ",")))


# ── 8. Compute EAU area and baseline total cost (2017 USD) ────────────────────####
# EAUs are equal-area by construction; area is fully determined by the raster
# cell size. No need to use area_km2 from eau_wmd (which records WMD-level area).
# 1 hectare = 10,000 m²

eau_area_m2 <- prod(res(eau_r))        # m²  (res[1] × res[2])
eau_area_ha <- eau_area_m2 / 1e4       # ha
eau_area_km2 <- eau_area_m2 / 1e6      # km² (for labelling)

cat(sprintf("\nEAU area: %.2f km² = %.0f ha\n", eau_area_km2, eau_area_ha))

fmv_by_eau <- fmv_by_eau %>%
  mutate(cost_2017_usd = fmv_per_ha_2017 * eau_area_ha)

cat(sprintf("  Total cost range (2017 USD):  $%s  to  $%s\n",
            format(round(min(fmv_by_eau$cost_2017_usd, na.rm = TRUE)),
                   big.mark = ","),
            format(round(max(fmv_by_eau$cost_2017_usd, na.rm = TRUE)),
                   big.mark = ",")))
cat(sprintf("  Mean total cost (2017 USD):   $%s\n",
            format(round(mean(fmv_by_eau$cost_2017_usd, na.rm = TRUE)),
                   big.mark = ",")))


# ── 9. Project costs to all decision years ────────────────────────────────────####
# cost[i, t] = cost_2017[i] × (1 + r)^(t − base_year)
#
# base_year = 2017 (dollar denomination year of PLACES data)
# r         = 0.02 (annual CPI inflation rate)
#
# Cost is scenario-invariant: the same value is assigned to all rcp/gcm rows
# within a given eau_id × year group. This is appropriate because land market
# prices are driven by local economic conditions, not by climate scenario.

decision_years <- sort(unique(data_panel$year))   # 2020, 2030, ..., 2100

cost_table <- fmv_by_eau %>% #note that even "Stationary" costs here are affected by inflation
  select(eau_id, cost_2017_usd) %>%
  crossing(year = decision_years) %>%
  mutate(
    inflation_multiplier = (1 + inflation_rate)^(year - base_year),
    cost                 = cost_2017_usd * inflation_multiplier
  ) %>%
  select(eau_id, year, cost)

cat("\nProjected mean total acquisition cost by decision year:\n\n")
cost_table %>%
  group_by(year) %>%
  summarise(
    mean_cost   = round(mean(cost)),
    median_cost = round(median(cost)),
    min_cost    = round(min(cost)),
    max_cost    = round(max(cost)),
    .groups = "drop"
  ) %>%
  print(n = Inf)


# ── 10. Join cost into data_panel ─────────────────────────────────────────────####
# Join on eau_id × year. The same cost value propagates across all rcp/gcm rows
# for each eau_id × year — this is intentional (cost is scenario-invariant).

n_rows_before <- nrow(data_panel)

data_panel <- data_panel %>%
  select(-cost) %>%
  left_join(cost_table, by = c("eau_id", "year"))

cat("\nJoined cost into data_panel.\n")
cat("  Rows before join: ", n_rows_before, "\n")
cat("  Rows after join:  ", nrow(data_panel), "\n")


# ── 11. Logic checks ──────────────────────────────────────────────────────────####

cat("\n========================================\n")
cat("  LOGIC CHECK: cost data\n")
cat("========================================\n")

# A1: No NAs in cost (every EAU × year combination should have a value)
n_na_cost <- sum(is.na(data_panel$cost))

# A2: All cost values are strictly positive
n_nonpos <- sum(data_panel$cost <= 0, na.rm = TRUE)

# A3: Cost is scenario-invariant within each eau_id × year group
#     (all 5 rcp/gcm rows should carry the same cost value)
cost_scenario_var <- data_panel %>%
  group_by(eau_id, year) %>%
  summarise(n_unique_cost = n_distinct(cost), .groups = "drop") %>%
  filter(n_unique_cost > 1)

# A4: Cost increases strictly over time for every EAU
#     (2% compounding means each decade's cost exceeds the previous)
cost_not_increasing <- data_panel %>%
  distinct(eau_id, year, cost) %>%
  arrange(eau_id, year) %>%
  group_by(eau_id) %>%
  mutate(cost_prev = lag(cost)) %>%
  filter(!is.na(cost_prev), cost <= cost_prev) %>%
  ungroup()

# A5: Row count is unchanged after the cost join
n_rows_expected <- n_rows_before

checks <- list(
  "No NAs in cost"                               = n_na_cost == 0,
  "All cost values strictly positive"            = n_nonpos == 0,
  "Cost scenario-invariant within eau_id × year" = nrow(cost_scenario_var) == 0,
  "Cost increases monotonically over time"       = nrow(cost_not_increasing) == 0,
  "Row count unchanged after cost join"          = nrow(data_panel) == n_rows_expected
)

for (nm in names(checks)) {
  cat(sprintf("  %s  %s\n", if (checks[[nm]]) "[PASS]" else "[FAIL]", nm))
}

failures <- names(checks)[!unlist(checks)]

if (length(failures) > 0) {
  
  cat("\n  --- Diagnostic detail ---\n")
  
  if (!checks[["No NAs in cost"]]) {
    cat("  NA count in cost:", n_na_cost, "\n")
    # Report which EAUs are missing cost data
    missing_eaus <- data_panel %>%
      filter(is.na(cost)) %>%
      distinct(eau_id, wmd_id)
    cat("  EAUs with NA cost:", nrow(missing_eaus), "\n")
    print(head(missing_eaus, 20))
  }
  if (!checks[["All cost values strictly positive"]]) {
    cat("  Non-positive cost count:", n_nonpos, "\n")
    print(data_panel %>% filter(cost <= 0) %>% distinct(eau_id, wmd_id, year, cost))
  }
  if (!checks[["Cost scenario-invariant within eau_id × year"]]) {
    cat("  Groups with varying cost:", nrow(cost_scenario_var), "\n")
    print(head(cost_scenario_var, 10))
  }
  if (!checks[["Cost increases monotonically over time"]]) {
    cat("  EAU-year pairs where cost did not increase:", nrow(cost_not_increasing), "\n")
    print(head(cost_not_increasing, 10))
  }
  if (!checks[["Row count unchanged after cost join"]]) {
    cat("  Expected:", n_rows_expected, "| Found:", nrow(data_panel), "\n")
  }
  
  cat("========================================\n\n")
  stop("Logic check FAILED: cost data has errors. ",
       "Investigate before proceeding to ILP scripts.\n",
       "Failed checks: ", paste(failures, collapse = ", "))
  
} else {
  cat(sprintf(
    "\n  All checks passed. EAUs: %d | Years: %d | Total rows: %d\n",
    n_distinct(data_panel$eau_id),
    n_distinct(data_panel$year),
    nrow(data_panel)
  ))
  cat(sprintf(
    "  Inflation: %.0f%% annual | Base year: %d | Area: %.0f ha per EAU\n",
    inflation_rate * 100, base_year, eau_area_ha
  ))
  cat("========================================\n\n")
}


# ── 12. Summary statistics ────────────────────────────────────────────────────####

cat("--- Cost distribution by WMD (2020 baseline, USD) ---\n\n")

data_panel %>%
  filter(year == 2020, rcp == "baseline") %>%
  group_by(wmd_id) %>%
  summarise(
    n_eaus      = n_distinct(eau_id),
    mean_cost   = round(mean(cost)),
    median_cost = round(median(cost)),
    min_cost    = round(min(cost)),
    max_cost    = round(max(cost)),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_cost)) %>%
  print(n = Inf)


cat("--- Cost distribution total (2020 baseline, USD) ---\n\n")

data_panel %>%
  filter(year == 2020, rcp == "baseline") %>%
  summarise(
    n_eaus      = n_distinct(eau_id),
    mean_cost   = round(mean(cost)),
    median_cost = round(median(cost)),
    min_cost    = round(min(cost)),
    max_cost    = round(max(cost)),
    total_cost  = round(sum(cost)),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_cost)) %>%
  print(n = Inf)


# ── 13. Save ──────────────────────────────────────────────────────────────────####

saveRDS(data_panel,  "derived_data/data_panel.rds")
write_csv(data_panel, "derived_data/data_panel.csv")

cat("\n✓ data_panel saved with cost column populated.\n")
cat(sprintf(
  "  Source: PLACES FMV: vacant (Nolte 2020) | Base: 2017 USD | Inflation: %.0f%%/yr\n",
  inflation_rate * 100
))
