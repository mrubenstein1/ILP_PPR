###Construct Data Panel for Analysis ###
# Build a data panel to contain data needed for model analysis
    # takes the outputs of EAU_prop_suitable (the datatable of percent suitable 
    # habitat per year per EAU), and uses it to construct the final data table. 
    # The final data panel (i.e., final product of data import/manipulation) will be a 
    # dataframe of EAU per time period per RCP per GCM.

# Note that only LC data from 2030 onwards is used from the 4.5 or 8.5 scenarios. 
    # LC estimates for 2020 will come from the "mean clim" analysis

#------ Setup ________
if (!exists(".SETUP_DONE")) source("00_setup.R")
 
# 1. Load EAU–WMD crosswalk ####
eau_wmd <- read_csv(file.path(DIR_COMMON, "eau_wmd_lookup.csv"), show_col_types = FALSE)

 
# 2. Load suitable habitat tables for both RCPs ####
lu_45 <- read_csv(file.path(DIR_DERIVED, "Lu_prop_45.csv"), show_col_types = FALSE)
lu_85 <- read_csv(file.path(DIR_DERIVED, "Lu_prop_85.csv"), show_col_types = FALSE)
keep45 <- readRDS(file.path(DIR_DERIVED, "Lu_prop_keep_rows_45.rds"))
keep85 <- readRDS(file.path(DIR_DERIVED, "Lu_prop_keep_rows_85.rds"))
keepmc <- readRDS(file.path(DIR_DERIVED, "Lu_prop_keep_rows_meanclim.rds"))
stopifnot(identical(keep45, keep85), identical(keep45, keepmc))
lu_cells <- which(keep45)
stopifnot(nrow(lu_45) == length(lu_cells), nrow(lu_85) == length(lu_cells))



# Columns look like rcp45_2014, rcp45_2020, ...; similarly for rcp85_* [file:1]

 
# 3. Identify common years (2030–2100) across both RCP tables ####
 
cols_45 <- grep("^rcp45_[0-9]{4}$", names(lu_45), value = TRUE)
cols_85 <- grep("^rcp85_[0-9]{4}$", names(lu_85), value = TRUE)

years_45 <- as.integer(sub("rcp45_", "", cols_45))
years_85 <- as.integer(sub("rcp85_", "", cols_85))

years_common <- sort(intersect(years_45, years_85))
years_use <- years_common[years_common >= 2030 & years_common <= 2100]

if (length(years_use) == 0) {
  stop("No overlapping years between 2030 and 2100 found in Lu_prop_45 and Lu_prop_85.")
}

cols_45_use <- paste0("rcp45_", years_use)
cols_85_use <- paste0("rcp85_", years_use)

 
# 4. Build long table for RCP 4.5: EAU × year ####
 
panel_45 <- lu_45 %>%
  select(all_of(cols_45_use)) %>%
  mutate(cell_id = lu_cells) %>%
  pivot_longer(
    cols      = starts_with("rcp45_"),
    names_to  = "time_period",
    values_to = "prop_suitable"
  ) %>%
  mutate(
    year         = as.integer(str_extract(time_period, "[0-9]{4}")),
    rcp          = "45"
  ) %>%
  select(cell_id, year, rcp, prop_suitable)

 
# 5. Build long table for RCP 8.5: EAU × year ####
 
panel_85 <- lu_85 %>%
  select(all_of(cols_85_use)) %>%
  mutate(cell_id = lu_cells)  %>%
  pivot_longer(
    cols      = starts_with("rcp85_"),
    names_to  = "time_period",
    values_to = "prop_suitable"
  ) %>%
  mutate(
    year         = as.integer(str_extract(time_period, "[0-9]{4}")),
    rcp          = "85"
  ) %>%
  select(cell_id, year, rcp, prop_suitable)

 
# 6. Stack both RCPs ####
 
panel_suitable <- bind_rows(panel_45, panel_85)

# At this point: one row per "retained EAU row" × year × RCP [file:1]



# 7. Attach EAU IDs and WMD IDs — BY PLACE (raster cell number) ####

eau_meta <- eau_wmd %>% select(cell_id, eau_id, wmd_id)

# Log both directions of the mismatch before joining.
cells_no_eau <- setdiff(lu_cells, eau_meta$cell_id)   # land cover outside any WMD
eaus_no_data <- setdiff(eau_meta$cell_id, lu_cells)   # EAU with no FOREsce coverage

cat(sprintf("  Land-cover cells not in any WMD : %d\n", length(cells_no_eau)))
cat(sprintf("  EAUs with no land-cover data    : %d   (expect 41)\n", length(eaus_no_data)))
if (length(eaus_no_data)) {
  cat("  Dropped EAUs by WMD:\n")
  eau_wmd %>% filter(cell_id %in% eaus_no_data) %>%
    count(wmd_id, name = "n_dropped") %>% arrange(desc(n_dropped)) %>%
    as.data.frame() %>% print(row.names = FALSE)
}

panel_suitable <- panel_suitable %>% inner_join(eau_meta, by = "cell_id")
stopifnot(n_distinct(panel_suitable$eau_id) == 1169L)

 
# 8. Duplicate rows for wet vs dry GCMs. ####
    #Note: THIS WOULD NEED TO CHANGE IF WE FIND A DIFFERENT LC SCENARIO TO MATCH THE DIFFERENT GCMS
 
gcm_levels <- c("wet", "dry")

panel_suitable_gcm <- panel_suitable %>%
  tidyr::crossing(gcm = gcm_levels)

 
# 9. Add empty placeholders for densities, transitions, and cost ####
 
eau_panel <- panel_suitable_gcm %>%
  mutate(
    abs_abundance   = NA_real_,  # to be filled later with WMD-level data
    scaled_abundance = NA_real_, # will depend on abs_abundance & prop_suitable
    trans_prob    = NA_real_,  # placeholder
    cost          = NA_real_   # placeholder
  ) %>%
  select(
    eau_id,
    wmd_id,
    year,
    rcp,
    gcm,
    prop_suitable,
    abs_abundance,
    scaled_abundance,
    trans_prob,
    cost,
    everything()
  )


 
# 10. Add rows for baseline ("meanclim") data from 2020 ####
#Note that LC data for the 2020 time step comes from the "mean climate" GCAM reference scenario,
  #unlike the 2030-2100 data which comes from either RCP4.5 or 8.5. 
#There will only be one 2020 row per EAU, since there is no applicable RCP or GCM (wet/dry) combinations in this timestep
 

# 1. Load meanclim 2020 LC data and attach EAU IDs via row-order alignment
#    Row order aligns with lu_45/lu_85 because keep_rows vectors are identical across all scenarios.
lu_meanclim_2020 <- read_csv(file.path(DIR_DERIVED, "Lu_prop_meanclim_2020.csv"),
                             show_col_types = FALSE)

baseline_2020 <- lu_meanclim_2020 %>%
  mutate(cell_id = lu_cells) %>%
  inner_join(eau_meta, by = "cell_id") %>%
  select(eau_id, prop_suitable_baseline = meanclim_2020)


# 2. Build one 2020 row per EAU.
#    RCP/GCM columns will be filled with "baseline/baseline" 
rows_2020 <- baseline_2020 %>%
  left_join(eau_wmd %>% select(eau_id, wmd_id), by = "eau_id") %>%
  mutate(
    year             = 2020L,
    rcp              = "baseline",
    gcm              = "baseline",
    prop_suitable    = prop_suitable_baseline,
    abs_abundance    = NA_real_,
    scaled_abundance = NA_real_,
    trans_prob       = NA_real_,
    cost             = NA_real_
  ) %>%
  select(eau_id, wmd_id, year, rcp, gcm, prop_suitable,
         abs_abundance, scaled_abundance, trans_prob, cost)

# 3. Bind 2020 rows into the panel
eau_panel <- bind_rows(eau_panel, rows_2020)


# 4. Add 5th row per EAU-Timestep combo for a "stationary" scenario (from 2030-2100)

# Identify the key columns that define an EAU–time step
key_cols <- c("eau_id", "wmd_id", "year")

# Create a stationary row for each EAU–year
stationary_rows <- eau_panel %>%
  distinct(across(all_of(key_cols))) %>%
  filter(year >= 2030) %>%
  mutate(
    rcp = "stationary",
    gcm = "stationary"
  )

# Make sure stationary_rows has ALL columns from eau_panel
#    For any column that's missing, add it as NA of the correct type
for (col in setdiff(names(eau_panel), names(stationary_rows))) {
  # Infer type from the existing column in eau_panel
  prototype <- eau_panel[[col]]
  if (is.numeric(prototype)) {
    stationary_rows[[col]] <- NA_real_
  } else if (is.character(prototype)) {
    stationary_rows[[col]] <- NA_character_
  } else if (is.logical(prototype)) {
    stationary_rows[[col]] <- NA
  } else {
    # Fallback: NA of same class
    stationary_rows[[col]] <- NA
  }
}

# Reorder columns to match eau_panel
stationary_rows <- stationary_rows %>%
  select(all_of(names(eau_panel)))

# Bind stationary rows to the existing panel
eau_panel <- bind_rows(
  eau_panel,
  stationary_rows
)

# 5. Copy 2020 land cover data to all "stationary" rows. 

eau_panel <- eau_panel %>%
  left_join(baseline_2020, by = "eau_id") %>%
  mutate(
    prop_suitable = if_else(
      rcp == "stationary" & year > 2020,
      prop_suitable_baseline,
      prop_suitable
    )
  ) %>%
  select(-prop_suitable_baseline)





 
# 10b. Attach gross habitat flows (loss and gain) — BY PLACE ####
#
# Script 05 needs the GROSS conversion rate, L / prop_suitable, where L is the
# area converting suitable -> unsuitable in the decade. Differencing
# prop_suitable gives only the NET change, (L - G) / prop_suitable, which
# conceals  conversion across the modelled EAUs and inverts the sign of the hazard in four mid-century
# decade-RCP cells. L and G are therefore measured at native FOREsce resolution
# in script 02 section 5 and read here.
#
# WHY A SEPARATE OBJECT AND NOT A COLUMN ON eau_panel:
# A flow row describes a TRANSITION, year -> year + 10, so the flows run
# 2020..2090 while the panel's rows run 2020..2100. They do not align. In
# particular the 2020 flows exist separately for RCP 4.5 and 8.5, while the
# panel carries a single 2020 row per EAU with rcp = "baseline" (script 05
# pools the two into it). Forcing these onto the panel would mean either
# inventing 2020 RCP rows -- which breaks the row-count checks below and every
# downstream assumption about panel shape -- or widening to prop_lost_45 /
# prop_lost_85 columns that are NA on almost every row. Keeping flows in their
# own long table costs script 05 one extra read and keeps the panel schema
# untouched.
#
# What this section DOES own is the cell_id -> eau_id translation, because that
# is this script's job and the by-place join is verified here (check 7 below).
# Script 05 then works purely in eau_id and never sees a cell number.

flows_45 <- read_csv(file.path(DIR_DERIVED, "Lu_flows_45.csv"), show_col_types = FALSE)
flows_85 <- read_csv(file.path(DIR_DERIVED, "Lu_flows_85.csv"), show_col_types = FALSE)

eau_flows <- bind_rows(
  flows_45 %>% mutate(rcp = "45"),
  flows_85 %>% mutate(rcp = "85")
) %>%
  # read_csv parses whole numbers as double. Coerce the two key columns to
  # integer: cell_id must match eau_meta's type for the join below, and year must
  # match the panel's integer year for the joins script 05 does against this
  # table. Leaving them double works by coercion in dplyr but breaks any
  # identical() comparison, including the year check further down.
  mutate(cell_id = as.integer(cell_id), year = as.integer(year)) %>%
  inner_join(eau_meta, by = "cell_id") %>%      # by raster cell number, as above
  filter(!is.na(prop_lost), !is.na(prop_gained),
         !is.na(prop_lost_unprot), !is.na(prop_suit_unprot)) %>%
  select(eau_id, wmd_id, year, rcp, prop_lost, prop_gained,
         prop_lost_unprot, prop_suit_unprot) %>%
  arrange(rcp, year, eau_id)

cat(sprintf("\n  Gross flows loaded: %d rows | %d EAUs | years %d-%d | RCPs %s\n",
            nrow(eau_flows), n_distinct(eau_flows$eau_id),
            min(eau_flows$year), max(eau_flows$year),
            paste(sort(unique(eau_flows$rcp)), collapse = ", ")))

# The flow table must cover exactly the EAUs the panel covers, and every
# non-terminal period, or script 05 will silently produce NA hazards. Checking
# here rather than in 05 keeps the failure next to the join that caused it.
flow_years    <- sort(unique(eau_flows$year))
expected_flow_years <- as.integer(c(2020L, head(years_use, -1)))  # 2020..2090
flow_eaus     <- sort(unique(eau_flows$eau_id))
panel_eaus    <- sort(unique(eau_panel$eau_id))
flow_complete <- eau_flows %>% count(year, rcp, name = "n_eau")

# setequal() rather than identical() on the set comparisons: identical() is
# strict about storage mode, so an integer-vs-double difference in a column read
# from CSV fails it even when every value matches.
flow_checks <- list(
  "Flow EAU set matches panel EAU set" =
    setequal(flow_eaus, panel_eaus) && length(flow_eaus) == length(panel_eaus),
  "Flow years are 2020..2090 (2100 terminal, no flow)" =
    setequal(flow_years, expected_flow_years) &&
    length(flow_years) == length(expected_flow_years),
  "Every year x RCP has the full EAU set" =
    all(flow_complete$n_eau == length(panel_eaus)),
  "Flows are non-negative" =
    min(eau_flows$prop_lost) >= 0 && min(eau_flows$prop_gained) >= 0 &&
    min(eau_flows$prop_lost_unprot) >= 0 && min(eau_flows$prop_suit_unprot) >= 0
)

for (nm in names(flow_checks))
  cat(sprintf("  %s  %s\n", if (flow_checks[[nm]]) "[PASS]" else "[FAIL]", nm))

if (any(!unlist(flow_checks))) {
  cat("\n  --- Diagnostic detail ---\n")
  cat(sprintf("  Flow EAUs: %d | panel EAUs: %d | in flows only: %d | in panel only: %d\n",
              length(flow_eaus), length(panel_eaus),
              length(setdiff(flow_eaus, panel_eaus)),
              length(setdiff(panel_eaus, flow_eaus))))
  cat(sprintf("  Flow years: %s\n", paste(flow_years, collapse = ", ")))
  cat("  EAU count per year x RCP:\n")
  print(as.data.frame(flow_complete), row.names = FALSE)
  stop("Gross flow table does not line up with the panel. Re-run script 02 ",
       "section 5 before proceeding.\n  Failed: ",
       paste(names(flow_checks)[!unlist(flow_checks)], collapse = ", "))
}

saveRDS(eau_flows, file.path(DIR_DERIVED, "eau_flows.rds"))


# 11. Save panel for later steps ####
 
saveRDS(eau_panel,  file.path(DIR_DERIVED, "eau_panel.rds"))
write_csv(eau_panel, file.path(DIR_DERIVED, "eau_panel.csv"))

 
# 12. Logic Check ####
 

# Pre-compute expected dimensions from upstream objects
n_future_years         <- length(years_use)            # number of 2030–2100 time steps
expected_years_per_eau <- n_future_years + 1L          # +1 for 2020
expected_rows_per_eau  <- (5L * n_future_years) + 1L   # 5 rows per future year, 1 for 2020

# ── Pre-compute all diagnostics before running checks ───────────────────────

# Check 1: 2020 should have 1 row per EAU; future years should have 5
rows_per_eau_year <- eau_panel %>%
  count(eau_id, year, name = "n_rows")

bad_2020_rows   <- rows_per_eau_year %>% filter(year == 2020, n_rows != 1)
bad_future_rows <- rows_per_eau_year %>% filter(year >= 2030, n_rows != 5)

# Check 2: Each EAU should have the correct number of distinct years
years_per_eau <- eau_panel %>%
  distinct(eau_id, year) %>%
  count(eau_id, name = "n_years")

bad_year_counts <- years_per_eau %>% filter(n_years != expected_years_per_eau)

# Check 3: Each EAU should have the correct total row count
rows_per_eau <- eau_panel %>%
  count(eau_id, name = "n_rows_eau")

bad_row_counts <- rows_per_eau %>% filter(n_rows_eau != expected_rows_per_eau)

# Check 4: EAU counts per WMD should match the lookup table
eaus_per_wmd_panel <- eau_panel %>%
  distinct(eau_id, wmd_id) %>%
  count(wmd_id, name = "n_eaus_panel")

eaus_per_wmd_lookup <- eau_wmd %>%
  filter(cell_id %in% lu_cells) %>%
  count(wmd_id, name = "n_eaus_lookup")

bad_wmd_counts <- eaus_per_wmd_lookup %>%
  left_join(eaus_per_wmd_panel, by = "wmd_id") %>%
  filter(n_eaus_lookup != n_eaus_panel | is.na(n_eaus_panel))

# Check 5: prop_suitable should be identical across wet/dry within each EAU–year–RCP
suitable_violations <- eau_panel %>%
  group_by(eau_id, year, rcp) %>%
  summarise(n_unique_prop = n_distinct(prop_suitable), .groups = "drop") %>%
  filter(n_unique_prop > 1)

# Check 6: prop_suitable should be constant over time for stationary scenario
stationary_violations <- eau_panel %>%
  filter(rcp == "stationary") %>%
  group_by(eau_id) %>%
  summarise(
    min_prop = min(prop_suitable, na.rm = TRUE),
    max_prop = max(prop_suitable, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!dplyr::near(min_prop, max_prop))

# Check 7: verify that spatial join has completed properly
verify_by_place <- function(lu_tbl, rcp_label) {
  truth_key <- eau_meta %>% select(cell_id, eau_id)
  worst <- 0
  for (y in intersect(years_use, unique(eau_panel$year))) {
    truth <- tibble(cell_id = lu_cells,
                    truth = lu_tbl[[paste0("rcp", rcp_label, "_", y)]]) %>%
      inner_join(truth_key, by = "cell_id")
    stored <- eau_panel %>% filter(year == y, rcp == rcp_label) %>%
      distinct(eau_id, prop_suitable)
    cmp <- inner_join(truth, stored, by = "eau_id")
    if (nrow(cmp) != nrow(stored)) return(list(ok = FALSE, worst = NA_real_))
    worst <- max(worst, max(abs(cmp$truth - cmp$prop_suitable)))
  }
  list(ok = worst <= 1e-9, worst = worst)
}
v45 <- verify_by_place(lu_45, "45"); v85 <- verify_by_place(lu_85, "85")
cat(sprintf("  %s prop_suitable on own EAU, RCP 4.5 (max dev %.2e)\n",
            if (v45$ok) "[PASS]" else "[FAIL]", v45$worst))
cat(sprintf("  %s prop_suitable on own EAU, RCP 8.5 (max dev %.2e)\n",
            if (v85$ok) "[PASS]" else "[FAIL]", v85$worst))
if (!v45$ok || !v85$ok) stop("BY-PLACE CHECK FAILED — habitat values misassigned.")

# ── Run all checks ───────────────────────────────────────────────────────────
checks <- list(
  "2020 has exactly 1 row per EAU"                   = nrow(bad_2020_rows) == 0,
  "Future years have exactly 5 rows per EAU"         = nrow(bad_future_rows) == 0,
  "Each EAU has correct number of time periods"      = nrow(bad_year_counts) == 0,
  "Each EAU has correct total row count"             = nrow(bad_row_counts) == 0,
  "EAU counts per WMD match lookup table"            = nrow(bad_wmd_counts) == 0,
  "prop_suitable identical across wet/dry GCMs"      = nrow(suitable_violations) == 0,
  "Stationary prop_suitable constant over time"      = nrow(stationary_violations) == 0
)

for (nm in names(checks)) {
  cat(sprintf("  %s  %s\n", if (checks[[nm]]) "[PASS]" else "[FAIL]", nm))
}

failures <- names(checks)[!unlist(checks)]

if (length(failures) > 0) {
  
  cat("\n  --- Diagnostic detail ---\n")
  
  if (!checks[["2020 has exactly 1 row per EAU"]]) {
    cat("  EAU-years with wrong 2020 row count:", nrow(bad_2020_rows), "\n")
    print(head(bad_2020_rows, 10))
  }
  if (!checks[["Future years have exactly 5 rows per EAU"]]) {
    cat("  EAU-years with wrong future row count:", nrow(bad_future_rows), "\n")
    print(head(bad_future_rows, 10))
  }
  if (!checks[["Each EAU has correct number of time periods"]]) {
    cat("  Expected years per EAU:", expected_years_per_eau, "\n")
    cat("  EAUs with wrong year count:", nrow(bad_year_counts), "\n")
    print(bad_year_counts %>% count(n_years))
  }
  if (!checks[["Each EAU has correct total row count"]]) {
    cat("  Expected rows per EAU:", expected_rows_per_eau, "\n")
    cat("  EAUs with wrong row count:", nrow(bad_row_counts), "\n")
    print(bad_row_counts %>% count(n_rows_eau))
  }
  if (!checks[["EAU counts per WMD match lookup table"]]) {
    cat("  WMDs with mismatched EAU counts:\n")
    print(bad_wmd_counts)
  }
  if (!checks[["prop_suitable identical across wet/dry GCMs"]]) {
    cat("  EAU-year-RCP combos where wet != dry:", nrow(suitable_violations), "\n")
    print(head(suitable_violations, 10))
  }
  if (!checks[["Stationary prop_suitable constant over time"]]) {
    cat("  EAUs where stationary prop_suitable varies over time:", nrow(stationary_violations), "\n")
    print(head(stationary_violations, 10))
  }
  
  cat("========================================\n\n")
  stop("Logic check FAILED: eau_panel has structural errors. ",
       "Investigate before proceeding to downstream analysis.\n",
       "Failed checks: ", paste(failures, collapse = ", "))
  
} else {
  cat(sprintf(
    "\n  All checks passed. EAUs: %d | Years: %d | Total rows: %d\n",
    n_distinct(eau_panel$eau_id),
    n_distinct(eau_panel$year),
    nrow(eau_panel)
  ))
  cat("========================================\n\n")
}
