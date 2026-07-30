
##### Benefit (Duck abundance) Data ####
    # Import and clean benefit data
    # add to the established data panel

#------ Setup ________
if (!exists(".SETUP_DONE")) source("00_setup.R")

#### 1. Import data #####
# ── 0. Load inputs ────────────────────────────────────────────────────────────
eau_panel <- readRDS("derived_data/eau_panel.rds")
benefit <- read.csv("input_data/benefit.csv")

# restructure around gcm/rcp/time step structure.
# Data provided (McKenna et al, 2026) reflects estimates for 3 periods: observation; midcentury; 
  #and end of century. Each period actually represents a ~30 year period; we have chosen to anchor that estimate
  #to the mid-year of each window in order to map to our decision timeframe. 

# anchor year: Obs (2007-2017) -> 2012, "3059" (2030-2059) -> 2045, "7099" (2070-2099) -> 2085.
benefit_new <- benefit %>%
  mutate(
    year = case_when(
      scenario == "Obs" ~ 2012L,             # observation window midpoint
      str_detect(scenario, "3059") ~ 2045L,  # mid-century window midpoint
      str_detect(scenario, "7099") ~ 2085L,  # end-century window midpoint
      TRUE ~ NA_integer_
    ),
    rcp = case_when(
      str_detect(scenario, "4.5") ~ "45",
      str_detect(scenario, "8.5") ~ "85",
      TRUE ~ NA_character_
    ),
    gcm = case_when(
      str_detect(str_to_lower(scenario), "dry") ~ "dry",
      str_detect(str_to_lower(scenario), "wet") ~ "wet",
      TRUE ~ NA_character_
    ),
    wmd_id = WMD,
    ben    = total_sum
  ) %>%
  select(wmd_id, year, rcp, gcm, ben)


#### 2. Interpolate Decadal Data ####
# One linear interpolation per WMD × RCP × GCM over the three window-midpoint
# anchors (2012, 2045, 2085), evaluated at the decision years 2020-2100 (by 10).
# The Obs anchor is SHARED across combos (one observed value per WMD); the mid/end
# anchors are combo-specific, so the combos diverge over the horizon. `rule = 2`
# holds the last (2085) anchor value constant for decision years beyond it, so
# 2090 and 2100 sit at the 2085 value.

decision_years <- as.integer(seq(2020, 2100, by = 10))

# ── Step 1: Extract the single observed anchor per WMD (rcp/gcm are NA at Obs) ─
ben_obs <- benefit_new %>%
  filter(is.na(rcp)) %>%
  select(wmd_id, ben)

# ── Step 2: Extract combo-specific future anchors (mid + end) ─────────────────
future_data <- benefit_new %>%
  filter(!is.na(rcp)) %>%
  select(wmd_id, rcp, gcm, year, ben)

# ── Step 3: Build full anchor set for interpolation ───────────────────────────
# Attach the shared observed value to every RCP × GCM combo at the Obs anchor year
anchor_obs <- future_data %>%
  distinct(wmd_id, rcp, gcm) %>%
  left_join(ben_obs, by = "wmd_id") %>%
  mutate(year = 2012L)

anchors <- bind_rows(anchor_obs, future_data) %>%
  arrange(wmd_id, rcp, gcm, year)          # approx() needs x sorted within group

# ── Step 4: Interpolate each combo at the decision years (flat tail via rule=2) ─
benefit_decadal_future <- anchors %>%
  group_by(wmd_id, rcp, gcm) %>%
  group_modify(~ tibble(
    year = decision_years,
    ben  = approx(x = .x$year, y = .x$ben, xout = decision_years, rule = 2)$y
  )) %>%
  ungroup()

# ── Step 5: Shared 2020 baseline = mean across the four combos at 2020 ────────
# 2020 is a modelled decision baseline: each combo is interpolated on the
# 2012->2045 segment, then averaged into one shared period-0 row per WMD
# (rcp/gcm = NA at baseline, matching the join convention used downstream).
ben_2020_shared <- benefit_decadal_future %>%
  filter(year == 2020L) %>%
  group_by(wmd_id) %>%
  summarise(ben = mean(ben), .groups = "drop") %>%
  mutate(year = 2020L)

# ── Step 6: Bind the shared 2020 baseline back onto the future decadal rows ───
benefit_decadal <- bind_rows(
  ben_2020_shared,
  benefit_decadal_future %>% filter(year > 2020)   # combo-specific 2030..2100
) %>%
  arrange(wmd_id, year, rcp, gcm)


##### 3. Match Benefit data into EAU_Panel ####

# ── WMD exclusions ───────────────────────────────────────────────────────────
# Hydrological: the data underpinning the duck abundance estimates performed
# poorly in these three Montana districts.
# Coverage: districts where >20% of EAUs fall outside the FOREsce extent are
# excluded, because allocating the district's full abundance across a heavily
# truncated denominator inflates apparent parcel-level density.
COVERAGE_MAX_MISSING_PCT <- 20

eau_wmd <- read_csv(file.path(DIR_DERIVED, "eau_wmd_lookup.csv"), show_col_types = FALSE)

coverage <- eau_wmd %>%
  mutate(has_data = eau_id %in% unique(eau_panel$eau_id)) %>%
  group_by(wmd_id) %>%
  summarise(n_total = n(), n_ok = sum(has_data),
            pct_missing = 100 * (1 - n_ok / n_total), .groups = "drop") %>%
  arrange(desc(pct_missing))

wmd_exclude_hydro    <- c("Benton Lake", "Bowdoin", "Northeast Montana")
wmd_exclude_coverage <- coverage %>%
  filter(pct_missing > COVERAGE_MAX_MISSING_PCT) %>% pull(wmd_id)
wmd_exclude <- union(wmd_exclude_hydro, wmd_exclude_coverage)

cat("  Excluded (hydrological):", paste(wmd_exclude_hydro, collapse = ", "), "\n")
cat("  Excluded (coverage >", COVERAGE_MAX_MISSING_PCT, "% missing):",
    paste(wmd_exclude_coverage, collapse = ", "), "\n")
print(coverage %>% filter(pct_missing > 0) %>% as.data.frame(), row.names = FALSE)

stopifnot(setequal(wmd_exclude_coverage, "Fergus Falls"))

# Remove excluded WMDs from the base panel before joining benefit data
eau_panel <- eau_panel %>%
  filter(!wmd_id %in% wmd_exclude)

# ── Extract 2020 baseline (one row per WMD, used for 2020 and stationary) ────
ben_2020_for_join <- benefit_decadal %>%
  filter(year == 2020) %>%
  distinct(wmd_id, .keep_all = TRUE) %>%
  select(wmd_id, ben)

# ── 2020: join on wmd_id only (rcp/gcm are NA and meaningless at baseline) ───
panel_2020_with_ben <- eau_panel %>%
  filter(year == 2020) %>%
  left_join(ben_2020_for_join, by = "wmd_id") %>%
  mutate(abs_abundance = ben) %>%
  select(-ben)

# ── Future stationary: hold 2020 (modelled) abundance constant across all years ─
panel_future_stationary_with_ben <- eau_panel %>%
  filter(year > 2020, rcp == "stationary") %>%
  left_join(ben_2020_for_join, by = "wmd_id") %>%
  mutate(abs_abundance = ben) %>%
  select(-ben)

# ── Future non-stationary: join on all four keys ──────────────────────────────
ben_future_for_join <- benefit_decadal %>%
  filter(year > 2020) %>%
  select(wmd_id, year, rcp, gcm, ben)

panel_future_nonstationary_with_ben <- eau_panel %>%
  filter(year > 2020, rcp != "stationary") %>%
  left_join(ben_future_for_join, by = c("wmd_id", "year", "rcp", "gcm")) %>%
  mutate(abs_abundance = ben) %>%
  select(-ben)

# ── Recombine ─────────────────────────────────────────────────────────────────
eau_panel_with_ben <- bind_rows(
  panel_2020_with_ben,
  panel_future_stationary_with_ben,
  panel_future_nonstationary_with_ben
) %>%
  arrange(eau_id, year, rcp, gcm)

#### 4. Distribute benefit by proportional habitat #####
eau_panel_alloc <- eau_panel_with_ben %>%
  group_by(wmd_id, year, rcp, gcm) %>%
  mutate(
    total_prop_suitable = sum(prop_suitable, na.rm = TRUE),
    habitat_share = ifelse(
      total_prop_suitable > 0,
      prop_suitable / total_prop_suitable,
      NA_real_
    ),
    scaled_abundance = habitat_share * abs_abundance
  ) %>%
  ungroup()


#### 5. Rename final object for clarity ####

data_panel <- eau_panel_alloc

### 6. Logic check ####

# Number of EAUs after WMD exclusion, derived from the filtered base panel
n_eaus_expected <- n_distinct(eau_panel$eau_id)

"Expected EAU count after exclusions" = n_distinct(eau_panel$eau_id) == 879L
"Expected WMD count after exclusions" = n_distinct(eau_panel$wmd_id) == 20L

# Rows per EAU breaks down as follows:
#   - 1 baseline row for 2020 (rcp = "baseline", gcm = "baseline")
#   - For each future decade (2030, 2040, ..., 2100 = 8 time steps):
#       4 non-stationary scenario rows (rcp45 × wet, rcp45 × dry,
#                                       rcp85 × wet, rcp85 × dry)
#     + 1 stationary row (rcp = "stationary", gcm = "stationary")
#     = 5 rows per future time step
n_rcp_gcm_combos     <- 4L                              # RCP45/85 × wet/dry
n_stationary         <- 1L                              # stationary scenario
n_scenarios_per_step <- n_rcp_gcm_combos + n_stationary # 5 rows per future year
n_future_years       <- length(seq(2030, 2100, by = 10))# 8 decades
n_baseline_rows      <- 1L                              # 2020 row

rows_per_eau    <- (n_scenarios_per_step * n_future_years) + n_baseline_rows
n_rows_expected <- n_eaus_expected * rows_per_eau

cat(sprintf("  Expected EAUs:         %d\n", n_eaus_expected))
cat(sprintf("  Expected rows per EAU: %d  (%d scenarios x %d decades + %d baseline)\n",
            rows_per_eau, n_scenarios_per_step, n_future_years, n_baseline_rows))
cat(sprintf("  Expected total rows:   %d\n", n_rows_expected))
cat("\n")


# ── Pre-compute all diagnostics ──────────────────────────────────────────────

# A1: All original anchor rows present in the interpolation input with matching
#     values. Checked against `anchors` (the approx() input) rather than
#     benefit_decadal, because the window-midpoint anchor years (2012/2045/2085) are
#     not decision years and never appear in benefit_decadal directly. Exactness of
#     the interpolation AT the knots is then guaranteed by approx().
future_src       <- benefit_new %>% filter(!is.na(rcp)) %>%
  select(wmd_id, rcp, gcm, year, ben)
obs_src          <- benefit_new %>% filter(is.na(rcp)) %>% distinct(wmd_id, ben)

unmatched        <- anti_join(future_src, anchors,
                              by = c("wmd_id", "rcp", "gcm", "year"))

value_mismatches <- inner_join(future_src, anchors,
                               by     = c("wmd_id", "rcp", "gcm", "year"),
                               suffix = c("_new", "_decadal")) %>%
  filter(!near(ben_new, ben_decadal))

obs_mismatches   <- obs_src %>%
  left_join(ben_obs %>% rename(ben_anchor = ben), by = "wmd_id") %>%
  filter(is.na(ben_anchor) | !near(ben, ben_anchor))

# A2: benefit_decadal row structure
#     2020: 1 row per WMD (observational anchor; rcp/gcm = NA)
#     Future decades: 4 rows per WMD x year (rcp45/85 x wet/dry)
bad_2020_ben <- benefit_decadal %>%
  filter(year == 2020) %>%
  count(wmd_id) %>%
  filter(n != 1)

bad_future_ben <- benefit_decadal %>%
  filter(year > 2020) %>%
  count(wmd_id, year) %>%
  filter(n != n_rcp_gcm_combos)

# B1: No NAs in abs_abundance in eau_panel_with_ben (pre-downscaling)
n_na_abs <- sum(is.na(eau_panel_with_ben$abs_abundance))

# B2: abs_abundance constant across all EAUs within each WMD x year x scenario
abs_inconsistent <- eau_panel_with_ben %>%
  group_by(wmd_id, year, rcp, gcm) %>%
  summarise(
    min_abs = min(abs_abundance, na.rm = TRUE),
    max_abs = max(abs_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!near(min_abs, max_abs))

# B3: 2020 observed abundance matches stationary abundance for all WMDs
stationary_mismatch <- eau_panel_with_ben %>%
  filter(year == 2020) %>%
  distinct(wmd_id, abs_abundance) %>%
  rename(ben_2020 = abs_abundance) %>%
  left_join(
    eau_panel_with_ben %>%
      filter(rcp == "stationary") %>%
      distinct(wmd_id, abs_abundance) %>%
      rename(ben_stationary = abs_abundance),
    by = "wmd_id"
  ) %>%
  filter(!near(ben_2020, ben_stationary))

# C1: scaled_abundance sums back to abs_abundance within each WMD x year x scenario
alloc_errors <- data_panel %>%
  group_by(wmd_id, year, rcp, gcm) %>%
  summarise(
    total_abs  = first(abs_abundance),
    sum_scaled = sum(scaled_abundance, na.rm = TRUE),
    diff       = sum_scaled - total_abs,
    .groups    = "drop"
  ) %>%
  filter(!is.na(diff), abs(diff) > 1e-6)

# C2: No NAs in abs_abundance or scaled_abundance in the final panel
n_na_abs_final    <- sum(is.na(data_panel$abs_abundance))
n_na_scaled_final <- sum(is.na(data_panel$scaled_abundance))

# C3: abs_abundance and scaled_abundance are constant across all stationary
#     time steps (2030-2100) and the 2020 baseline within each EAU.
#     These should all match because the stationary scenario holds 2020
#     landcover and abundance fixed throughout the planning horizon.
stationary_vary <- data_panel %>%
  filter(rcp == "stationary" | year == 2020) %>%
  group_by(eau_id) %>%
  summarise(
    n_distinct_abs    = n_distinct(abs_abundance),
    n_distinct_scaled = n_distinct(scaled_abundance),
    .groups = "drop"
  ) %>%
  filter(n_distinct_abs > 1 | n_distinct_scaled > 1)

# C4: Row count per EAU and total rows in final panel
bad_eau_rowcounts <- data_panel %>%
  count(eau_id, name = "n_rows") %>%
  filter(n_rows != rows_per_eau)


# ── Run checks ───────────────────────────────────────────────────────────────

checks <- list(
  "All benefit_new anchor rows present in interpolation input"          = nrow(unmatched) == 0,
  "Anchor values preserved correctly after interpolation"               = nrow(value_mismatches) == 0 & nrow(obs_mismatches) == 0,
  "benefit_decadal: 1 row per WMD at 2020"                             = nrow(bad_2020_ben) == 0,
  "benefit_decadal: 4 rows per WMD at future decades"                  = nrow(bad_future_ben) == 0,
  "No NAs in abs_abundance after panel join"                            = n_na_abs == 0,
  "abs_abundance constant within WMD x year x scenario"                = nrow(abs_inconsistent) == 0,
  "2020 and stationary abs_abundance match for all WMDs"               = nrow(stationary_mismatch) == 0,
  "scaled_abundance sums to abs_abundance within each group"           = nrow(alloc_errors) == 0,
  "No NAs in abs_abundance or scaled_abundance in final panel"         = n_na_abs_final == 0 & n_na_scaled_final == 0,
  "abs/scaled_abundance constant across stationary time steps per EAU" = nrow(stationary_vary) == 0,
  "Each EAU has correct row count in final panel"                      = nrow(bad_eau_rowcounts) == 0,
  "Final panel has expected total row count"                           = nrow(data_panel) == n_rows_expected
)

for (nm in names(checks)) {
  cat(sprintf("  %s  %s\n", if (checks[[nm]]) "[PASS]" else "[FAIL]", nm))
}

failures <- names(checks)[!unlist(checks)]

if (length(failures) > 0) {
  
  cat("\n  --- Diagnostic detail ---\n")
  
  if (!checks[["All benefit_new anchor rows present in interpolation input"]]) {
    cat("  Unmatched rows:", nrow(unmatched), "\n")
    print(unmatched)
  }
  if (!checks[["Anchor values preserved correctly after interpolation"]]) {
    cat("  Future anchor value mismatches:", nrow(value_mismatches), "\n")
    print(value_mismatches %>% select(wmd_id, year, rcp, gcm, ben_new, ben_decadal))
    cat("  Obs anchor value mismatches:", nrow(obs_mismatches), "\n")
    print(obs_mismatches)
  }
  if (!checks[["benefit_decadal: 1 row per WMD at 2020"]]) {
    cat("  WMDs with wrong 2020 row count:\n")
    print(bad_2020_ben)
  }
  if (!checks[["benefit_decadal: 4 rows per WMD at future decades"]]) {
    cat("  WMD-years with wrong future row count:\n")
    print(bad_future_ben)
  }
  if (!checks[["No NAs in abs_abundance after panel join"]]) {
    cat("  NA count in abs_abundance:", n_na_abs, "\n")
    print(eau_panel_with_ben %>%
            filter(is.na(abs_abundance)) %>%
            distinct(wmd_id, year, rcp, gcm))
  }
  if (!checks[["abs_abundance constant within WMD x year x scenario"]]) {
    cat("  Groups with inconsistent abs_abundance:", nrow(abs_inconsistent), "\n")
    print(head(abs_inconsistent, 10))
  }
  if (!checks[["2020 and stationary abs_abundance match for all WMDs"]]) {
    cat("  WMDs with mismatched 2020 vs stationary abundance:", nrow(stationary_mismatch), "\n")
    print(stationary_mismatch)
  }
  if (!checks[["scaled_abundance sums to abs_abundance within each group"]]) {
    cat("  Groups with allocation error > 1e-6:", nrow(alloc_errors), "\n")
    print(head(alloc_errors %>% arrange(desc(abs(diff))), 10))
  }
  if (!checks[["No NAs in abs_abundance or scaled_abundance in final panel"]]) {
    cat("  NAs in abs_abundance (final panel):   ", n_na_abs_final, "\n")
    cat("  NAs in scaled_abundance (final panel):", n_na_scaled_final, "\n")
    print(data_panel %>%
            filter(is.na(abs_abundance) | is.na(scaled_abundance)) %>%
            distinct(wmd_id, year, rcp, gcm))
  }
  if (!checks[["abs/scaled_abundance constant across stationary time steps per EAU"]]) {
    cat("  EAUs with varying stationary/baseline abundance:", nrow(stationary_vary), "\n")
    print(head(stationary_vary, 10))
  }
  if (!checks[["Each EAU has correct row count in final panel"]]) {
    cat("  Expected rows per EAU:", rows_per_eau, "\n")
    cat("  EAUs with wrong row count:", nrow(bad_eau_rowcounts), "\n")
    print(bad_eau_rowcounts %>% count(n_rows))
  }
  if (!checks[["Final panel has expected total row count"]]) {
    cat("  Expected:", n_rows_expected, "| Found:", nrow(data_panel), "\n")
  }
  
  cat("========================================\n\n")
  stop("Logic check FAILED: benefit data has errors. ",
       "Investigate before proceeding to downstream analysis.\n",
       "Failed checks: ", paste(failures, collapse = ", "))
  
} else {
  cat(sprintf(
    "\n  All checks passed. WMDs: %d | EAUs: %d | Total rows: %d\n",
    n_distinct(data_panel$wmd_id),
    n_distinct(data_panel$eau_id),
    nrow(data_panel)
  ))
  cat("========================================\n\n")
}

saveRDS(data_panel, "derived_data/panel_04_benefit.rds")


