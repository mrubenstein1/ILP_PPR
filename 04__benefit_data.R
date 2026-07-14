
##### Benefit (Duck abundance) Data ####
    # Import and clean benefit data
    # add to the established data panel

# ══════════════════════════════════════════════════════════════════════════════
# TEMPORAL ALIGNMENT OF THE BENEFIT ANCHORS  (ALIGNMENT switch, below)
# ──────────────────────────────────────────────────────────────────────────────
# The benefit source (benefit.csv) reports three windows, NOT single years:
#   "Obs"   = average observations   2007-2017   (window midpoint ~2012)
#   "3059"  = average projections    2030-2059   (window midpoint  2045)
#   "7099"  = average projections    2070-2099   (window midpoint  2085)
#
# Two alignments are supported via the ALIGNMENT switch:
#
#   "current"  — the original pipeline: anchor Obs/3059/7099 at 2020/2050/2100 and
#                use the raw Obs value directly as the 2020 baseline. Retained so
#                the validated production numbers are exactly reproducible.
#
#   "midpoint" — anchor each window at its TEMPORAL MIDPOINT: 2012 / 2045 / 2085.
#                The decision clock is UNCHANGED (period 0 is still 2020, first
#                decision enacted over 2020->2030). 2020 is therefore a MODELLED
#                decision baseline, produced by linear interpolation on the
#                2012->2045 segment, then AVERAGED across the four RCP x GCM combos
#                into one shared period-0 row (mirrors the shared-baseline
#                convention script 05 uses for the 2020 hazard). Decades after the
#                last anchor (2090, 2100) are HELD CONSTANT at the 2085 value.
#
# IMPORTANT (design): 2020 is a decision anchor, not an observation year. Under
# "midpoint" it is an estimate produced by projecting the early-2010s observation
# forward to 2020 — consistent with how land cover reaches 2020 via a modelled
# layer and cost reaches 2020 via inflation. Cost (06) and hazard (05) are benefit
# -independent, so switching ALIGNMENT changes ONLY the benefit trajectory; the
# hazard and cost columns are identical across alignments. Any downstream change is
# therefore attributable to the benefit alignment alone.
#
# Full-pipeline workflow to produce a solver-ready panel per alignment (for 08):
#   set ALIGNMENT; run 04 -> 05 -> 06 (05/06 fill trans_prob/cost and save
#   input_data/data_panel.rds); snapshot that file to data_panel_<alignment>.rds;
#   repeat for the other alignment; point 08 at each. This script additionally
#   writes a benefit-stage snapshot (input_data/data_panel_benefit_<alignment>.rds)
#   that the Layer-1 comparison (benefit_alignment_comparison.R) reads directly —
#   no solver required.
# ══════════════════════════════════════════════════════════════════════════════


# ── 0. Setup (safe to run standalone, once per alignment) ─────────────────────####
# Libraries are normally loaded by script 01; re-assert the ones this script uses so
# it can be re-run in a fresh session for each alignment without re-sourcing 01.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(readr); library(tibble)
})

# eau_panel is produced by script 03. Load from disk if it is not already in memory,
# so each alignment run starts from the same clean pre-benefit panel (reproducible
# and independent between runs).
if (!exists("eau_panel")) eau_panel <- readRDS("input_data/eau_panel.rds")


# ══ ALIGNMENT SWITCH ══════════════════════════════════════════════════════════####
# "current"  = 2020 / 2050 / 2100 anchors, raw Obs as 2020 baseline (production).
# "midpoint" = 2012 / 2045 / 2085 anchors, averaged interpolated 2020 baseline,
#              held-flat tail after 2085.
ALIGNMENT <- "midpoint"          # "current" | "midpoint"
# ══════════════════════════════════════════════════════════════════════════════

stopifnot(ALIGNMENT %in% c("current", "midpoint"))

# Anchor years for each window under the chosen alignment.
ANCHOR_YEARS <- if (identical(ALIGNMENT, "current")) {
  c(Obs = 2020L, mid = 2050L, end = 2100L)
} else {
  c(Obs = 2012L, mid = 2045L, end = 2085L)
}

# Decision-model periods (the clock is the SAME for both alignments).
DECISION_YEARS <- as.integer(seq(2020, 2100, by = 10))

cat(sprintf(
  "\n[04] Benefit alignment = '%s'  |  anchors: Obs=%d  mid=%d  end=%d\n",
  ALIGNMENT, ANCHOR_YEARS[["Obs"]], ANCHOR_YEARS[["mid"]], ANCHOR_YEARS[["end"]]))
if (identical(ALIGNMENT, "midpoint"))
  cat("     2020 = interpolated on the Obs->mid segment, averaged across RCPxGCM; ",
      "2090/2100 held flat at the ", ANCHOR_YEARS[["end"]], " value.\n", sep = "")


#### 1. Import data #####
benefit <- read.csv("input_data/benefit.csv")

# restructure around gcm/rcp/time step structure.
# NOTE: `year` is now the ANCHOR YEAR implied by the ALIGNMENT switch (the window
# midpoint under "midpoint"), not a literal observation year.
benefit_new <- benefit %>%
  mutate(
    year = case_when(
      scenario == "Obs"           ~ ANCHOR_YEARS[["Obs"]],
      str_detect(scenario, "3059") ~ ANCHOR_YEARS[["mid"]],  # mid-century window
      str_detect(scenario, "7099") ~ ANCHOR_YEARS[["end"]],  # end-century window
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


#### 2. Interpolate to decadal decision years ####
# One linear interpolation per WMD x RCP x GCM over the three anchors
# (Obs, mid, end), evaluated at the decision years 2020..2100. `rule = 2` holds the
# last anchor value constant for any decision year beyond it (the post-2085 tail
# under "midpoint"; a no-op under "current", whose last anchor is 2100). The Obs
# anchor is SHARED across combos (a single observed value per WMD); the mid/end
# anchors are combo-specific, so the combos diverge over the horizon.

# ── Step 1: shared Obs anchor (one value per WMD; rcp/gcm are NA at Obs) ──────
ben_obs <- benefit_new %>%
  filter(is.na(rcp)) %>%          # Obs rows carry no rcp/gcm
  select(wmd_id, ben)

# ── Step 2: combo-specific future anchors (mid + end) ────────────────────────
future_anchors <- benefit_new %>%
  filter(!is.na(rcp)) %>%
  select(wmd_id, rcp, gcm, year, ben)

# ── Step 3: attach the shared Obs anchor to every combo at the Obs anchor year ─
anchor_obs_per_combo <- future_anchors %>%
  distinct(wmd_id, rcp, gcm) %>%
  left_join(ben_obs, by = "wmd_id") %>%
  mutate(year = ANCHOR_YEARS[["Obs"]])

anchors <- bind_rows(anchor_obs_per_combo, future_anchors) %>%
  arrange(wmd_id, rcp, gcm, year)          # approx() needs x sorted within group

# ── Step 4: interpolate each combo at the decision years (flat tail via rule=2) ─
benefit_decadal_future <- anchors %>%
  group_by(wmd_id, rcp, gcm) %>%
  group_modify(~ tibble(
    year = DECISION_YEARS,
    ben  = approx(x = .x$year, y = .x$ben, xout = DECISION_YEARS, rule = 2)$y
  )) %>%
  ungroup()

# ── Step 5: shared 2020 baseline = mean across the four combos at 2020 ────────
# Under "current" every combo returns the raw Obs value at 2020, so this mean is
# exactly the Obs value (reproduces the original pipeline). Under "midpoint" the
# four combos differ at 2020 (interpolated on the Obs->mid segment), so this mean
# is the averaged modelled baseline. The 2020 row carries rcp/gcm = NA (baseline).
ben_2020_shared <- benefit_decadal_future %>%
  filter(year == 2020L) %>%
  group_by(wmd_id) %>%
  summarise(ben = mean(ben), .groups = "drop") %>%
  mutate(year = 2020L)

# benefit_decadal: shared 2020 baseline (1 row/WMD) + combo-specific 2030..2100.
# Same structure the downstream join expects.
benefit_decadal <- bind_rows(
  ben_2020_shared,
  benefit_decadal_future %>% filter(year > 2020L)
) %>%
  arrange(wmd_id, year, rcp, gcm)


##### 3. Match Benefit data into EAU_Panel ####

#FIRST: remove the 4 WMDs we are removing from the analysis.
  #3 from Montana (Benton Lake, Bowdoin, and Northeast Montana) because the
  #underlying hydrological data that underpinned the duck abundance estimates performed poorly
  #in that region; and a 4th (Windom) because of a mismatch of spatial extent with the FOREsce
  #landcover data.

# Drop excluded WMDs
wmd_exclude <- c("Windom", "Benton Lake", "Bowdoin", "Northeast Montana")

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

# A0: every benefit row was mapped to an anchor year (catches unexpected scenario
#     strings in benefit.csv that match none of Obs / 3059 / 7099).
n_unmapped_year <- sum(is.na(benefit_new$year))

# A1: the interpolation ANCHORS faithfully carry the source values. We check the
#     INPUT to approx() rather than benefit_decadal, because under "midpoint" the
#     anchor years (2012/2045/2085) are not decision years and never appear in
#     benefit_decadal directly (unlike "current", where they coincide). Exactness
#     of the interpolation AT the knots is guaranteed by approx(); the flat tail is
#     checked separately (A4).
bn_future <- benefit_new %>% filter(!is.na(rcp)) %>% select(wmd_id, rcp, gcm, year, ben)
future_unmatched <- anti_join(bn_future, anchors,
                              by = c("wmd_id", "rcp", "gcm", "year"))
future_valcheck <- inner_join(bn_future, anchors,
                              by = c("wmd_id", "rcp", "gcm", "year"),
                              suffix = c("_src", "_anc")) %>%
  filter(!near(ben_src, ben_anc))

bn_obs <- benefit_new %>% filter(is.na(rcp)) %>% distinct(wmd_id, ben)
obs_valcheck <- bn_obs %>%
  left_join(ben_obs %>% rename(ben_anc = ben), by = "wmd_id") %>%
  filter(is.na(ben_anc) | !near(ben, ben_anc))

# A2: benefit_decadal row structure
#     2020: 1 row per WMD (shared baseline; rcp/gcm = NA)
#     Future decades: 4 rows per WMD x year (rcp45/85 x wet/dry)
bad_2020_ben <- benefit_decadal %>%
  filter(year == 2020) %>%
  count(wmd_id) %>%
  filter(n != 1)

bad_future_ben <- benefit_decadal %>%
  filter(year > 2020) %>%
  count(wmd_id, year) %>%
  filter(n != n_rcp_gcm_combos)

# A4: held-flat tail. For any decision year beyond the last (end) anchor, the value
#     must equal that combo's end-anchor value. Vacuous under "current" (no decision
#     year exceeds the 2100 anchor); active under "midpoint" (2090, 2100 == 2085).
end_year   <- ANCHOR_YEARS[["end"]]
tail_years <- DECISION_YEARS[DECISION_YEARS > end_year]
if (length(tail_years) > 0) {
  end_anchor_vals <- anchors %>%
    filter(year == end_year) %>%
    select(wmd_id, rcp, gcm, ben_end = ben)
  tail_check <- benefit_decadal_future %>%
    filter(year %in% tail_years) %>%
    left_join(end_anchor_vals, by = c("wmd_id", "rcp", "gcm")) %>%
    filter(is.na(ben_end) | !near(ben, ben_end))
} else {
  tail_check <- benefit_decadal_future[0, ]
}

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

# B3: 2020 (baseline) abundance matches stationary abundance for all WMDs
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
  "All benefit rows mapped to an anchor year"                          = n_unmapped_year == 0,
  "All future (mid/end) anchors carried into interpolation input"      = nrow(future_unmatched) == 0,
  "Future anchor values preserved in interpolation input"             = nrow(future_valcheck) == 0,
  "Observation anchor values preserved"                               = nrow(obs_valcheck) == 0,
  "benefit_decadal: 1 row per WMD at 2020"                            = nrow(bad_2020_ben) == 0,
  "benefit_decadal: 4 rows per WMD at future decades"                 = nrow(bad_future_ben) == 0,
  "Held-flat tail beyond last anchor == end-anchor value"             = nrow(tail_check) == 0,
  "No NAs in abs_abundance after panel join"                           = n_na_abs == 0,
  "abs_abundance constant within WMD x year x scenario"               = nrow(abs_inconsistent) == 0,
  "2020 and stationary abs_abundance match for all WMDs"              = nrow(stationary_mismatch) == 0,
  "scaled_abundance sums to abs_abundance within each group"          = nrow(alloc_errors) == 0,
  "No NAs in abs_abundance or scaled_abundance in final panel"        = n_na_abs_final == 0 & n_na_scaled_final == 0,
  "abs/scaled_abundance constant across stationary time steps per EAU"= nrow(stationary_vary) == 0,
  "Each EAU has correct row count in final panel"                     = nrow(bad_eau_rowcounts) == 0,
  "Final panel has expected total row count"                          = nrow(data_panel) == n_rows_expected
)

for (nm in names(checks)) {
  cat(sprintf("  %s  %s\n", if (checks[[nm]]) "[PASS]" else "[FAIL]", nm))
}

failures <- names(checks)[!unlist(checks)]

if (length(failures) > 0) {

  cat("\n  --- Diagnostic detail ---\n")

  if (!checks[["All benefit rows mapped to an anchor year"]]) {
    cat("  benefit.csv rows with unmapped scenario string:", n_unmapped_year, "\n")
    print(benefit %>% mutate(.mapped = benefit_new$year) %>%
            filter(is.na(.mapped)) %>% distinct(scenario) %>% head(10))
  }
  if (!checks[["All future (mid/end) anchors carried into interpolation input"]]) {
    cat("  Future anchor rows missing from interpolation input:", nrow(future_unmatched), "\n")
    print(head(future_unmatched, 10))
  }
  if (!checks[["Future anchor values preserved in interpolation input"]]) {
    cat("  Future anchor value mismatches:", nrow(future_valcheck), "\n")
    print(head(future_valcheck %>% select(wmd_id, rcp, gcm, year, ben_src, ben_anc), 10))
  }
  if (!checks[["Observation anchor values preserved"]]) {
    cat("  Obs anchor value mismatches:", nrow(obs_valcheck), "\n")
    print(head(obs_valcheck, 10))
  }
  if (!checks[["benefit_decadal: 1 row per WMD at 2020"]]) {
    cat("  WMDs with wrong 2020 row count:\n"); print(bad_2020_ben)
  }
  if (!checks[["benefit_decadal: 4 rows per WMD at future decades"]]) {
    cat("  WMD-years with wrong future row count:\n"); print(bad_future_ben)
  }
  if (!checks[["Held-flat tail beyond last anchor == end-anchor value"]]) {
    cat("  Tail rows not equal to end-anchor value:", nrow(tail_check), "\n")
    print(head(tail_check, 10))
  }
  if (!checks[["No NAs in abs_abundance after panel join"]]) {
    cat("  NA count in abs_abundance:", n_na_abs, "\n")
    print(eau_panel_with_ben %>% filter(is.na(abs_abundance)) %>%
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
    "\n  All checks passed. Alignment: %s | WMDs: %d | EAUs: %d | Total rows: %d\n",
    ALIGNMENT,
    n_distinct(data_panel$wmd_id),
    n_distinct(data_panel$eau_id),
    nrow(data_panel)
  ))
  cat("========================================\n\n")
}


### 7. Save alignment-tagged benefit-stage panel ####
# This is the benefit-stage snapshot (trans_prob / cost still NA — those are filled
# by 05 / 06). The Layer-1 comparison reads these tagged files directly; no solver
# needed. The canonical input_data/data_panel.rds is intentionally NOT written here
# (05 owns the first full save), so the sequential 04->05->06 run is unchanged.
dir.create("input_data", showWarnings = FALSE, recursive = TRUE)
saveRDS(data_panel, sprintf("input_data/data_panel_benefit_%s.rds", ALIGNMENT))
write_csv(data_panel, sprintf("input_data/data_panel_benefit_%s.csv", ALIGNMENT))
cat(sprintf("  Wrote input_data/data_panel_benefit_%s.{rds,csv}\n\n", ALIGNMENT))
