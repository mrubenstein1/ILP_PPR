### Transition Probability (Risk of Loss) ####
#
# Operationalization:
#   trans_prob[i, r, t] = max(epsilon, prop_lost_unprot[i, r, t] / prop_suit_unprot[i, r, t])
#
# This is the GROSS conversion rate: of the UNPROTECTED habitat standing at t,
# what fraction converts to unsuitable by t+1. Protected land (PAD-US GAP 1-3,
# script 02 section 4b) cannot convert and is excluded from numerator and denominator.
#
# ── WHY GROSS AND NOT NET ─────────────────────────────────────────────────────
# This script previously differenced prop_suitable:
#
#   (ps_t - ps_t1) / ps_t
#
# Writing L for area converting suitable -> unsuitable and G for the reverse,
# ps_t1 = ps_t - L + G, so that expression is (L - G) / ps_t -- the NET rate,
# where L / ps_t was intended. The two agree only where G = 0.
#
# Diagnostic testing (see archive) measured L and G separately at native FOREsce resolution
# and found, across the 879 modelled EAUs:
#   * the net metric conceals 20.4% of all habitat conversion
#   * it understates the landscape hazard by 25.6%
#   * 8.0% of EAU-years lose habitat but show net <= 0, so they are floored to
#     epsilon and enter the model as effectively risk-free
#   * in four mid-century decade-RCP cells the mean net rate is NEGATIVE while
#     real conversion continues -- netting inverts the sign of the hazard
#
# L cannot be recovered from prop_suitable, because aggregating to EAU scale nets
# the two flows irreversibly. It is therefore measured pixel-by-pixel in script
# 02 section 5, translated to eau_id in script 03 section 10b, and read here.
#
# Consequences for this script: trans_raw is now non-negative by construction, so
# the "recovering" category disappears, epsilon binds only where loss is exactly
# zero (about 18.7% of EAU-years rather than 22.5%), and no forward difference of
# prop_suitable is taken anywhere.
#
# epsilon: minimum background risk of conversion applied to all EAUs with
#          prop_suitable > 0. Reflects that even stable parcels face some
#          real-world conversion risk if unprotected. Set to 0 to disable.
#
# Special cases:
#   prop_suitable[t] = 0  : trans_prob = 0 (no habitat, no conservation urgency)
#   prop_suit_unprot[t] = 0 < prop_suitable[t] : trans_prob = epsilon (all habitat protected)
#   t = 2100 (terminal)   : trans_prob = 0 (no subsequent decision period)
#   Stationary scenario   : trans_prob = epsilon (background risk applies)
#   2020 baseline row     : mean of RCP 4.5 and 8.5 2020→2030 transitions
#
# INPUTS:  data_panel (from script 04)
#          derived_data/eau_flows.rds (from script 03 section 10b)
# OUTPUTS: data_panel with trans_prob column populated; saved to derived_data/

#------ Setup ________
if (!exists(".SETUP_DONE")) source("00_setup.R")

# ── 0. Load data  ──────────────────────────────────────####
data_panel <- readRDS(file.path(DIR_DERIVED, "panel_04_benefit.rds"))

# Gross habitat flows, keyed by eau_id x year x rcp. `year` is the period the
# transition STARTS from, so these run 2020..2090; 2100 is terminal and has no
# row. Script 03 has already verified that this covers exactly the panel's EAU
# set for every non-terminal period.
flows_path <- file.path(DIR_DERIVED, "eau_flows.rds")
if (!file.exists(flows_path))
  stop("Required input not found: ", flows_path,
       "\n  Run script 02 (section 5, gross flows) and then script 03 first. ",
       "This script can no longer derive the hazard from prop_suitable alone -- ",
       "see WHY GROSS AND NOT NET at the top of this file.")
eau_flows <- readRDS(flows_path)

# eau_flows is built in script 03 against the FULL covered-EAU set (1,169: every
# EAU with land-cover data). This panel has already had the WMD exclusion applied
# and carries fewer. Left unrestricted, the surplus EAUs join to a missing
# prop_suitable, produce NA loss rates, and those NAs then propagate silently
# through the epsilon derivation into every trans_prob in the panel.
#
# Restrict here, once, and report it -- rather than relying on each downstream
# join to drop them by accident.
panel_eaus    <- sort(unique(data_panel$eau_id))
n_flow_before <- n_distinct(eau_flows$eau_id)
eau_flows     <- eau_flows %>% filter(eau_id %in% panel_eaus)

cat(sprintf("\n  Gross flows: %d EAUs available, %d in this panel (%d dropped by the\n",
            n_flow_before, n_distinct(eau_flows$eau_id),
            n_flow_before - n_distinct(eau_flows$eau_id)))
cat("  WMD exclusion upstream of this script).\n")

# The reverse direction is the one that would be fatal: a panel EAU with no flow
# data becomes an NA hazard. Check it explicitly rather than discovering it as an
# NA count in section 9.
missing_flow_eaus <- setdiff(panel_eaus, unique(eau_flows$eau_id))
if (length(missing_flow_eaus))
  stop(length(missing_flow_eaus), " EAU(s) in the panel have no gross flow data ",
       "(e.g. ", paste(head(missing_flow_eaus, 5), collapse = ", "), "). ",
       "eau_flows does not cover the panel -- re-run script 02 section 5 and ",
       "script 03 section 10b.")


# ══ PARAMETER ════════════════════════════════════════════════════════════════####
# epsilon is derived from the data in Section 6: it is set to one order of
# magnitude below the smallest observed positive loss rate, so that all real
# loss signal is preserved while stable/recovering EAUs get a non-zero floor.
#
# To override with a fixed value for sensitivity analyses, set this to a number:
epsilon_override <- NULL   # e.g., epsilon_override <- 0.005
# ══════════════════════════════════════════════════════════════════════════════


# ── 1. Calculate raw proportional loss (GROSS) ───────────####
#Note: this is computed as a habitat erosion rate, but will be
  # interpreted as a per-period Bernoulli conversion probability in the final ILP.

  # From t to t+1, what is the expected % of an EAU's suitable habitat that is
  # expected to convert to unsuitable habitat?
  # in other words, what proportion of currently suitable habitat will be lost?
  # how much urgency is there to acquire this parcel in this time period?
  # parcels with high proportional loss are those expected to experience habitat loss and
  # should be urgently considered for acquisition; those that have little expected loss may be less urgent.

  # prop_lost is the area converting suitable -> unsuitable, measured directly
  # rather than inferred from the change in stock. Habitat converting the other
  # way (prop_gained) is NOT subtracted: restoration elsewhere in an EAU does not
  # reduce the risk facing habitat currently standing there, and acquisition does
  # not cause the gain. See WHY GROSS AND NOT NET at the top of this file.

  # prop_lost pixels are by construction a subset of the pixels suitable at t, so
  # the ratio is bounded in [0, 1] and can never be negative. That is a change
  # from the net formulation, which could go negative wherever gain outran loss.

gross_proportional_loss <- function(ps_start, prop_lost) {
  case_when(
    ps_start == 0        ~ 0,                     # no habitat to lose
    TRUE                 ~ prop_lost / ps_start   # gross rate, in [0, 1]
  )
}


# ── 2. Extract unique prop_suitable time series per EAU × RCP ─────────────────####

prop_ts <- data_panel %>%
  filter(rcp %in% c("baseline", "45", "85")) %>%
  distinct(eau_id, year, rcp, prop_suitable) %>%
  arrange(eau_id, rcp, year)


# ── 3. Compute raw 2020 → 2030 transition (cross-RCP boundary step) ───────────####
  #Note: At 2020 decision point, decision maker must estimate expected loss by 2030.
  # because no 2030 RCP has materialized, this calculation takes the average of the
  # 4.5 and 8.5 gross loss rates (pooled in section 5 below).
  #
  # ps_start is the meanclim 2020 stock, and the 2020 flows in eau_flows were
  # measured against that same meanclim raster, so numerator and denominator
  # refer to the same landscape.

ps_2020 <- prop_ts %>%
  filter(year == 2020, rcp == "baseline") %>%
  select(eau_id, ps_start = prop_suitable)

trans_2020_raw <- eau_flows %>%
  filter(year == 2020L, rcp %in% c("45", "85")) %>%
  select(eau_id, year, rcp, prop_lost_unprot, prop_suit_unprot) %>%
  mutate(trans_raw = gross_proportional_loss(prop_suit_unprot, prop_lost_unprot)) %>%
  select(eau_id, year, rcp, trans_raw)


# ── 4. Compute raw 2030 → 2100 transitions (within-RCP forward differences) ─── ####
  # No forward difference of prop_suitable is taken any more: the loss for the
  # period starting at `year` is read directly from eau_flows. 2100 has no flow
  # row (it is terminal), so the left_join leaves prop_lost NA there and the
  # terminal branch sets the hazard to 0.

trans_future_raw <- prop_ts %>%
  filter(rcp %in% c("45", "85"), year >= 2030) %>%
  left_join(eau_flows %>% select(eau_id, year, rcp, prop_lost_unprot, prop_suit_unprot),
            by = c("eau_id", "year", "rcp")) %>%
  mutate(
    trans_raw = case_when(
      year == 2100 ~ 0,    # terminal period: no next decision
      TRUE         ~ gross_proportional_loss(prop_suit_unprot, prop_lost_unprot)
    )
  ) %>%
  select(eau_id, year, rcp, trans_raw)

# A missing flow on a NON-terminal row would silently become an NA hazard and
# then an NA objective coefficient in the ILP. Fail here instead, where the
# cause is still visible.
n_missing_flow <- sum(is.na(trans_future_raw$trans_raw))
if (n_missing_flow > 0)
  stop(n_missing_flow, " non-terminal EAU-year(s) have no gross flow data. ",
       "eau_flows does not cover the panel. Re-run script 02 section 5 and ",
       "script 03 section 10b.")


# ── 5. Pool raw estimates for the 2020 baseline row ──────────####

trans_baseline_2020_raw <- trans_2020_raw %>%
  group_by(eau_id) %>%
  summarise(trans_raw = mean(trans_raw), .groups = "drop") %>%
  mutate(year = 2020L, rcp = "baseline")


# ── 6. DIAGNOSTIC: examine raw distribution and derive epsilon ──────####

raw_all <- bind_rows(trans_2020_raw, trans_future_raw) %>%
  filter(year < 2100)    # exclude terminal year (always 0 by definition)

# ── Derive epsilon from the data ──────────────────────────────────────────────
# epsilon = one order of magnitude below the smallest observed positive loss,
# so every real loss signal retains its actual value and is distinguishable
# from the stable/recovering floor.
min_positive_loss <- min(raw_all$trans_raw[raw_all$trans_raw > 0])
epsilon <- if (!is.null(epsilon_override)) epsilon_override else min_positive_loss / 10

# A non-finite epsilon is silent and catastrophic: apply_floor() computes
# pmax(NA, trans_raw), so EVERY non-terminal trans_prob in the panel becomes NA
# and the failure only surfaces as an NA count in section 9, far from its cause.
# The usual reason is an NA in trans_raw, which means some EAU-year has loss data
# but no matching prop_suitable (or vice versa). Name it here instead.
if (!is.finite(epsilon)) {
  n_raw_na <- sum(is.na(raw_all$trans_raw))
  stop("epsilon is not finite (", epsilon, "). ",
       n_raw_na, " of ", nrow(raw_all), " raw loss rates are NA.\n",
       "  This means prop_lost and prop_suitable do not cover the same EAU-years. ",
       "Check that eau_flows was restricted to the panel EAU set in section 0.")
}

#6a. Describe transition prob distribution (raw) ####
cat("\n=======================================================\n")
cat("  RAW trans_prob DISTRIBUTION (before floor applied)\n")
cat("=======================================================\n\n")

# Table 1: counts across all EAU-years, broken into the three trajectory types.
# NOTE: under gross accounting "Recovering (n)" must be 0 in every row -- loss is
# a count of converting pixels and cannot be negative. The column is retained
# precisely because a non-zero entry would be a visible signal that the flows are
# wrong or mis-joined. Under the old net metric it was routinely large.
cat("-- Table 1: All EAU-years by trajectory type --\n\n")

summary_counts <- raw_all %>%
  group_by(rcp) %>%
  summarise(
    n_total      = n(),
    n_recovering = sum(trans_raw < 0),
    pct_recov    = round(100 * n_recovering / n_total, 1),
    n_stable     = sum(trans_raw == 0),
    pct_stable   = round(100 * n_stable / n_total, 1),
    n_losing     = sum(trans_raw > 0),
    pct_losing   = round(100 * n_losing / n_total, 1),
    .groups      = "drop"
  ) %>%
  rename(
    "RCP"                          = rcp,
    "Total EAU-years"              = n_total,
    "Recovering (n)"               = n_recovering,
    "Recovering (%)"               = pct_recov,
    "Stable (n)"                   = n_stable,
    "Stable (%)"                   = pct_stable,
    "Losing habitat (n)"           = n_losing,
    "Losing habitat (%)"           = pct_losing
  )

print(summary_counts, width = Inf)

# Table 2: distributional stats for losing EAU-years only
cat("\n-- Table 2: Loss rate distribution (losing EAU-years only) --\n\n")

summary_loss <- raw_all %>%
  filter(trans_raw > 0) %>%
  group_by(rcp) %>%
  summarise(
    n            = n(),
    min_loss     = round(min(trans_raw),              8),
    p25          = round(quantile(trans_raw, 0.25),   6),
    median_loss  = round(median(trans_raw),           6),
    mean_loss    = round(mean(trans_raw),             6),
    p75          = round(quantile(trans_raw, 0.75),   6),
    p90          = round(quantile(trans_raw, 0.90),   6),
    max_loss     = round(max(trans_raw),              6),
    .groups      = "drop"
  ) %>%
  rename(
    "RCP"                          = rcp,
    "n (losing EAU-years)"         = n,
    "Min loss rate"                = min_loss,
    "25th percentile"              = p25,
    "Median loss rate"             = median_loss,
    "Mean loss rate"               = mean_loss,
    "75th percentile"              = p75,
    "90th percentile"              = p90,
    "Max loss rate (worst parcel)" = max_loss
  )

print(summary_loss, width = Inf)

cat(sprintf("\n  Smallest positive loss rate:   %.2e\n", min_positive_loss))
cat(sprintf("  epsilon (floor = above / 10):  %.2e\n", epsilon))
if (!is.null(epsilon_override)) {
  cat(sprintf("  *** epsilon_override active: using %.2e instead of data-derived value ***\n",
              epsilon_override))
}


# ── 7. Apply epsilon floor and assemble final trans_prob table ────────────────

apply_floor <- function(trans_raw, ps, is_terminal, eps) {
  case_when(
    is_terminal  ~ 0,          # terminal period: no next decision
    ps == 0      ~ 0,          # no habitat: no conservation urgency
    TRUE         ~ pmax(eps, trans_raw)   # floor at epsilon
  )
}

# RCP 4.5 and 8.5 (2020 step)
trans_2020 <- trans_2020_raw %>%
  left_join(ps_2020 %>% rename(ps_start_val = ps_start), by = "eau_id") %>%
  mutate(trans_prob = apply_floor(trans_raw, ps_start_val,
                                  is_terminal = FALSE, eps = epsilon)) %>%
  select(eau_id, year, rcp, trans_prob)

# RCP 4.5 and 8.5 (2030–2100)
trans_future <- trans_future_raw %>%
  left_join(prop_ts %>% select(eau_id, year, rcp, prop_suitable),
            by = c("eau_id", "year", "rcp")) %>%
  mutate(trans_prob = apply_floor(trans_raw, prop_suitable,
                                  is_terminal = (year == 2100), eps = epsilon)) %>%
  select(eau_id, year, rcp, trans_prob)

# Stationary: background risk applies (except terminal year)
trans_stationary <- data_panel %>%
  filter(rcp == "stationary") %>%
  distinct(eau_id, year, rcp) %>%
  mutate(trans_prob = if_else(year == 2100, 0, epsilon))

# 2020 baseline: mean of floored RCP 4.5 and 8.5 values
trans_baseline_2020 <- trans_2020 %>%
  group_by(eau_id) %>%
  summarise(trans_prob = mean(trans_prob), .groups = "drop") %>%
  mutate(year = 2020L, rcp = "baseline")

# Combine
trans_prob_table <- bind_rows(
  trans_2020,
  trans_future,
  trans_stationary,
  trans_baseline_2020
) %>%
  arrange(eau_id, rcp, year)


# ── 8. Join back into data_panel ──────────────────────────────────────────────

data_panel <- data_panel %>%
  select(-trans_prob) %>%
  left_join(trans_prob_table, by = c("eau_id", "year", "rcp"))


# ── 9. Logic checks ───────────────────────────────────────────────────────────

cat("\n========================================\n")
cat("  LOGIC CHECK: transition probability\n")
cat("========================================\n")

n_na       <- sum(is.na(data_panel$trans_prob))
n_oob      <- sum(data_panel$trans_prob < 0 | data_panel$trans_prob > 1, na.rm = TRUE)

terminal_nonzero <- data_panel %>%
  filter(year == 2100, trans_prob != 0)

below_floor <- data_panel %>%
  filter(prop_suitable > 0, year < 2100, trans_prob < epsilon)

gcm_mismatch <- data_panel %>%
  filter(rcp %in% c("45", "85")) %>%
  group_by(eau_id, year, rcp) %>%
  summarise(n_unique = n_distinct(trans_prob), .groups = "drop") %>%
  filter(n_unique > 1)

stationary_wrong <- data_panel %>%
  filter(rcp == "stationary", year < 2100, trans_prob != epsilon)

# CHECK: gross loss cannot be negative. prop_lost counts pixels converting
# suitable -> unsuitable, so the rate is bounded in [0, 1] by construction. A
# negative value means the flows are corrupt or joined to the wrong EAU. Under
# the old net metric negative values were routine, which is what epsilon existed
# to absorb; under gross they are impossible, so this is a real check.
negative_raw <- raw_all %>% filter(trans_raw < 0)

# CHECK: the stock identity, ps_t1 == ps_t - L + G, through the full join chain.
# Script 02 section 6 already verifies this at cell level. What it cannot see is
# whether the values survived cell_id -> eau_id (script 03) and the join into the
# panel intact. That is what this tests: a mis-join would leave the residual
# non-zero even though the raster arithmetic was correct.
#
# lead() appears here only as a check. The estimator above no longer differences
# prop_suitable anywhere.
IDENTITY_TOL <- 1e-6

identity_resid <- prop_ts %>%
  filter(rcp %in% c("45", "85"), year >= 2030) %>%
  arrange(eau_id, rcp, year) %>%
  group_by(eau_id, rcp) %>%
  mutate(ps_next = lead(prop_suitable)) %>%
  ungroup() %>%
  filter(!is.na(ps_next)) %>%          # drops 2100, which has no successor
  left_join(eau_flows %>% select(eau_id, year, rcp, prop_lost, prop_gained,
                                 prop_suit_unprot),
            by = c("eau_id", "year", "rcp")) %>%
  mutate(resid = ps_next - (prop_suitable - prop_lost + prop_gained))

max_identity_resid <- suppressWarnings(max(abs(identity_resid$resid), na.rm = TRUE))
n_identity_na      <- sum(is.na(identity_resid$resid))

checks <- list(
  "No NAs in trans_prob"                               = n_na == 0,
  "All trans_prob values in [0, 1]"                    = n_oob == 0,
  "Terminal year (2100) trans_prob always 0"           = nrow(terminal_nonzero) == 0,
  "No non-terminal EAU (prop > 0) below epsilon floor" = nrow(below_floor) == 0,
  "trans_prob identical for wet/dry within RCP"        = nrow(gcm_mismatch) == 0,
  "Stationary non-terminal trans_prob = epsilon"       = nrow(stationary_wrong) == 0,
  "Gross loss is never negative"                       = nrow(negative_raw) == 0,
  "Stock identity survives the join chain"             =
    is.finite(max_identity_resid) && max_identity_resid < IDENTITY_TOL &&
    n_identity_na == 0,
  "Unprotected habitat <= prop_suitable (join chain)" =
    !any(identity_resid$prop_suit_unprot > identity_resid$prop_suitable + 1e-3)
)

for (nm in names(checks)) {
  cat(sprintf("  %s  %s\n", if (checks[[nm]]) "[PASS]" else "[FAIL]", nm))
}

failures <- names(checks)[!unlist(checks)]

if (length(failures) > 0) {
  
  cat("\n  --- Diagnostic detail ---\n")
  if (!checks[["No NAs in trans_prob"]])
    cat("  NA count:", n_na, "\n")
  if (!checks[["All trans_prob values in [0, 1]"]])
    cat("  Out-of-range count:", n_oob, "\n")
  if (!checks[["Terminal year (2100) trans_prob always 0"]])
    cat("  Non-zero 2100 rows:", nrow(terminal_nonzero), "\n")
  if (!checks[["No non-terminal EAU (prop > 0) below epsilon floor"]])
    cat("  Below-floor rows:", nrow(below_floor), "\n")
  if (!checks[["trans_prob identical for wet/dry within RCP"]])
    cat("  GCM-mismatched groups:", nrow(gcm_mismatch), "\n")
  if (!checks[["Stationary non-terminal trans_prob = epsilon"]])
    cat("  Wrong stationary rows:", nrow(stationary_wrong), "\n")
  if (!checks[["Gross loss is never negative"]]) {
    cat("  Negative gross loss rows:", nrow(negative_raw), "\n")
    cat("  This should be impossible. The flows are corrupt or mis-joined --\n")
    cat("  re-run script 02 section 5 and check its logic check output.\n")
    print(head(as.data.frame(negative_raw), 10), row.names = FALSE)
  }
  if (!checks[["Stock identity survives the join chain"]]) {
    cat(sprintf("  Max identity residual: %.3e (tolerance %.0e), NA rows: %d\n",
                max_identity_resid, IDENTITY_TOL, n_identity_na))
    cat("  ps_t1 != ps_t - L + G. Script 02's cell-level check passed, so this\n")
    cat("  points at the cell_id -> eau_id join in script 03, not the rasters.\n")
  }
  
  cat("========================================\n\n")
  stop("Logic check FAILED. Investigate before proceeding.\n",
       "Failed checks: ", paste(failures, collapse = ", "))
  
} else {
  cat(sprintf(
    "\n  All checks passed. epsilon = %.2e | Rows: %d\n",
    epsilon, nrow(data_panel)
  ))
  cat("========================================\n\n")
}


# ── 10. Post-floor summary ────────────────────────────────────────────────────

cat("--- Final trans_prob summary by RCP × year (excluding terminal 2100) ---\n\n")

data_panel %>%
  filter(rcp %in% c("45", "85"), year < 2100) %>%
  distinct(eau_id, year, rcp, trans_prob) %>%
  group_by(rcp, year) %>%
  summarise(
    mean_tp   = round(mean(trans_prob),   4),
    median_tp = round(median(trans_prob), 4),
    min_tp    = round(min(trans_prob),    4),
    max_tp    = round(max(trans_prob),    4),
    .groups   = "drop"
  ) %>%
  print(n = Inf)


# ── 11. Save ──────────────────────────────────────────────────────────────────

saveRDS(data_panel,  file.path(DIR_DERIVED, "panel_05_risk.rds"))
write_csv(data_panel, file.path(DIR_DERIVED, "panel_05_risk.csv"))

cat(sprintf("\n✓ data_panel saved with trans_prob (epsilon = %.2e).\n", epsilon))
