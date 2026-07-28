### Transition Probability (Risk of Loss) ####
#
# Operationalization:
#   trans_prob[i, r, t] = max(epsilon, (prop_suitable_t - prop_suitable_t+1) / prop_suitable_t)
#
# epsilon: minimum background risk of conversion applied to all EAUs with
#          prop_suitable > 0. Reflects that even stable or recovering parcels
#          face some real-world conversion risk if unprotected. Set to 0 to
#          disable. 
#
# Special cases:
#   prop_suitable[t] = 0  : trans_prob = 0 (no habitat, no conservation urgency)
#   t = 2100 (terminal)   : trans_prob = 0 (no subsequent decision period)
#   Stationary scenario   : trans_prob = epsilon (background risk applies)
#   2020 baseline row     : mean of RCP 4.5 and 8.5 2020→2030 transitions
#
# INPUTS:  data_panel (from script 04)
# OUTPUTS: data_panel with trans_prob column populated; saved to derived_data/

#------ Setup ________
if (!isTRUE(.SETUP_DONE)) source("00_setup.R")

# ── 0. Load data  ──────────────────────────────────────####
data_panel <- readRDS("derived_data/panel_04_benefit.rds")


# ══ PARAMETER ════════════════════════════════════════════════════════════════####
# epsilon is derived from the data in Section 6: it is set to one order of
# magnitude below the smallest observed positive loss rate, so that all real
# loss signal is preserved while stable/recovering EAUs get a non-zero floor.
#
# To override with a fixed value for sensitivity analyses, set this to a number:
epsilon_override <- NULL   # e.g., epsilon_override <- 0.005
# ══════════════════════════════════════════════════════════════════════════════


# ── 1. Calculate raw proportional loss ───────────────────####
#Note: this is computed as a habitat erosion rate, but will be 
  # interpreted as a per-period Bernoulli conversion probability in the final ILP. 

  # From t to t+1, what is the expected % of an EAU's suitable habitat that is 
  # expected to convert to unsuitable habitat?
  # in other words, what proportion of currently suitable habitat will be lost?
  # how much urgency is there to acquire this parcel in this time period?
  # parcels with high proportional loss are those expected to experience habitat loss and
  # should be urgently considered for acquisition; those that have little expected loss may be less urgent. 


raw_proportional_loss <- function(ps_start, ps_end) {
  case_when(
    ps_start == 0        ~ 0,                                 # no habitat to lose
    TRUE                 ~ (ps_start - ps_end) / ps_start    # can be negative (recovery)
  )
}


# ── 2. Extract unique prop_suitable time series per EAU × RCP ─────────────────####

prop_ts <- data_panel %>%
  filter(rcp %in% c("baseline", "45", "85")) %>%
  distinct(eau_id, year, rcp, prop_suitable) %>%
  arrange(eau_id, rcp, year)


# ── 3. Compute raw 2020 → 2030 transition (cross-RCP boundary step) ───────────####
  #Note: At 2020 decision point, decision maker must estimate expected loss by 2030.
  # because no 2030 RCP has materialized, this calculations takes the average of 4.5 and 8.5 net loss 

ps_2020 <- prop_ts %>%
  filter(year == 2020, rcp == "baseline") %>%
  select(eau_id, ps_start = prop_suitable)

trans_2020_raw <- prop_ts %>%
  filter(year == 2030, rcp %in% c("45", "85")) %>%
  select(eau_id, rcp, ps_end = prop_suitable) %>%
  left_join(ps_2020, by = "eau_id") %>%
  mutate(
    year      = 2020L,
    trans_raw = raw_proportional_loss(ps_start, ps_end)
  ) %>%
  select(eau_id, year, rcp, trans_raw)


# ── 4. Compute raw 2030 → 2100 transitions (within-RCP forward differences) ─── ####

trans_future_raw <- prop_ts %>%
  filter(rcp %in% c("45", "85"), year >= 2030) %>%
  arrange(eau_id, rcp, year) %>%
  group_by(eau_id, rcp) %>%
  mutate(
    ps_next   = lead(prop_suitable),
    trans_raw = case_when(
      is.na(ps_next) ~ 0,    # terminal period (2100): no next decision
      TRUE           ~ raw_proportional_loss(prop_suitable, ps_next)
    )
  ) %>%
  ungroup() %>%
  select(eau_id, year, rcp, trans_raw)


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

#6a. Describe transition prob distribution (raw) ####
cat("\n=======================================================\n")
cat("  RAW trans_prob DISTRIBUTION (before floor applied)\n")
cat("=======================================================\n\n")

# Table 1: counts across all EAU-years, broken into the three trajectory types
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

checks <- list(
  "No NAs in trans_prob"                               = n_na == 0,
  "All trans_prob values in [0, 1]"                    = n_oob == 0,
  "Terminal year (2100) trans_prob always 0"           = nrow(terminal_nonzero) == 0,
  "No non-terminal EAU (prop > 0) below epsilon floor" = nrow(below_floor) == 0,
  "trans_prob identical for wet/dry within RCP"        = nrow(gcm_mismatch) == 0,
  "Stationary non-terminal trans_prob = epsilon"       = nrow(stationary_wrong) == 0
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

saveRDS(data_panel,  "derived_data/panel_05_risk.rds")
write_csv(data_panel, "derived_data/panel_05_risk.csv")

cat(sprintf("\n✓ data_panel saved with trans_prob (epsilon = %.2e).\n", epsilon))
