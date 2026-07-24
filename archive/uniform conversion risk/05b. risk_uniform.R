# ══════════════════════════════════════════════════════════════════════════════
# 05b__risk_uniform.R  —  Uniform literature-derived conversion risk (EXPERIMENT)
# ══════════════════════════════════════════════════════════════════════════════
#
# Drop-in replacement for the trans_prob column produced by 05__risk_of_loss.R.
# Swaps the data-derived, EAU- and scenario-varying habitat erosion rate for a
# single UNIFORM per-decade conversion probability derived from published annual
# loss rates for the PPR.
#
# ── PURPOSE ───────────────────────────────────────────────────────────────────
# This is an experiment, not the primary analysis. Script 05 derives trans_prob
# from decade-over-decade change in FOREsce prop_suitable, floors it at a
# data-derived epsilon, and discards recovery. Those hazards are small and noisy,
# and they carry BOTH of the signals the rolling policy exploits: the future
# benefit path AND the future risk path. Making risk uniform shuts off the second
# channel — so whatever rolling-vs-myopic gap survives is attributable to BENEFIT
# foresight alone. Running the sweep below decomposes the effect and maps its
# sensitivity to risk magnitude across a ~14x span.
#
# Three structural consequences of uniformity, for interpreting the output:
#   (1) Survival collapses to a function of time alone: S[i,t] = S[t] for all i.
#       Risk no longer DIFFERENTIATES parcels; it becomes a common time-weighting
#       on the benefit trajectory. All targeting signal comes from b and cost.
#   (2) The myopic hazard belief becomes CORRECT (hazard genuinely is constant),
#       so myopic's only remaining error is on benefit.
#   (3) The stationary-null patch in 07__ilp_core.R (lam[,1] <- lam[,2], ~line 208)
#       becomes a no-op — the input is already flat. The near-degenerate objective
#       that 09's P3 works around should disappear, and P3 may pass on its
#       stronger original form.
# Greedy never reads trans_prob, so its SCHEDULE is unchanged across the sweep;
# only its evaluated J moves, because the evaluation S changed. Useful control.
#
# ── RATES ─────────────────────────────────────────────────────────────────────
# Kemink et al. (2023, Conservation Science and Practice 5:e12939) report PPR
# wetland drainage rates of 0.09–0.57%/yr, summarizing Oslund et al. (2010),
# Doherty et al. (2013), and Dahl (2014). Note this is a SECONDARY attribution;
# Doherty et al. (2013, Wildlife Society Bulletin 37:546–563) itself gives the
# wetland range as 0.05–0.57%/yr, and additionally reports grassland loss of
# 0.4–1.3%/yr.
#
# DENOMINATOR CAVEAT. prop_suitable (script 02) is a GRASS + WETLAND composite
# (FOREsce classes 19 perennial grass, 20 open water, 26 grassland, 28 woody
# wetland, 29 herbaceous wetland) — not wetland alone. The wetland rates are
# therefore conservative for this denominator; the grassland rates are the more
# apt bound for the (likely dominant) grass fraction. The sweep brackets both
# rather than committing to either.
#
# INTERPRETIVE CAVEAT. A published annual AREAL loss rate is being read as a
# per-parcel per-decade Bernoulli conversion probability. That units gap is the
# same one script 05 carries; it does not vanish here, it only becomes uniform
# and citable rather than data-derived and noisy.
#
# ── CONVERSION ────────────────────────────────────────────────────────────────
#   p_decade = 1 - (1 - p_annual)^10
# Compounding rather than 10 * p_annual, to match the multiplicative survival
# product in 07 (S[t] = S[t-1] * (1 - lam[t-1])). The two agree to ~1e-5 at the
# wetland bounds and diverge meaningfully only at the grassland high rate.
#
# ── POSITION IN PIPELINE ──────────────────────────────────────────────────────
# Runs AFTER 06, despite the name. It reads the COMPLETE panel (cost already
# populated) and swaps one independent column, so re-running 06's raster work
# would be wasted time. Overwrites input_data/data_panel.rds in place.
#
# IDEMPOTENT: replaces trans_prob wholesale, so it can be re-run at a new rate
# directly on top of an already-swapped panel. Re-run 05 + 06 only to restore the
# data-derived version.
#
# INPUTS:  input_data/data_panel.rds  (post-06: trans_prob AND cost populated)
# OUTPUTS: input_data/data_panel.rds  (trans_prob replaced; `uniform_rate` attr)
#          input_data/data_panel.csv
#          input_data/risk_provenance.txt   (which risk model the panel carries)
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# ══ PARAMETERS ════════════════════════════════════════════════════════════════####
# The sweep. Set RATE_LABEL to one of these and run; repeat for each. One rate per
# run, because there is one canonical data_panel.rds. Archive output_data/ between
# 08 runs (see the banner this script prints at the end).
#
#   wetland_low   0.0009/yr  Kemink et al. 2023 (via Oslund 2010 / Doherty 2013 / Dahl 2014)
#   wetland_high  0.0057/yr  as above
#   grass_low     0.0040/yr  Doherty et al. 2013
#   grass_high    0.0130/yr  Doherty et al. 2013
RATE_LIBRARY <- c(
  wetland_low  = 0.0009,
  wetland_high = 0.0057,
  grass_low    = 0.0040,
  grass_high   = 0.0130
)

RATE_LABEL <- "wetland_low"          # <-- change this per run

ANNUAL_LOSS_RATE <- RATE_LIBRARY[[RATE_LABEL]]
# To use an off-library rate for a one-off, override both lines directly, e.g.:
#   RATE_LABEL <- "custom_0.20pct"; ANNUAL_LOSS_RATE <- 0.0020

PANEL_PATH <- "input_data/data_panel.rds"
CSV_PATH   <- "input_data/data_panel.csv"
PROV_PATH  <- "input_data/risk_provenance.txt"
# ══════════════════════════════════════════════════════════════════════════════


# ── 1. Load panel and guard against running before 06 ─────────────────────────####

if (!file.exists(PANEL_PATH))
  stop("Panel not found at: ", PANEL_PATH,
       "\n  Run 01-06 first. This script swaps a column in the completed panel.")

data_panel <- readRDS(PANEL_PATH)
n_rows_before <- nrow(data_panel)

required_cols <- c("eau_id", "wmd_id", "year", "rcp", "gcm",
                   "prop_suitable", "scaled_abundance", "trans_prob", "cost")
missing_cols <- setdiff(required_cols, names(data_panel))
if (length(missing_cols))
  stop("Panel is missing column(s): ", paste(missing_cols, collapse = ", "),
       "\n  This script must run AFTER 06__cost_data.R, on the complete panel.")

# cost fully populated? (the tell-tale of a pre-06 panel: cost is all NA placeholder)
if (all(is.na(data_panel$cost)))
  stop("cost is entirely NA — this panel predates 06__cost_data.R.\n",
       "  Run 06 first, then re-run this script.")
if (anyNA(data_panel$cost))
  stop("cost has ", sum(is.na(data_panel$cost)), " NA value(s). ",
       "Re-run 06 and confirm its logic checks pass before proceeding.")
if (anyNA(data_panel$scaled_abundance))
  stop("scaled_abundance has ", sum(is.na(data_panel$scaled_abundance)),
       " NA value(s). Re-run 04 and confirm its logic checks pass.")

# Report what risk model the panel is currently carrying (attr set by a prior run
# of this script; absent => data-derived from 05).
prior_rate <- attr(data_panel, "uniform_rate")
prior_lab  <- attr(data_panel, "uniform_rate_label")

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  05b — UNIFORM CONVERSION RISK (experiment)\n")
cat("══════════════════════════════════════════════════════════════════\n\n")
cat("  Panel loaded:        ", PANEL_PATH, "\n")
cat("  Rows:                ", n_rows_before, "\n")
cat("  EAUs:                ", n_distinct(data_panel$eau_id), "\n")
if (is.null(prior_rate)) {
  cat("  Incoming risk model:  data-derived (script 05)\n")
} else {
  cat(sprintf("  Incoming risk model:  ALREADY UNIFORM — '%s' (%.4f%%/yr)\n",
              prior_lab, 100 * prior_rate))
  cat("                        (swap is idempotent; overwriting)\n")
}


# ── 2. Convert annual rate to per-decade conversion probability ───────────────####

STEP_YEARS <- 10L
p_decade   <- 1 - (1 - ANNUAL_LOSS_RATE)^STEP_YEARS

years_sorted <- sort(unique(data_panel$year))
n_t           <- length(years_sorted)
terminal_year <- max(years_sorted)

# Sanity: the panel should be on a decadal step matching STEP_YEARS.
step_check <- unique(diff(years_sorted))
if (length(step_check) != 1L || step_check != STEP_YEARS)
  stop("Panel year steps are ", paste(step_check, collapse = "/"),
       " but STEP_YEARS = ", STEP_YEARS,
       ". The annual->decadal conversion assumes a uniform ", STEP_YEARS,
       "-year step.")

cat(sprintf("\n  Rate label:          %s\n", RATE_LABEL))
cat(sprintf("  Annual loss rate:    %.4f%%/yr\n", 100 * ANNUAL_LOSS_RATE))
cat(sprintf("  Decadal trans_prob:  %.6f  (%.4f%% per decade)\n",
            p_decade, 100 * p_decade))
cat(sprintf("  Terminal year:       %d  (trans_prob forced to 0)\n", terminal_year))


# ── 3. BEFORE: distribution of the incoming trans_prob ────────────────────────####
# Captured before the overwrite. This is the comparison that gives the experiment
# its context: how far the uniform rate sits from what the FOREsce deltas implied.

cat("\n──────────────────────────────────────────────────────────────────\n")
cat("  BEFORE — incoming trans_prob by RCP x decade (non-terminal)\n")
cat("──────────────────────────────────────────────────────────────────\n\n")

before_summary <- data_panel %>%
  filter(year < terminal_year) %>%
  distinct(eau_id, year, rcp, trans_prob) %>%
  group_by(rcp, year) %>%
  summarise(
    n         = dplyr::n(),
    mean_tp   = mean(trans_prob),
    median_tp = median(trans_prob),
    min_tp    = min(trans_prob),
    max_tp    = max(trans_prob),
    .groups   = "drop"
  ) %>%
  arrange(rcp, year)

print(as.data.frame(before_summary %>%
        mutate(across(ends_with("_tp"), ~ signif(.x, 4)))),
      row.names = FALSE)

before_overall <- data_panel %>%
  filter(year < terminal_year) %>%
  distinct(eau_id, year, rcp, trans_prob) %>%
  pull(trans_prob)

cat(sprintf("\n  Pooled non-terminal trans_prob: mean %.3e | median %.3e | max %.3e\n",
            mean(before_overall), median(before_overall), max(before_overall)))
cat(sprintf("  Uniform replacement:            %.3e\n", p_decade))
cat(sprintf("  Ratio (uniform / incoming median): %.1fx\n",
            p_decade / max(median(before_overall), .Machine$double.eps)))


# ── 4. Replace trans_prob — uniform, every parcel, every scenario ─────────────####
# Uniform means uniform: no prop_suitable == 0 exception, no per-RCP or per-GCM
# variation, baseline and stationary rows included.
#
# The ONE retained exception is structural, not substantive: the terminal period's
# hazard is forced to 0. 07's survival recursion (S[t] = S[t-1] * (1 - lam[t-1]))
# never reads the terminal hazard, so a non-zero value there would be silently
# ignored — the panel would claim a risk the model does not use. Zeroing it keeps
# the panel honest about what was actually run, and preserves 05's convention so
# 09's checks stay meaningful.

data_panel <- data_panel %>%
  mutate(trans_prob = if_else(year == terminal_year, 0, p_decade))


# ── 5. Implied survival under the uniform rate ────────────────────────────────####
# Mirrors compute_survival_matrix() in 07 exactly: S[1] = 1, S[t] = S[t-1] * (1 - lam[t-1]),
# with lam flat at p_decade for all non-terminal periods. Identical for every EAU
# by construction — that is the point of the experiment.

lam_vec <- c(rep(p_decade, n_t - 1), 0)
S <- numeric(n_t); S[1] <- 1
if (n_t >= 2) for (t in 2:n_t) S[t] <- S[t - 1] * (1 - lam_vec[t - 1])

survival_curve <- data.frame(
  year            = years_sorted,
  period_idx      = seq_len(n_t) - 1L,
  survival_S      = round(S, 6),
  cumulative_loss = round(1 - S, 6),
  pct_lost        = sprintf("%.2f%%", 100 * (1 - S))
)

cat("\n──────────────────────────────────────────────────────────────────\n")
cat("  AFTER — implied survival of an UNPROTECTED parcel\n")
cat("──────────────────────────────────────────────────────────────────\n\n")
print(survival_curve, row.names = FALSE)
cat(sprintf("\n  Cumulative loss probability by %d: %.2f%%\n",
            terminal_year, 100 * (1 - S[n_t])))
cat("  (1 - S is the ONLY weight through which risk enters V in 07. Identical\n")
cat("   for every parcel here, so it sets the SCALE and TIME PROFILE of V but\n")
cat("   contributes no cross-parcel targeting signal.)\n")


# ── 6. Logic checks ───────────────────────────────────────────────────────────####

cat("\n========================================\n")
cat("  LOGIC CHECK: uniform trans_prob\n")
cat("========================================\n")

n_na  <- sum(is.na(data_panel$trans_prob))
n_oob <- sum(data_panel$trans_prob < 0 | data_panel$trans_prob > 1, na.rm = TRUE)

terminal_nonzero <- data_panel %>% filter(year == terminal_year, trans_prob != 0)

nonterm_vals <- data_panel %>% filter(year < terminal_year) %>% pull(trans_prob)
n_distinct_nonterm <- n_distinct(nonterm_vals)

# uniformity across every stratum we care about
strata_var <- data_panel %>%
  filter(year < terminal_year) %>%
  group_by(rcp, gcm, year) %>%
  summarise(n_unique = n_distinct(trans_prob), .groups = "drop") %>%
  filter(n_unique > 1)

checks <- list(
  "No NAs in trans_prob"                          = n_na == 0,
  "All trans_prob values in [0, 1]"               = n_oob == 0,
  "Terminal year trans_prob always 0"             = nrow(terminal_nonzero) == 0,
  "Exactly one distinct non-terminal value"       = n_distinct_nonterm == 1L,
  "Non-terminal value equals p_decade"            = isTRUE(all.equal(unique(nonterm_vals), p_decade)),
  "Uniform within every rcp x gcm x year stratum" = nrow(strata_var) == 0,
  "Row count unchanged"                           = nrow(data_panel) == n_rows_before,
  "cost still fully populated"                    = !anyNA(data_panel$cost),
  "scaled_abundance still fully populated"        = !anyNA(data_panel$scaled_abundance)
)

for (nm in names(checks)) {
  cat(sprintf("  %s  %s\n", if (isTRUE(checks[[nm]])) "[PASS]" else "[FAIL]", nm))
}

failures <- names(checks)[!vapply(checks, isTRUE, logical(1))]

if (length(failures) > 0) {

  cat("\n  --- Diagnostic detail ---\n")
  if (!isTRUE(checks[["No NAs in trans_prob"]]))
    cat("  NA count:", n_na, "\n")
  if (!isTRUE(checks[["All trans_prob values in [0, 1]"]]))
    cat("  Out-of-range count:", n_oob, "\n")
  if (!isTRUE(checks[["Terminal year trans_prob always 0"]]))
    cat("  Non-zero terminal rows:", nrow(terminal_nonzero), "\n")
  if (!isTRUE(checks[["Exactly one distinct non-terminal value"]])) {
    cat("  Distinct non-terminal values:", n_distinct_nonterm, "\n")
    print(head(sort(unique(nonterm_vals)), 10))
  }
  if (!isTRUE(checks[["Non-terminal value equals p_decade"]]))
    cat("  Found:", unique(nonterm_vals), "| expected:", p_decade, "\n")
  if (!isTRUE(checks[["Uniform within every rcp x gcm x year stratum"]])) {
    cat("  Strata with >1 distinct value:", nrow(strata_var), "\n")
    print(head(strata_var, 10))
  }
  if (!isTRUE(checks[["Row count unchanged"]]))
    cat("  Before:", n_rows_before, "| after:", nrow(data_panel), "\n")

  cat("========================================\n\n")
  stop("Logic check FAILED: uniform trans_prob swap has errors. ",
       "Panel NOT saved.\n",
       "Failed checks: ", paste(failures, collapse = ", "))
}

cat(sprintf("\n  All checks passed. trans_prob = %.6f | Rows: %d\n",
            p_decade, nrow(data_panel)))
cat("========================================\n")


# ── 7. Save ───────────────────────────────────────────────────────────────────####
# Provenance is stamped as an RDS attribute (survives readRDS in 08/09, invisible
# to every downstream select()) and mirrored to a plain-text sidecar, because the
# CSV cannot carry it and the panel filename does not change.

attr(data_panel, "uniform_rate")        <- ANNUAL_LOSS_RATE
attr(data_panel, "uniform_rate_label")  <- RATE_LABEL
attr(data_panel, "uniform_p_decade")    <- p_decade
attr(data_panel, "risk_model")          <- "uniform_literature_05b"
attr(data_panel, "risk_stamped_at")     <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

saveRDS(data_panel, PANEL_PATH)
write_csv(data_panel, CSV_PATH)

writeLines(c(
  "PANEL RISK PROVENANCE",
  "---------------------",
  sprintf("risk_model:    uniform_literature (05b__risk_uniform.R)"),
  sprintf("rate_label:    %s", RATE_LABEL),
  sprintf("annual_rate:   %.6f  (%.4f%%/yr)", ANNUAL_LOSS_RATE, 100 * ANNUAL_LOSS_RATE),
  sprintf("p_decade:      %.8f  (%.4f%% per decade)", p_decade, 100 * p_decade),
  sprintf("terminal_year: %d (trans_prob = 0)", terminal_year),
  sprintf("cum_loss_%d:  %.4f", terminal_year, 1 - S[n_t]),
  sprintf("stamped_at:    %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Source: Kemink et al. 2023 (Conserv Sci Pract 5:e12939), summarizing",
  "Oslund et al. 2010, Doherty et al. 2013, Dahl 2014. Grassland rates from",
  "Doherty et al. 2013 (Wildl Soc Bull 37:546-563).",
  "",
  "To restore the data-derived risk model: re-run 05__risk_of_loss.R then",
  "06__cost_data.R, and delete this file."
), PROV_PATH)

cat(sprintf("\n  Saved: %s  (uniform_rate attr = %.4f)\n", PANEL_PATH, ANNUAL_LOSS_RATE))
cat(sprintf("         %s\n", CSV_PATH))
cat(sprintf("         %s\n", PROV_PATH))


# ── 8. Run reminder ───────────────────────────────────────────────────────────####

cat("\n")
cat("  ╔════════════════════════════════════════════════════════════════╗\n")
cat("  ║  ARCHIVE output_data/ BEFORE RUNNING 08                         ║\n")
cat("  ╠════════════════════════════════════════════════════════════════╣\n")
cat("  ║  08 writes model_results.csv / model_trajectories.csv /         ║\n")
cat("  ║  acquisition_schedule_spatial.csv to fixed paths. Running it    ║\n")
cat("  ║  now WILL overwrite your data-derived results, and those are    ║\n")
cat("  ║  NOT cheap to regenerate (Gurobi, 5 scenarios x 3 models).      ║\n")
cat("  ║                                                                 ║\n")
cat("  ║    file.rename('output_data', 'output_data_<label>')            ║\n")
cat("  ║                                                                 ║\n")
cat("  ║  ...then run 08. Repeat per rate in the sweep.                  ║\n")
cat("  ╚════════════════════════════════════════════════════════════════╝\n")
cat(sprintf("\n  Suggested archive name for the CURRENT contents: output_data_%s\n",
            if (is.null(prior_lab)) "data_derived" else prior_lab))
cat(sprintf("  This run's results will belong in:               output_data_%s\n\n",
            RATE_LABEL))

cat("  Next: [archive] -> 09__validation.R (P3 should now be non-degenerate) -> 08__run_models.R\n\n")
