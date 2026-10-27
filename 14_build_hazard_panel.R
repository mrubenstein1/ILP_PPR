# ══════════════════════════════════════════════════════════════════════════════
# 14_build_hazard_panel.R — rescaled conversion hazard for the high-hazard run
# ══════════════════════════════════════════════════════════════════════════════
#
# USAGE (from the project root; the canonical run must already exist):
#   PPR_HAZARD=1.33 Rscript 14_build_hazard_panel.R
#   PPR_HAZARD=1.33 Rscript 00_run_all.R 08
#
# PPR_HAZARD is a target ANNUAL conversion rate in percent (1.33 = Rashford et al.
# 2011). This script reads the CANONICAL data_panel.rds (output of script 06) and
# writes a copy in which only trans_prob is rescaled, into the high-hazard run's own
# derived_data/. No land-cover processing (01-06) is rerun; scripts 08-10 then read
# the scaled panel exactly as they read the canonical one.
#
# ── THE TRANSFORM ─────────────────────────────────────────────────────────────
#   trans_prob' = 1 - (1 - trans_prob)^k
# This multiplies the underlying (continuous-time) conversion rate by k in every
# EAU, decade and scenario. The relative spatial pattern and the decade-to-decade
# trajectory are preserved on that rate scale; on the probability scale they are
# preserved almost exactly for typical values and compressed only near 1. Unlike
# trans_prob * k, it can never exceed 1. Zeros (no habitat, terminal 2100) stay 0.
#
# ── HOW k IS CHOSEN ───────────────────────────────────────────────────────────
# k is solved so that the HABITAT-WEIGHTED mean of the 2020 trans_prob equals the
# target compounded to a decade, 1 - (1 - r)^10. "2020 trans_prob" is the shared
# baseline row (rcp = gcm = "baseline"), which is period 0 of every scenario in the
# ILP. Weights are each EAU's unprotected suitable habitat at the start of the 2020
# step (eau_flows.rds, prop_suit_unprot), so the mean is a landscape rate: habitat
# lost / habitat present. EAUs are equal-area, so proportion weights = area weights.
# The unweighted mean (and the k it would imply) is printed and recorded for
# comparison only.
#
# OUTPUTS (runs/<LC_SCENARIO>_haz<tag>/derived_data/):
#   data_panel.rds        canonical panel with trans_prob rescaled
#   hazard_scaling.json   target, solved k, both means, source panel md5
# ══════════════════════════════════════════════════════════════════════════════

if (!exists(".SETUP_DONE")) source("00_setup.R")

if (!HAZARD_RUN)
  stop("PPR_HAZARD is not set. Usage:  PPR_HAZARD=1.33 Rscript 14_build_hazard_panel.R",
       call. = FALSE)
if (normalizePath(DIR_DERIVED, mustWork = FALSE) ==
    normalizePath(DIR_CANON_DERIVED, mustWork = FALSE))
  stop("Output folder is the canonical derived_data/; refusing to overwrite it.",
       call. = FALSE)

SRC_PANEL <- file.path(DIR_CANON_DERIVED, "data_panel.rds")
SRC_FLOWS <- file.path(DIR_CANON_DERIVED, "eau_flows.rds")
OUT_PANEL <- file.path(DIR_DERIVED, "data_panel.rds")
OUT_JSON  <- file.path(DIR_DERIVED, "hazard_scaling.json")

for (f in c(SRC_PANEL, SRC_FLOWS))
  if (!file.exists(f))
    stop("Canonical input not found: ", f, "\n  Run the canonical pipeline first ",
         "(PPR_HAZARD unset): Rscript 00_run_all.R", call. = FALSE)

# ── 1. Load canonical inputs ──────────────────────────────────────────────────####
src   <- readRDS(SRC_PANEL)
flows <- readRDS(SRC_FLOWS)
if (anyNA(src$trans_prob))
  stop("The canonical panel has NA trans_prob values; fix script 05 first.", call. = FALSE)

# ── 2. Target ─────────────────────────────────────────────────────────────────####
r_annual       <- HAZARD_TARGET_ANNUAL_PCT / 100
target_decadal <- 1 - (1 - r_annual)^10

# ── 3. 2020 hazard and habitat weights ────────────────────────────────────────####
base20 <- src %>%
  filter(year == 2020, rcp == "baseline", gcm == "baseline") %>%
  select(eau_id, prop_suitable, trans_prob)
if (anyDuplicated(base20$eau_id))
  stop("The 2020 baseline slice has more than one row per EAU.", call. = FALSE)

# Weights: unprotected suitable habitat at the start of the 2020 step. Both RCPs'
# 2020 flows were measured against the same meanclim 2020 raster, so the weights
# must agree across RCPs. Check that rather than assume it.
w20 <- flows %>%
  filter(year == 2020L, rcp %in% c("45", "85")) %>%
  select(eau_id, rcp, prop_suit_unprot) %>%
  tidyr::pivot_wider(names_from = rcp, values_from = prop_suit_unprot,
                     names_prefix = "w_")
w_diff <- max(abs(w20$w_45 - w20$w_85))
if (!is.finite(w_diff) || w_diff > 1e-9)
  stop(sprintf("2020 unprotected habitat differs between RCP 4.5 and 8.5 (max |diff| = %.3e).",
               w_diff), " Expected identical weights; check script 02 section 5.",
       call. = FALSE)

base20 <- base20 %>%
  left_join(w20 %>% select(eau_id, w = w_45), by = "eau_id")
if (anyNA(base20$w))
  stop(sum(is.na(base20$w)), " panel EAU(s) have no 2020 flow weight.", call. = FALSE)
if (sum(base20$w) <= 0)
  stop("Total 2020 unprotected habitat is zero; cannot form a weighted mean.", call. = FALSE)

# ── 4. Transform and solve for k ──────────────────────────────────────────────####
# 1 - exp(k * log1p(-p)) is 1 - (1 - p)^k, computed accurately for small p.
scale_hazard <- function(p, k) 1 - exp(k * log1p(-p))

w_mean  <- function(p, w) sum(w * p) / sum(w)
hab     <- base20$prop_suitable > 0          # EAUs with habitat (unweighted mean)

solve_k <- function(mean_fun) {
  f  <- function(k) mean_fun(k) - target_decadal
  lo <- 1e-6; hi <- 1e6
  if (f(lo) > 0 || f(hi) < 0)
    stop("Target ", signif(target_decadal, 4), " per decade cannot be reached by ",
         "rescaling this panel.", call. = FALSE)
  uniroot(f, c(lo, hi), tol = 1e-12)$root
}

k_weighted   <- solve_k(function(k) w_mean(scale_hazard(base20$trans_prob, k), base20$w))
k_unweighted <- solve_k(function(k) mean(scale_hazard(base20$trans_prob[hab], k)))
k <- k_weighted                                   # the run uses the weighted solution

src_w_mean  <- w_mean(base20$trans_prob, base20$w)
src_uw_mean <- mean(base20$trans_prob[hab])

# ── 5. Apply to every row ─────────────────────────────────────────────────────####
out <- src %>% mutate(trans_prob = scale_hazard(trans_prob, k))

out20 <- out %>%
  filter(year == 2020, rcp == "baseline", gcm == "baseline") %>%
  select(eau_id, trans_prob) %>%
  left_join(base20 %>% select(eau_id, w, prop_suitable), by = "eau_id")
out_w_mean  <- w_mean(out20$trans_prob, out20$w)
out_uw_mean <- mean(out20$trans_prob[out20$prop_suitable > 0])

# ── 6. Checks ─────────────────────────────────────────────────────────────────####
stat_flat <- out %>%
  filter(rcp == "stationary", year < 2100) %>%
  group_by(eau_id) %>%
  summarise(n = n_distinct(trans_prob), .groups = "drop")

k1_dev <- max(abs(scale_hazard(src$trans_prob, 1) - src$trans_prob))

checks <- list(
  "Scaled weighted 2020 mean equals target"     = abs(out_w_mean - target_decadal) < 1e-8,
  "No NAs in scaled trans_prob"                 = !anyNA(out$trans_prob),
  "All scaled trans_prob in [0, 1]"             = all(out$trans_prob >= 0 & out$trans_prob <= 1),
  "Zeros stay zero (no habitat, terminal)"      = all(out$trans_prob[src$trans_prob == 0] == 0),
  "Terminal year (2100) trans_prob always 0"    = all(out$trans_prob[out$year == 2100] == 0),
  "Stationary hazard still flat within EAU"     = all(stat_flat$n == 1),
  "Rescaling moves every value the same way"    = if (k >= 1) all(out$trans_prob >= src$trans_prob)
                                                  else all(out$trans_prob <= src$trans_prob),
  "All other columns identical to canonical"    = identical(select(out, -trans_prob),
                                                            select(src, -trans_prob)),
  "k = 1 reproduces the source (to 1e-12)"      = k1_dev < 1e-12
)

cat("\n========================================\n")
cat("  LOGIC CHECK: rescaled hazard panel\n")
cat("========================================\n")
for (nm in names(checks))
  cat(sprintf("  %s  %s\n", if (checks[[nm]]) "[PASS]" else "[FAIL]", nm))
failures <- names(checks)[!unlist(checks)]
if (length(failures))
  stop("Logic check FAILED. Nothing was written.\n  Failed checks: ",
       paste(failures, collapse = ", "), call. = FALSE)

# ── 7. Report ─────────────────────────────────────────────────────────────────####
cat(sprintf("\n  Target: %.4g%% per year  ->  %.4f per decade\n",
            HAZARD_TARGET_ANNUAL_PCT, target_decadal))
cat(sprintf("  Canonical 2020 mean:  weighted %.5f | unweighted %.5f (EAUs with habitat)\n",
            src_w_mean, src_uw_mean))
cat(sprintf("  Solved k:             weighted %.4f (USED) | unweighted %.4f\n",
            k_weighted, k_unweighted))
cat(sprintf("  Scaled 2020 mean:     weighted %.5f | unweighted %.5f\n",
            out_w_mean, out_uw_mean))

cat("\n--- Trajectory: mean and max trans_prob by RCP x year, canonical vs scaled ---\n\n")
traj <- bind_rows(src %>% mutate(run = "canonical"), out %>% mutate(run = "scaled")) %>%
  filter(rcp %in% c("45", "85"), year < 2100) %>%
  distinct(run, eau_id, year, rcp, trans_prob) %>%
  group_by(rcp, year, run) %>%
  summarise(mean = mean(trans_prob), max = max(trans_prob), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = run, values_from = c(mean, max)) %>%
  mutate(ratio_of_means = mean_scaled / mean_canonical) %>%
  mutate(across(where(is.double), ~ round(.x, 4)))
print(as.data.frame(traj), row.names = FALSE)

# ── 8. Save ───────────────────────────────────────────────────────────────────####
saveRDS(out, OUT_PANEL)
jsonlite::write_json(
  list(generated                 = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
       lc_scenario               = LC_SCENARIO,
       run_id                    = RUN_ID,
       target_annual_pct         = HAZARD_TARGET_ANNUAL_PCT,
       target_decadal            = target_decadal,
       transform                 = "trans_prob' = 1 - (1 - trans_prob)^k",
       weighting                 = "unprotected suitable habitat at start of 2020 step (eau_flows prop_suit_unprot)",
       k                         = k_weighted,
       k_unweighted              = k_unweighted,
       canonical_mean_2020_weighted   = src_w_mean,
       canonical_mean_2020_unweighted = src_uw_mean,
       scaled_mean_2020_weighted      = out_w_mean,
       scaled_mean_2020_unweighted    = out_uw_mean,
       source_panel              = SRC_PANEL,
       source_panel_md5          = unname(tools::md5sum(SRC_PANEL))),
  OUT_JSON, auto_unbox = TRUE, pretty = TRUE, digits = NA)

cat(sprintf("\n✓ Scaled panel (k = %.4f) written to %s\n  Next:  PPR_HAZARD=%s Rscript 00_run_all.R 08\n",
            k, OUT_PANEL, .hazard_env))
