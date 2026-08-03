# ══════════════════════════════════════════════════════════════════════════════
# 12_diagnostic_gross_vs_net.R  —  Does trans_prob measure what it should?
# ══════════════════════════════════════════════════════════════════════════════
#
# STANDALONE, READ-ONLY DIAGNOSTIC. Writes only to derived_data/ and output_figs/
# under names prefixed "diag_". Touches nothing the pipeline reads. Safe to run at
# any time; safe to delete afterwards.
#
# ── THE QUESTION ──────────────────────────────────────────────────────────────
# Script 05 computes the per-period hazard as a NET change in habitat stock:
#
#       trans_raw = (ps_t - ps_{t+1}) / ps_t
#
# The intended quantity is the GROSS transition rate: the fraction of habitat
# suitable at t that is unsuitable at t+1. Writing L for area converting
# suitable -> unsuitable and G for area converting unsuitable -> suitable,
#
#       ps_{t+1} = ps_t - L + G     =>    (ps_t - ps_{t+1}) / ps_t = (L - G) / ps_t
#
# so script 05 computes (L - G)/ps_t where L/ps_t was intended. The two agree only
# where G = 0. This script recovers L and G separately from the FOREsce rasters --
# which requires going back to pixel level, because the aggregation in script 02
# nets them irreversibly -- and quantifies the difference.
#
# ── WHAT IT REPORTS ───────────────────────────────────────────────────────────
#   A. FIDELITY   — does the recomputed prop_suitable reproduce the panel's? Does
#                   the recomputed net reproduce the panel's pre-floor trans_prob?
#                   (If A fails, nothing downstream in this script is trustworthy.)
#   B. CENSORING  — how many EAU-years have gross loss > 0 but net <= 0, and are
#                   therefore floored to epsilon and enter the model as risk-free?
#   C. MAGNITUDE  — gross vs net distributions, by decade and RCP.
#   D. LITERATURE — both metrics converted to the per-decade basis the model uses,
#                   set against published PPR conversion rates.
#   E. ANOMALY    — probe for the undiagnosed 2050 dip in per-period risk.
#
# ── INPUTS ────────────────────────────────────────────────────────────────────
#   input_data/prairie_potholes_gcam_ref_rcp45/*.tif
#   input_data/prairie_potholes_gcam_ref_rcp85/*.tif
#   input_data/prairie_potholes_gcam_ref_meanclim/..._meanclim_2020.tif
#   derived_data/foresce_eau_shared_mask.tif      (script 02)
#   derived_data/eau_wmd_lookup.csv               (script 01)
#   derived_data/data_panel.rds                   (script 06)  [for section A only]
#
# ── OUTPUTS ───────────────────────────────────────────────────────────────────
#   derived_data/diag_gross_vs_net.csv    long: eau_id x year x rcp, all metrics
#   derived_data/diag_summary.csv         the per-decade summary table
#   output_figs/diag_gross_vs_net_*.png   optional figures (DRAW_FIGS)
#
# ── RUNTIME ───────────────────────────────────────────────────────────────────
# Script 02 does 10 resamples per RCP in ~4 min. This does ~26 (one per year for
# the denominator, plus lost/gained per year-pair), so budget ~10 min per RCP,
# ~20 min total. Set QUICK <- TRUE for a 2-pair smoke test first -- do that before
# committing to a full run.
# ══════════════════════════════════════════════════════════════════════════════

if (!exists(".SETUP_DONE")) source("00_setup.R")

if (!HAS_SPATIAL)
  stop("This diagnostic needs terra + sf. Install them, or run it on a machine ",
       "that has the spatial stack.")

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(tidyr); library(readr); library(stringr)
})


# ══ PARAMETERS ══════════════════════════════════════════════════════════════####

# Smoke test: process only the first two year-pairs of RCP 4.5. Confirms the
# machinery and the section-A fidelity checks in ~2 min before the full run.
QUICK <- TRUE

# Suitable-habitat classes. MUST match script 02 exactly, or section A will fail.
#   19 = Perennial Grass, 20 = Open Water, 26 = Grassland,
#   28 = Woody Wetland,   29 = Herbaceous Wetland
SUITABLE_CLASSES <- c(19, 20, 26, 28, 29)

# Tolerance for the fidelity checks in section A. prop_suitable is a ratio of
# pixel counts, so agreement should be near-exact; 1e-6 leaves room for float
# accumulation across ~3e5 pixels per EAU without masking a real discrepancy.
FIDELITY_TOL <- 1e-6

# Published PPR conversion rates, ANNUAL, for the section-D comparison.
LIT_RATES <- tibble::tribble(
  ~source,                        ~habitat,    ~annual_rate,
  "Doherty et al. 2018",          "wetland",   0.0057,
  "Doherty et al. 2018",          "grassland", 0.0133,
  "Rashford et al. 2011",         "grassland", 0.0133,
  "Kemink et al. 2023 (upper)",   "wetland",   0.0057,
  "Kemink et al. 2023 (lower)",   "wetland",   0.0009
)

DRAW_FIGS <- TRUE     # write diagnostic PNGs to output_figs/

# Paths (mirroring scripts 02 and 06)
LU_DIRS <- list(`45` = "input_data/prairie_potholes_gcam_ref_rcp45",
                `85` = "input_data/prairie_potholes_gcam_ref_rcp85")
MEANCLIM_FILE <- file.path("input_data/prairie_potholes_gcam_ref_meanclim",
                           "prairie_potholes_gcam_ref_meanclim_2020.tif")
SHARED_MASK <- file.path(DIR_DERIVED, "foresce_eau_shared_mask.tif")
EAU_LOOKUP  <- file.path(DIR_DERIVED, "eau_wmd_lookup.csv")
PANEL_RDS   <- file.path(DIR_DERIVED, "data_panel.rds")
# ══════════════════════════════════════════════════════════════════════════════


hdr <- function(txt) cat(sprintf("\n%s\n  %s\n%s\n",
                                 strrep("=", 78), txt, strrep("=", 78)))

hdr("12_diagnostic_gross_vs_net.R")
if (QUICK) cat("  *** QUICK MODE: two year-pairs, RCP 4.5 only. Not a full run. ***\n")


# ── 0. Load the EAU grid and lookup ───────────────────────────────────────####

for (f in c(SHARED_MASK, EAU_LOOKUP)) if (!file.exists(f))
  stop("Required input not found: ", f,
       "\n  Run scripts 01 and 02 first.")

aoi_mask <- rast(SHARED_MASK)
eau_wmd  <- read_csv(EAU_LOOKUP, show_col_types = FALSE)

# Scale conversion factor: number of native FOREsce pixels per EAU cell.
# Identical construction to script 02, so the denominators match.
scf     <- (res(aoi_mask)[1]^2) / (30^2)
n_cells <- ncell(aoi_mask)

cat(sprintf("  EAU grid: %d x %d cells at %.0f m (%.0f native pixels per EAU)\n",
            nrow(aoi_mask), ncol(aoi_mask), res(aoi_mask)[1], scf))
cat(sprintf("  EAUs in lookup: %d\n", nrow(eau_wmd)))

# cell_id -> eau_id. Script 03 joins land cover to the panel BY CELL NUMBER
# (the post-7/30 by-place join); we reproduce that keying exactly.
cell_to_eau <- eau_wmd %>% select(cell_id, eau_id, wmd_id)


# ── 1. Per-year-pair flow decomposition ───────────────────────────────────####
#
# Everything happens at native FOREsce resolution BEFORE aggregation, which is the
# whole point -- script 02's resample-then-difference destroys the information that
# separates L from G.
#
#   suit_t  = 1 where suitable at t,   0 otherwise,  NA outside coverage
#   lost    = suit_t * (1 - suit_t1)   -> 1 exactly where suitable -> unsuitable
#   gained  = (1 - suit_t) * suit_t1   -> 1 exactly where unsuitable -> suitable
#
# Each is then summed to EAU resolution and divided by scf, giving the fraction of
# the EAU in each category -- the same units as prop_suitable, so the three are
# directly comparable and the stock identity below can be checked exactly.
#
# NOTE ON `others`: script 02 uses others = NA; we use others = 0. Summing with
# NA-skip and summing zeros give the same count of 1s, so the two agree EXCEPT for
# an EAU with no suitable pixels at all, which yields NA under 02's encoding and 0
# here. Section A detects any such case -- it is a candidate explanation for part
# of the 41-EAU coverage drop, and worth knowing about either way.

rcl_suitable <- cbind(SUITABLE_CLASSES, 1)

classify_suitable <- function(path) {
  r <- rast(path)
  r <- crop(r, aoi_mask)
  classify(r, rcl = rcl_suitable, others = 0)
}

# Aggregate a native-resolution 0/1 raster to EAU proportions.
to_eau_prop <- function(r_native) {
  agg <- resample(r_native, aoi_mask, method = "sum")
  mask(agg / scf, aoi_mask)
}

# Returns a data.frame keyed by cell_id with the three flow proportions.
decompose_pair <- function(path_t, path_t1, label) {
  t0 <- Sys.time()
  
  suit_t  <- classify_suitable(path_t)
  suit_t1 <- classify_suitable(path_t1)
  
  prop_t      <- to_eau_prop(suit_t)
  prop_lost   <- to_eau_prop(suit_t  * (1 - suit_t1))
  prop_gained <- to_eau_prop((1 - suit_t) * suit_t1)
  prop_t1     <- to_eau_prop(suit_t1)
  
  out <- data.frame(
    cell_id     = seq_len(n_cells),
    ps_t        = values(prop_t,      na.rm = FALSE)[, 1],
    ps_t1       = values(prop_t1,     na.rm = FALSE)[, 1],
    prop_lost   = values(prop_lost,   na.rm = FALSE)[, 1],
    prop_gained = values(prop_gained, na.rm = FALSE)[, 1]
  )
  
  cat(sprintf("    %-28s %5.1f s\n", label,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  out
}


# ── 2. Enumerate the year-pairs to process ────────────────────────────────####
#
# Mirrors script 05's structure:
#   * 2020 -> 2030 crosses the meanclim/RCP boundary (the 2020 decision is made
#     before any RCP has materialised), so ps_start comes from meanclim.
#   * 2030 -> 2100 are within-RCP forward differences.
#   * 2100 is terminal: no subsequent period, hazard is 0 by construction, excluded.

list_rcp_files <- function(dir, suffix) {
  fn <- list.files(dir, pattern = paste0("prairie_potholes_gcam_ref_rcp", suffix,
                                         "_[0-9]{4}\\.tif$"), full.names = TRUE)
  if (!length(fn)) stop("No rasters matched in ", dir)
  yr <- as.integer(str_extract(basename(fn), "[0-9]{4}(?=\\.tif$)"))
  o  <- order(yr)
  data.frame(path = fn[o], year = yr[o], stringsAsFactors = FALSE)
}

build_pairs <- function(suffix) {
  fl <- list_rcp_files(LU_DIRS[[suffix]], suffix)
  fl <- fl[fl$year >= 2030, ]          # 2014/2020 RCP layers unused by the panel
  
  pairs <- data.frame(
    rcp        = suffix,
    year       = head(fl$year, -1),    # the period the hazard is attached to
    path_t     = head(fl$path, -1),
    path_t1    = tail(fl$path, -1),
    stringsAsFactors = FALSE
  )
  
  # prepend the 2020 -> 2030 step, anchored on meanclim
  if (file.exists(MEANCLIM_FILE)) {
    pairs <- rbind(
      data.frame(rcp = suffix, year = 2020L,
                 path_t = MEANCLIM_FILE, path_t1 = fl$path[fl$year == 2030],
                 stringsAsFactors = FALSE),
      pairs)
  } else {
    warning("meanclim raster not found -- skipping the 2020 -> 2030 step.")
  }
  pairs
}

pair_tbl <- bind_rows(lapply(names(LU_DIRS), build_pairs))
if (QUICK) pair_tbl <- pair_tbl %>% filter(rcp == "45") %>% head(2)

cat(sprintf("\n  Year-pairs to process: %d\n", nrow(pair_tbl)))


# ── 3. Run the decomposition ──────────────────────────────────────────────####

hdr("Processing rasters (this is the slow part)")
t_all <- Sys.time()

flows <- vector("list", nrow(pair_tbl))
for (k in seq_len(nrow(pair_tbl))) {
  p <- pair_tbl[k, ]
  cat(sprintf("  RCP %s  %d -> next\n", p$rcp, p$year))
  d <- decompose_pair(p$path_t, p$path_t1, sprintf("rcp%s %d", p$rcp, p$year))
  d$rcp  <- p$rcp
  d$year <- p$year
  flows[[k]] <- d
}

flows <- bind_rows(flows) %>%
  inner_join(cell_to_eau, by = "cell_id") %>%     # by-place join, as in script 03
  filter(!is.na(ps_t))

cat(sprintf("\n  Total: %.1f min | EAU-year rows: %d | distinct EAUs: %d\n",
            as.numeric(difftime(Sys.time(), t_all, units = "mins")),
            nrow(flows), n_distinct(flows$eau_id)))


# ── 4. Derive the competing metrics ───────────────────────────────────────####

metrics <- flows %>%
  mutate(
    # stock identity: ps_{t+1} should equal ps_t - lost + gained, exactly
    ident_resid = ps_t1 - (ps_t - prop_lost + prop_gained),
    
    gross_loss  = ifelse(ps_t > 0, prop_lost   / ps_t, 0),  # what we want
    gross_gain  = ifelse(ps_t > 0, prop_gained / ps_t, 0),  # what netting hides
    net_loss    = ifelse(ps_t > 0, (ps_t - ps_t1) / ps_t, 0),  # what 05 computes
    
    censored    = gross_loss > 0 & net_loss <= 0,
    understated = gross_loss > net_loss + 1e-12
  )


# ══ A. FIDELITY ══════════════════════════════════════════════════════════════####
hdr("A. FIDELITY -- does this reproduce the pipeline's own numbers?")

# A1: the stock identity must hold to floating-point precision. If it does not, the
# flow decomposition is wrong and everything below is void.
max_resid <- max(abs(metrics$ident_resid), na.rm = TRUE)
a1 <- max_resid < FIDELITY_TOL
cat(sprintf("  [%s] stock identity  ps_t1 == ps_t - L + G   (max resid %.2e)\n",
            ifelse(a1, "PASS", "FAIL"), max_resid))

# A2 / A3: agreement with the panel, if it is available.
a2 <- a3 <- NA
if (file.exists(PANEL_RDS)) {
  panel <- readRDS(PANEL_RDS)
  
  # A2 -- recomputed prop_suitable vs the panel's, for the RCP rows (2030+).
  # A mismatch here means this script's raster handling has drifted from script 02
  # (most likely the others = NA vs others = 0 encoding on all-unsuitable EAUs).
  panel_ps <- panel %>%
    filter(rcp %in% c("45", "85"), year >= 2030) %>%
    distinct(eau_id, year, rcp, panel_ps = prop_suitable)
  
  cmp_ps <- metrics %>%
    filter(year >= 2030) %>%
    distinct(eau_id, year, rcp, ps_t) %>%
    inner_join(panel_ps, by = c("eau_id", "year", "rcp")) %>%
    mutate(dev = abs(ps_t - panel_ps))
  
  a2 <- nrow(cmp_ps) > 0 && max(cmp_ps$dev, na.rm = TRUE) < FIDELITY_TOL
  cat(sprintf("  [%s] prop_suitable matches panel      (n = %d, max dev %.2e)\n",
              ifelse(isTRUE(a2), "PASS", "FAIL"), nrow(cmp_ps),
              if (nrow(cmp_ps)) max(cmp_ps$dev, na.rm = TRUE) else NA_real_))
  
  # A3 -- recomputed net vs the panel's trans_prob. These agree only where the
  # floor did not bind, so restrict to rows above the floor. Disagreement there
  # would indicate a genuine discrepancy in how 05 differences the series.
  panel_tp <- panel %>%
    filter(rcp %in% c("45", "85"), year >= 2030, year < 2100) %>%
    distinct(eau_id, year, rcp, trans_prob)
  
  cmp_tp <- metrics %>%
    filter(year >= 2030, year < 2100) %>%
    inner_join(panel_tp, by = c("eau_id", "year", "rcp")) %>%
    filter(net_loss > 1e-4)                       # comfortably above any epsilon
  a3 <- nrow(cmp_tp) > 0 &&
    max(abs(cmp_tp$net_loss - cmp_tp$trans_prob), na.rm = TRUE) < 1e-6
  cat(sprintf("  [%s] net reproduces trans_prob        (n = %d, max dev %.2e)\n",
              ifelse(isTRUE(a3), "PASS", "FAIL"), nrow(cmp_tp),
              if (nrow(cmp_tp)) max(abs(cmp_tp$net_loss - cmp_tp$trans_prob),
                                    na.rm = TRUE) else NA_real_))
} else {
  cat("  [SKIP] panel not found -- A2/A3 skipped. Run script 06 to enable them.\n")
}

if (!a1)
  stop("Stock identity failed. The flow decomposition is not sound -- stop and ",
       "investigate before reading any result below.")
if (isFALSE(a2))
  warning("Recomputed prop_suitable does not match the panel. Sections B-D are ",
          "still internally consistent, but reconcile this before citing them ",
          "against pipeline numbers.")


# ══ B. CENSORING ═════════════════════════════════════════════════════════════####
hdr("B. CENSORING -- real loss that the net metric erases")

nonterm <- metrics %>% filter(year < 2100)

cens <- nonterm %>%
  summarise(
    n_eau_years     = n(),
    n_any_loss      = sum(gross_loss > 0),
    n_censored      = sum(censored),
    pct_censored    = 100 * n_censored / n_eau_years,
    pct_of_losing   = 100 * n_censored / pmax(n_any_loss, 1),
    n_understated   = sum(understated),
    pct_understated = 100 * n_understated / n_eau_years
  )

cat(sprintf("  EAU-years (non-terminal):                     %8d\n", cens$n_eau_years))
cat(sprintf("  ... with ANY real habitat loss (gross > 0):   %8d  (%.1f%%)\n",
            cens$n_any_loss, 100 * cens$n_any_loss / cens$n_eau_years))
cat(sprintf("  ... CENSORED (gross > 0 but net <= 0):        %8d  (%.1f%% of all,\n",
            cens$n_censored, cens$pct_censored))
cat(sprintf("                                                           %.1f%% of losing)\n",
            cens$pct_of_losing))
cat(sprintf("  ... understated to any degree (gross > net):  %8d  (%.1f%%)\n",
            cens$n_understated, cens$pct_understated))

cat("\n  Magnitude of the loss inside the censored rows -- every one of these\n")
cat("  enters the model at the epsilon floor, i.e. as effectively risk-free:\n\n")
if (cens$n_censored > 0) {
  cz <- nonterm %>% filter(censored)
  print(cz %>% summarise(
    n        = n(),
    mean     = round(mean(gross_loss), 5),
    median   = round(median(gross_loss), 5),
    p90      = round(quantile(gross_loss, 0.90), 5),
    max      = round(max(gross_loss), 5)) %>% as.data.frame(), row.names = FALSE)
  cat(sprintf("\n  Worst censored parcel loses %.1f%% of its habitat in a decade\n",
              100 * max(cz$gross_loss)))
  cat("  and is modelled as facing the background floor.\n")
} else {
  cat("  (none -- gain is nowhere large enough to flip the sign)\n")
}

# Landscape-scale accounting: what share of total real loss is netted away?
tot <- nonterm %>%
  summarise(L = sum(prop_lost, na.rm = TRUE), G = sum(prop_gained, na.rm = TRUE)) %>%
  mutate(net = L - G, hidden_pct = 100 * G / pmax(L, 1e-12))
cat(sprintf("\n  Landscape totals (EAU-fraction units, summed over all EAU-years):\n"))
cat(sprintf("    gross loss L      = %10.2f\n", tot$L))
cat(sprintf("    gross gain G      = %10.2f\n", tot$G))
cat(sprintf("    net (L - G)       = %10.2f\n", tot$net))
cat(sprintf("    => the net metric conceals %.1f%% of all habitat conversion\n",
            tot$hidden_pct))


# ══ C. MAGNITUDE ═════════════════════════════════════════════════════════════####
hdr("C. MAGNITUDE -- gross vs net by decade and RCP")

by_decade <- nonterm %>%
  group_by(rcp, year) %>%
  summarise(
    n            = n(),
    gross_mean   = mean(gross_loss),
    gross_median = median(gross_loss),
    gross_max    = max(gross_loss),
    gain_mean    = mean(gross_gain),
    net_mean     = mean(net_loss),
    net_median   = median(net_loss),
    ratio        = gross_mean / pmax(net_mean, 1e-12),
    pct_censored = 100 * mean(censored),
    .groups = "drop"
  )

print(by_decade %>%
        mutate(across(c(gross_mean, gross_median, gross_max,
                        gain_mean, net_mean, net_median), ~ round(.x, 5)),
               ratio = round(ratio, 2), pct_censored = round(pct_censored, 1)) %>%
        as.data.frame(), row.names = FALSE)

overall <- nonterm %>%
  summarise(gross_mean = mean(gross_loss), gross_median = median(gross_loss),
            gross_max = max(gross_loss), net_mean = mean(net_loss),
            net_median = median(net_loss), net_max = max(net_loss))
cat(sprintf("\n  Overall  gross: mean %.5f  median %.5f  max %.5f\n",
            overall$gross_mean, overall$gross_median, overall$gross_max))
cat(sprintf("           net:   mean %.5f  median %.5f  max %.5f\n",
            overall$net_mean, overall$net_median, overall$net_max))
cat(sprintf("           mean gross is %.1fx the mean net\n",
            overall$gross_mean / max(overall$net_mean, 1e-12)))


# ══ D. LITERATURE ════════════════════════════════════════════════════════════####
hdr("D. LITERATURE -- like-for-like, on the model's decadal clock")

cat("  The model's time step is a DECADE, so published ANNUAL rates must be\n")
cat("  compounded before comparison:  decadal = 1 - (1 - annual)^10\n\n")

lit <- LIT_RATES %>%
  mutate(decadal_rate = 1 - (1 - annual_rate)^10) %>%
  mutate(across(c(annual_rate, decadal_rate), ~ round(.x, 4)))
print(as.data.frame(lit), row.names = FALSE)

cat("\n  Model, same units (per decade):\n\n")
model_row <- data.frame(
  metric = c("GROSS (intended)", "NET (script 05, current)"),
  mean   = round(c(overall$gross_mean,   overall$net_mean),   5),
  median = round(c(overall$gross_median, overall$net_median), 7),
  max    = round(c(overall$gross_max,    overall$net_max),    4)
)
print(model_row, row.names = FALSE)

lit_lo <- min(lit$decadal_rate); lit_hi <- max(lit$decadal_rate)
cat(sprintf("\n  Published decadal range: %.4f - %.4f\n", lit_lo, lit_hi))
cat(sprintf("  Model GROSS mean is %.2fx the low end, %.2fx the high end.\n",
            overall$gross_mean / lit_lo, overall$gross_mean / lit_hi))
cat(sprintf("  Model NET   mean is %.2fx the low end, %.2fx the high end.\n",
            overall$net_mean / lit_lo, overall$net_mean / lit_hi))
cat("\n  Read this next to the .qmd 'Conversion Risk' section, which compares the\n")
cat("  model's per-DECADE rates against per-YEAR literature values. Recomputing on\n")
cat("  a common basis is the first step; whether to then CALIBRATE the landscape\n")
cat("  total to the literature (the Kemink et al. approach) is a separate decision.\n")


# ══ E. THE 2050 ANOMALY ══════════════════════════════════════════════════════####
hdr("E. ANOMALY PROBE -- the undiagnosed 2050 dip")

cat("  If the dip is real, gross loss dips at 2050 too. If gross is smooth while\n")
cat("  net dips, the dip is a netting artefact -- a pulse of restoration masking a\n")
cat("  normal amount of conversion.\n\n")

print(by_decade %>%
        select(rcp, year, gross_mean, gain_mean, net_mean) %>%
        mutate(across(where(is.numeric), ~ round(.x, 5))) %>%
        as.data.frame(), row.names = FALSE)

# Duplicate-layer check: identical consecutive rasters would produce exactly this
# signature and is a cheap thing to rule out.
cat("\n  Consecutive-year correlation of prop_suitable (1.0000 => identical layers):\n\n")
for (sfx in unique(metrics$rcp)) {
  wide <- metrics %>% filter(rcp == sfx) %>%
    distinct(eau_id, year, ps_t) %>%
    pivot_wider(names_from = year, values_from = ps_t) %>%
    select(-eau_id)
  yrs <- names(wide)
  for (j in seq_len(ncol(wide) - 1)) {
    cc <- suppressWarnings(cor(wide[[j]], wide[[j + 1]], use = "complete.obs"))
    flag <- if (!is.na(cc) && cc > 0.99999) "  <-- effectively identical" else ""
    cat(sprintf("    rcp%s  %s -> %s   r = %.5f%s\n", sfx, yrs[j], yrs[j + 1], cc, flag))
  }
}


# ══ SAVE ═════════════════════════════════════════════════════════════════════####
hdr("Saving")

out_long <- metrics %>%
  select(eau_id, wmd_id, cell_id, year, rcp,
         ps_t, ps_t1, prop_lost, prop_gained,
         gross_loss, gross_gain, net_loss, censored, understated) %>%
  arrange(rcp, year, eau_id)

write_csv(out_long,   file.path(DIR_DERIVED, "diag_gross_vs_net.csv"))
write_csv(by_decade,  file.path(DIR_DERIVED, "diag_summary.csv"))
cat(sprintf("  %s  (%d rows)\n", file.path(DIR_DERIVED, "diag_gross_vs_net.csv"),
            nrow(out_long)))
cat(sprintf("  %s\n", file.path(DIR_DERIVED, "diag_summary.csv")))


# ── Figures (optional) ────────────────────────────────────────────────────####
if (DRAW_FIGS && requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  if (!dir.exists(DIR_FIGS)) dir.create(DIR_FIGS, recursive = TRUE)
  
  # F1 -- gross vs net scatter. Points below the diagonal are understated; points
  # at or left of zero on the x-axis are censored to the floor entirely.
  p1 <- ggplot(nonterm, aes(net_loss, gross_loss, colour = censored)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_point(alpha = 0.25, size = 0.6) +
    scale_colour_manual(values = c(`FALSE` = "#2E7D5B", `TRUE` = "#B0413E"),
                        labels = c("Retained", "Censored to floor"), name = NULL) +
    facet_wrap(~ rcp, labeller = as_labeller(c(`45` = "RCP 4.5", `85` = "RCP 8.5"))) +
    labs(x = "Net loss rate (script 05, current)",
         y = "Gross loss rate (intended)",
         title = "Every point below the dashed line understates conversion risk",
         subtitle = "Points left of the vertical line are floored to epsilon despite real habitat loss") +
    theme_minimal(base_size = 11) + theme(legend.position = "top")
  ggsave(file.path(DIR_FIGS, "diag_gross_vs_net_scatter.png"), p1,
         width = 9, height = 4.8, dpi = 200, bg = "white")
  
  # F2 -- the two metrics through time, with the literature band for scale
  p2 <- by_decade %>%
    select(rcp, year, GROSS = gross_mean, NET = net_mean) %>%
    pivot_longer(c(GROSS, NET), names_to = "metric", values_to = "rate") %>%
    ggplot(aes(year, rate, colour = metric)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = lit_lo, ymax = lit_hi,
             fill = "grey80", alpha = 0.45) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
    scale_colour_manual(values = c(GROSS = "#B0413E", NET = "#2E7D5B"), name = NULL) +
    facet_wrap(~ rcp, labeller = as_labeller(c(`45` = "RCP 4.5", `85` = "RCP 8.5"))) +
    labs(x = "Decision year", y = "Mean per-decade loss rate",
         title = "Mean conversion rate: gross vs net",
         subtitle = "Grey band = published PPR rates, compounded to a decadal basis") +
    theme_minimal(base_size = 11) + theme(legend.position = "top")
  ggsave(file.path(DIR_FIGS, "diag_gross_vs_net_time.png"), p2,
         width = 9, height = 4.8, dpi = 200, bg = "white")
  
  cat(sprintf("  %s\n  %s\n",
              file.path(DIR_FIGS, "diag_gross_vs_net_scatter.png"),
              file.path(DIR_FIGS, "diag_gross_vs_net_time.png")))
}

hdr("Done -- no pipeline file was modified")
if (QUICK) cat("  (QUICK mode: partial coverage. Set QUICK <- FALSE for the real numbers.)\n")