# ══════════════════════════════════════════════════════════════════════════════
# 11__compare_risk_models.R  —  Compare model results across the 5 risk models
# ══════════════════════════════════════════════════════════════════════════════
#
# Reads model_results.csv from the five archived run folders and produces tidy
# comparison tables. Two views, per the design of the 05 / 05b experiment:
#
#   FULL TABLE   every metric, every run x scenario x model, in one tidy file.
#   FORESIGHT    the rolling-vs-myopic gap — the quantity 05b was built to isolate.
#                Under a UNIFORM rate, risk carries no cross-parcel targeting signal
#                (survival S[i,t] = S[t] for all i), so whatever rolling-vs-myopic
#                gap survives is attributable to BENEFIT foresight alone. Sweeping
#                the rate maps how that gap scales with conversion magnitude, and
#                the data-derived run is the heterogeneous-risk reference it brackets.
#
# The five runs (folder = output_data_<run>):
#   data_derived   script 05  (EAU- and scenario-varying habitat-erosion hazard)
#   wetland_low    0.0009/yr  ┐
#   grass_low      0.0040/yr  │  script 05b uniform sweep, ascending magnitude
#   wetland_high   0.0057/yr  │
#   grass_high     0.0130/yr  ┘
#
# INPUTS:  output_data_<run>/model_results.csv   (one per run; from 08)
# OUTPUTS: output_data_comparison/
#            compare_run_manifest.csv         what loaded, with rate + p_decade
#            compare_all_metrics_long.csv      FULL tidy table (every row, every run)
#            compare_foresight_gap_long.csv    rolling/myopic/greedy J + gaps, per run x scenario
#            compare_foresight_gap_matrix.csv  HEADLINE: myopic value-added gap %, scenario x run
#            compare_value_added_matrix.csv    value_added, (scenario x model) x run
#
# Read-only on the run folders; writes only under output_data_comparison/.
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# ══ CONFIG ════════════════════════════════════════════════════════════════════####
# Where the output_data_<run> folders live (project root, same as 08/10 assume).
PROJECT_ROOT <- "."
OUT_DIR      <- file.path(PROJECT_ROOT, "output_data_comparison")

# Run registry. Order = reading order in every wide table: the data-derived
# primary analysis first, then the uniform sweep ascending in annual rate so the
# magnitude-sensitivity reads left-to-right. annual_rate = NA for the (heterogeneous)
# data-derived run. p_decade is the actual uniform trans_prob, 1-(1-rate)^10.
RUNS <- tibble::tribble(
  ~run,           ~annual_rate, ~risk_kind,
  "data_derived",   NA_real_,   "data-derived (heterogeneous)",
  "wetland_low",    0.0009,     "uniform",
  "grass_low",      0.0040,     "uniform",
  "wetland_high",   0.0057,     "uniform",
  "grass_high",     0.0130,     "uniform"
) %>%
  mutate(p_decade = 1 - (1 - annual_rate)^10)

RUN_LEVELS  <- RUNS$run
SCEN_LEVELS <- c("rcp45_wet", "rcp45_dry", "rcp85_wet", "rcp85_dry", "stationary")
MODEL_LEVELS <- c("greedy", "myopic", "rolling")

# Columns 08 always writes; scenario/model/J/J_baseline are load-critical, the rest
# are carried through if present.
CRITICAL_COLS <- c("scenario", "model", "J", "J_baseline")
CARRY_COLS    <- c("value_added", "n_acquired", "total_spend",
                   "gap_J_pct", "gap_value_added_pct")

TOL <- 1e-4   # relative tolerance for the internal cross-checks
# ══════════════════════════════════════════════════════════════════════════════


# ── 1. Load each run's model_results.csv (skip missing, don't die) ────────────####

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

load_run <- function(run) {
  path <- file.path(PROJECT_ROOT, paste0("output_data_", run), "model_results.csv")
  if (!file.exists(path))
    return(list(ok = FALSE, reason = "model_results.csv not found", path = path, data = NULL))
  df <- suppressMessages(read_csv(path, show_col_types = FALSE))
  miss <- setdiff(CRITICAL_COLS, names(df))
  if (length(miss))
    return(list(ok = FALSE,
                reason = paste("missing column(s):", paste(miss, collapse = ", ")),
                path = path, data = NULL))
  df$run <- run
  list(ok = TRUE, reason = "", path = path, data = df)
}

loaded <- lapply(RUN_LEVELS, load_run)
names(loaded) <- RUN_LEVELS

# Manifest: one row per expected run, whether it loaded, and basic shape.
# `loaded` is position-aligned with RUNS (both come from RUN_LEVELS), so build the
# columns with vapply over the list rather than indexing it inside a data-masked
# mutate (which mis-resolves loaded[[run]]).
manifest <- RUNS
manifest$loaded   <- vapply(loaded, function(x) isTRUE(x$ok), logical(1))
manifest$reason   <- vapply(loaded, function(x) x$reason,     character(1))
manifest$n_rows   <- vapply(loaded, function(x) if (isTRUE(x$ok)) nrow(x$data)                       else NA_integer_, integer(1))
manifest$n_scen   <- vapply(loaded, function(x) if (isTRUE(x$ok)) dplyr::n_distinct(x$data$scenario) else NA_integer_, integer(1))
manifest$n_models <- vapply(loaded, function(x) if (isTRUE(x$ok)) dplyr::n_distinct(x$data$model)    else NA_integer_, integer(1))
manifest$path     <- vapply(loaded, function(x) x$path,       character(1))

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  RISK-MODEL COMPARISON — run manifest\n")
cat("══════════════════════════════════════════════════════════════════\n\n")
print(as.data.frame(manifest %>%
                      mutate(annual_rate = ifelse(is.na(annual_rate), NA, round(annual_rate, 4)),
                             p_decade     = round(p_decade, 5)) %>%
                      select(run, risk_kind, annual_rate, p_decade, loaded, n_rows, n_scen, n_models, reason)),
      row.names = FALSE)

ok_runs <- manifest$run[manifest$loaded]
if (length(ok_runs) == 0) stop("No run folders loaded. Check PROJECT_ROOT and folder names.")
if (length(ok_runs) < nrow(RUNS))
  warning("Only ", length(ok_runs), " of ", nrow(RUNS),
          " runs loaded; comparison proceeds on those. See manifest 'reason'.")


# ── 2. Assemble the FULL tidy table ───────────────────────────────────────────####

all_long <- bind_rows(lapply(loaded[ok_runs], `[[`, "data")) %>%
  # keep critical + whatever optional columns are present
  select(run, any_of(c(CRITICAL_COLS, CARRY_COLS))) %>%
  left_join(RUNS, by = "run") %>%
  mutate(
    run      = factor(run,      levels = RUN_LEVELS),
    scenario = factor(scenario, levels = SCEN_LEVELS),
    model    = factor(model,    levels = MODEL_LEVELS),
    # recompute value_added from J so the table is self-consistent even if the
    # stored column is absent; keep the stored one (if any) for the cross-check below
    value_added_recomputed = J - J_baseline
  ) %>%
  arrange(run, scenario, model)


# ── 3. Internal cross-checks (report, don't stop) ─────────────────────────────####

cat("\n────────────────────────────────────────────────────────────────\n")
cat("  Internal consistency checks\n")
cat("────────────────────────────────────────────────────────────────\n")

flags <- list()

# 3a. Expected 5 scenarios x 3 models = 15 rows per loaded run
row_shape <- all_long %>% count(run, name = "n") %>% filter(n != 15L)
flags[["rows_per_run == 15"]] <- nrow(row_shape) == 0

# 3b. value_added matches J - J_baseline (if the stored column exists)
if ("value_added" %in% names(all_long)) {
  va_mismatch <- all_long %>%
    filter(abs(value_added - value_added_recomputed) >
             TOL * pmax(1, abs(value_added_recomputed)))
  flags[["value_added == J - J_baseline"]] <- nrow(va_mismatch) == 0
}

# 3c. J_baseline identical across models within run x scenario (08 assigns one ref)
base_var <- all_long %>%
  group_by(run, scenario) %>%
  summarise(n_base = dplyr::n_distinct(round(J_baseline, 6)), .groups = "drop") %>%
  filter(n_base > 1)
flags[["J_baseline constant across models"]] <- nrow(base_var) == 0

# 3d. Policy ordering rolling >= myopic >= greedy within run x scenario
ord <- all_long %>%
  select(run, scenario, model, J) %>%
  pivot_wider(names_from = model, values_from = J)
ord_bad <- ord %>%
  filter(rolling < myopic - TOL * pmax(1, abs(rolling)) |
           myopic  < greedy - TOL * pmax(1, abs(myopic)))
flags[["ordering rolling>=myopic>=greedy"]] <- nrow(ord_bad) == 0

for (nm in names(flags))
  cat(sprintf("  %s  %s\n", if (isTRUE(flags[[nm]])) "[PASS]" else "[WARN]", nm))

if (!isTRUE(flags[["ordering rolling>=myopic>=greedy"]])) {
  cat("\n  Ordering violations (rolling should dominate; small ones can be solver ties):\n")
  print(as.data.frame(ord_bad), row.names = FALSE)
}


# ── 4. FORESIGHT view: rolling vs myopic (and greedy as control) ──────────────####
# Rebuilt directly from J so it is independent of 08's stored gap columns:
#   foresight (myopic) gap, value-added terms = 100 * (J_r - J_m) / (J_r - J_b)
#   foresight (myopic) gap, total-J terms     = 100 * (J_r - J_m) / J_r
# The value-added gap is where model choice is visible; the total-J gap is the
# landscape headline (small, because acquisition touches a tiny share of EAUs).

pct <- function(num, den) ifelse(abs(den) < 1e-9, NA_real_, 100 * num / den)

wideJ <- all_long %>%
  select(run, annual_rate, p_decade, scenario, J_baseline, model, J) %>%
  # J_baseline is constant across models within run x scenario (checked in step 3),
  # so it stays an id column here and comes through correctly — no reassignment needed.
  pivot_wider(names_from = model, values_from = J,
              names_glue = "J_{model}")

foresight_long <- wideJ %>%
  transmute(
    run, annual_rate, p_decade, scenario,
    J_rolling, J_myopic, J_greedy, J_baseline,
    va_rolling = J_rolling - J_baseline,
    va_myopic  = J_myopic  - J_baseline,
    va_greedy  = J_greedy  - J_baseline,
    # myopic = the true foresight comparison (same information, differ only in horizon)
    myopic_gap_va_pct = pct(J_rolling - J_myopic, J_rolling - J_baseline),
    myopic_gap_J_pct  = pct(J_rolling - J_myopic, J_rolling),
    # greedy = control (its schedule is fixed across the sweep; only evaluated J moves)
    greedy_gap_va_pct = pct(J_rolling - J_greedy, J_rolling - J_baseline),
    greedy_gap_J_pct  = pct(J_rolling - J_greedy, J_rolling)
  ) %>%
  arrange(run, scenario)

# Cross-check the recomputed myopic VA gap against 08's stored gap_value_added_pct
if ("gap_value_added_pct" %in% names(all_long)) {
  stored_gap <- all_long %>%
    filter(model == "myopic") %>%
    select(run, scenario, stored = gap_value_added_pct)
  chk <- foresight_long %>%
    select(run, scenario, recomputed = myopic_gap_va_pct) %>%
    left_join(stored_gap, by = c("run", "scenario")) %>%
    filter(!is.na(stored) & !is.na(recomputed) &
             abs(stored - recomputed) > 1e-3 * pmax(1, abs(stored)))
  cat(sprintf("  %s  recomputed myopic VA gap == 08 stored gap_value_added_pct\n",
              if (nrow(chk) == 0) "[PASS]" else "[WARN]"))
  if (nrow(chk) > 0) print(as.data.frame(chk), row.names = FALSE)
}

# HEADLINE matrix: myopic value-added gap %, scenario (rows) x run (cols)
foresight_matrix <- foresight_long %>%
  select(scenario, run, myopic_gap_va_pct) %>%
  mutate(myopic_gap_va_pct = round(myopic_gap_va_pct, 3)) %>%
  pivot_wider(names_from = run, values_from = myopic_gap_va_pct) %>%
  arrange(scenario)


# ── 5. VALUE-ADDED matrix: (scenario x model) x run ───────────────────────────####

value_added_matrix <- all_long %>%
  select(run, scenario, model, value_added_recomputed) %>%
  mutate(value_added = round(value_added_recomputed, 1)) %>%
  select(-value_added_recomputed) %>%
  pivot_wider(names_from = run, values_from = value_added) %>%
  arrange(scenario, model)


# ── 6. Console read-out of the headline ───────────────────────────────────────####

cat("\n────────────────────────────────────────────────────────────────\n")
cat("  FORESIGHT gap (rolling vs myopic), value-added %  —  scenario x run\n")
cat("  Positive = rolling beats myopic. Higher = foresight matters more.\n")
cat("────────────────────────────────────────────────────────────────\n\n")
print(as.data.frame(foresight_matrix), row.names = FALSE)

# Does the foresight gap grow with conversion magnitude across the uniform sweep?
uni <- intersect(c("wetland_low", "grass_low", "wetland_high", "grass_high"), ok_runs)
if (length(uni) >= 2) {
  climate_mean <- foresight_long %>%
    filter(scenario != "stationary", run %in% uni) %>%
    group_by(run) %>%
    summarise(mean_myopic_gap_va_pct = round(mean(myopic_gap_va_pct, na.rm = TRUE), 3),
              .groups = "drop") %>%
    mutate(run = factor(run, levels = RUN_LEVELS)) %>%
    arrange(run)
  cat("\n  Climate-scenario mean foresight gap (VA %), uniform sweep ascending in rate:\n")
  print(as.data.frame(climate_mean), row.names = FALSE)
  mono <- all(diff(climate_mean$mean_myopic_gap_va_pct) >= -1e-6)
  cat(sprintf("    -> monotone increasing with rate: %s\n", if (mono) "yes" else "no"))
}
if ("data_derived" %in% ok_runs) {
  dd <- foresight_long %>%
    filter(run == "data_derived", scenario != "stationary")
  cat(sprintf("\n  data-derived run, climate scenarios: mean foresight gap (VA %%) = %.3f\n",
              mean(dd$myopic_gap_va_pct, na.rm = TRUE)))
}


# ── 7. Write CSVs ─────────────────────────────────────────────────────────────####

manifest_out <- manifest %>%
  select(run, risk_kind, annual_rate, p_decade, loaded, n_rows, n_scen, n_models, reason, path)

write_csv(manifest_out,        file.path(OUT_DIR, "compare_run_manifest.csv"))
write_csv(all_long %>% mutate(across(where(is.factor), as.character)),
          file.path(OUT_DIR, "compare_all_metrics_long.csv"))
write_csv(foresight_long %>% mutate(scenario = as.character(scenario)),
          file.path(OUT_DIR, "compare_foresight_gap_long.csv"))
write_csv(foresight_matrix %>% mutate(scenario = as.character(scenario)),
          file.path(OUT_DIR, "compare_foresight_gap_matrix.csv"))
write_csv(value_added_matrix %>% mutate(across(where(is.factor), as.character)),
          file.path(OUT_DIR, "compare_value_added_matrix.csv"))

cat("\n────────────────────────────────────────────────────────────────\n")
cat("  Wrote to", OUT_DIR, ":\n")
cat("    compare_run_manifest.csv\n")
cat("    compare_all_metrics_long.csv        (full table)\n")
cat("    compare_foresight_gap_long.csv      (foresight detail)\n")
cat("    compare_foresight_gap_matrix.csv    (headline: myopic VA gap, scenario x run)\n")
cat("    compare_value_added_matrix.csv\n")
cat("────────────────────────────────────────────────────────────────\n\n")