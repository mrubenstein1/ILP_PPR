# ══════════════════════════════════════════════════════════════════════════════
# test_C_spatial_alignment.R — Does prop_suitable sit on the right EAU?
# ══════════════════════════════════════════════════════════════════════════════

# STATUS: RESOLVED 2026-07-30. Sections 3-4 confirm the fix (100% assignment,
# exact membership). Sections 5-6 reconstruct the ORIGINAL positional mapping
# and will continue to report 0% / ~315 km by design — they measure the historic
# bug, not current state. The live check now lives at the end of 03.
#
# DIAGNOSTIC ONLY. Reads saved artifacts; writes nothing; changes nothing.
# Run from the project root (same working directory as 00_setup.R).
#
# THE QUESTION
#   Script 02 writes Lu_prop_*.csv in RASTER CELL ORDER (the order values() sweeps
#   the grid), keeping only cells with valid land cover.
#   Script 01 writes eau_wmd_lookup in ALPHABETICAL-WMD order, then re-indexes
#   eau_id = 1:N over that ordering.
#   Script 03 joins the two BY ROW POSITION, which is only correct if those two
#   orderings coincide. This script tests whether they do.
#
# THE BRIDGE
#   Script 02 saved Lu_prop_keep_rows_45.rds — a TRUE/FALSE flag per raster cell.
#   which(keep_rows) therefore recovers the RASTER CELL NUMBER for each row of
#   Lu_prop_45.csv: exactly the information script 03 discards.
#   cellFromXY() converts each EAU's stored centroid to its own cell number.
#   Joining on cell number = joining by PLACE instead of by ORDER.
#
# WHAT IT REPORTS
#   1. Geometry check    — is the cell numbering shared across the rasters?
#   2. Membership        — are the right EAUs in the panel at all?
#   3. Assignment        — does each EAU hold ITS OWN habitat value?  <- the point
#   4. Permutation test  — if not, is it a pure relabelling or real corruption?
#   5. Displacement      — how far did each value travel, in km?
# ══════════════════════════════════════════════════════════════════════════════

if (!exists(".SETUP_DONE")) source("00_setup.R")

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(terra)
})

TEST_YEAR <- 2030L   # any decade in 2030-2100; prop_suitable does not vary by gcm
TOL       <- 1e-9    # values are proportions; this is far below meaningful precision

hr <- function(txt) cat("\n", strrep("=", 78), "\n  ", txt, "\n",
                        strrep("=", 78), "\n", sep = "")
ok <- function(x) if (isTRUE(x)) "[PASS]" else "[FAIL]"


# ── 0. Locate inputs ──────────────────────────────────────────────────────────
hr("0. INPUTS")

P <- list(
  mask     = file.path(DIR_DERIVED, "foresce_eau_shared_mask.tif"),
  wmd_r    = file.path(DIR_DERIVED, "wmd_raster_equal_area.tif"),
  keep45   = file.path(DIR_DERIVED, "Lu_prop_keep_rows_45.rds"),
  lu45     = file.path(DIR_DERIVED, "Lu_prop_45.csv"),
  lookup   = file.path(DIR_DERIVED, "eau_wmd_lookup.csv"),
  panel    = file.path(DIR_DERIVED, "eau_panel.rds"),      # script 03 output
  panel_alt= file.path(DIR_DERIVED, "data_panel.rds")      # fallback (post-exclusion)
)

missing <- names(P)[!file.exists(unlist(P))]
missing <- setdiff(missing, c("panel", "panel_alt"))
if (length(missing))
  stop("Missing required file(s): ",
       paste(unlist(P)[missing], collapse = ", "))

# Prefer the script-03 panel: nothing has been dropped from it yet, so a
# membership failure below is unambiguously the row-order join and not the
# deliberate WMD exclusions in script 04.
use_full_panel <- file.exists(P$panel)
panel_path     <- if (use_full_panel) P$panel else P$panel_alt
if (!file.exists(panel_path))
  stop("Neither eau_panel.rds nor data_panel.rds found in ", DIR_DERIVED)

cat("  Panel under test : ", panel_path,
    if (use_full_panel) "  (pre-exclusion)" else "  (POST-exclusion — 4 WMDs already dropped)",
    "\n", sep = "")
cat("  Test year        : ", TEST_YEAR, " (RCP 4.5)\n", sep = "")


# ── 1. Geometry: do the rasters share a cell numbering? ───────────────────────
hr("1. GEOMETRY CHECK")

mask_r <- rast(P$mask)
wmd_r  <- rast(P$wmd_r)

geom_same <- isTRUE(all.equal(as.vector(ext(mask_r)), as.vector(ext(wmd_r)))) &&
             identical(dim(mask_r)[1:2], dim(wmd_r)[1:2])

cat(sprintf("  %s shared mask and EAU raster share extent + dimensions\n", ok(geom_same)))
cat(sprintf("       mask : %d rows x %d cols = %d cells\n",
            nrow(mask_r), ncol(mask_r), ncell(mask_r)))

if (!geom_same)
  stop("Rasters differ in geometry — cell numbers are not comparable. ",
       "Investigate before continuing; the rest of this test assumes a shared grid.")

keep45 <- readRDS(P$keep45)
len_ok <- length(keep45) == ncell(mask_r)
cat(sprintf("  %s keep_rows length (%d) equals raster cell count (%d)\n",
            ok(len_ok), length(keep45), ncell(mask_r)))
if (!len_ok)
  stop("keep_rows was computed against a different grid than the shared mask.")

lu45 <- read_csv(P$lu45, show_col_types = FALSE)
rows_ok <- nrow(lu45) == sum(keep45)
cat(sprintf("  %s Lu_prop_45.csv rows (%d) equals sum(keep_rows) (%d)\n",
            ok(rows_ok), nrow(lu45), sum(keep45)))
if (!rows_ok)
  stop("Row count mismatch — Lu_prop_45.csv was not produced by this keep_rows vector.")


# ── 2. Build the CORRECT cell -> eau_id key ───────────────────────────────────
hr("2. BUILDING THE CORRECT (BY-PLACE) KEY")

eau_wmd <- read_csv(P$lookup, show_col_types = FALSE)

eau_cells <- cellFromXY(mask_r, cbind(eau_wmd$x_coord, eau_wmd$y_coord))
if (anyNA(eau_cells))
  stop(sum(is.na(eau_cells)), " EAU centroid(s) fell outside the raster. ",
       "Coordinate/CRS problem — resolve before interpreting this test.")
if (anyDuplicated(eau_cells))
  stop("Two or more EAUs map to the same raster cell — lookup is not 1:1 with cells.")

eau_by_cell <- tibble(cell = eau_cells, eau_id = eau_wmd$eau_id)

# Row i of Lu_prop_45.csv describes raster cell which(keep45)[i].
col_name <- paste0("rcp45_", TEST_YEAR)
if (!col_name %in% names(lu45)) {
  avail <- grep("^rcp45_", names(lu45), value = TRUE)
  stop("Column '", col_name, "' not in Lu_prop_45.csv. Available: ",
       paste(avail, collapse = ", "), "\n  Edit TEST_YEAR at the top.")
}

lu_by_cell <- tibble(
  cell         = which(keep45),
  prop_correct = lu45[[col_name]]
)

correct <- inner_join(lu_by_cell, eau_by_cell, by = "cell") %>%
  select(eau_id, cell, prop_correct)

cat(sprintf("  Land-cover rows keyed to a cell : %d\n", nrow(lu_by_cell)))
cat(sprintf("  EAUs keyed to a cell            : %d\n", nrow(eau_by_cell)))
cat(sprintf("  EAUs that genuinely have land-cover data : %d\n", nrow(correct)))


# ── 3. Pull the panel's stored value for the same EAU-year ────────────────────
hr("3. MEMBERSHIP — are the right EAUs in the panel?")

panel <- readRDS(panel_path)

panel_slice <- panel %>%
  filter(year == TEST_YEAR, rcp == "45") %>%
  distinct(eau_id, prop_suitable) %>%
  rename(prop_panel = prop_suitable)

if (nrow(panel_slice) == 0)
  stop("No rows for year == ", TEST_YEAR, " & rcp == '45' in the panel.")
if (anyDuplicated(panel_slice$eau_id))
  stop("prop_suitable is not unique per EAU at this year/rcp — unexpected panel shape.")

in_panel_only   <- setdiff(panel_slice$eau_id, correct$eau_id)
in_correct_only <- setdiff(correct$eau_id, panel_slice$eau_id)

cat(sprintf("  EAUs in panel                              : %d\n", nrow(panel_slice)))
cat(sprintf("  EAUs that truly have land-cover data       : %d\n", nrow(correct)))
cat(sprintf("  In panel but WITHOUT real land-cover data  : %d\n", length(in_panel_only)))
cat(sprintf("  Has land-cover data but MISSING from panel : %d\n", length(in_correct_only)))

membership_ok <- length(in_panel_only) == 0 && length(in_correct_only) == 0
cat(sprintf("\n  %s panel membership matches true land-cover coverage\n", ok(membership_ok)))

if (!use_full_panel)
  cat("       (note: this panel has the 4 excluded WMDs already removed, so some\n",
      "        difference here is expected — read section 4 as the real verdict.)\n", sep = "")

if (length(in_panel_only)) {
  cat("\n  WMDs holding EAUs that were kept without real data:\n")
  eau_wmd %>% filter(eau_id %in% in_panel_only) %>%
    count(wmd_id, name = "n_eaus") %>% arrange(desc(n_eaus)) %>%
    as.data.frame() %>% print(row.names = FALSE)
}
if (length(in_correct_only)) {
  cat("\n  WMDs holding EAUs that were dropped despite having data:\n")
  eau_wmd %>% filter(eau_id %in% in_correct_only) %>%
    count(wmd_id, name = "n_eaus") %>% arrange(desc(n_eaus)) %>%
    as.data.frame() %>% print(row.names = FALSE)
}


# ── 4. ASSIGNMENT — does each EAU hold its OWN habitat value? ─────────────────
hr("4. ASSIGNMENT — THE DECISIVE TEST")

cmp <- inner_join(correct, panel_slice, by = "eau_id") %>%
  mutate(agrees = abs(prop_correct - prop_panel) <= TOL + 1e-6 * abs(prop_correct))

n_cmp   <- nrow(cmp)
n_agree <- sum(cmp$agrees)
pct     <- 100 * n_agree / n_cmp

cat(sprintf("  EAUs compared          : %d\n", n_cmp))
cat(sprintf("  Values that AGREE      : %d  (%.2f%%)\n", n_agree, pct))
cat(sprintf("  Values that DISAGREE   : %d  (%.2f%%)\n", n_cmp - n_agree, 100 - pct))

assignment_ok <- n_agree == n_cmp
cat(sprintf("\n  %s every EAU holds its own habitat value\n", ok(assignment_ok)))

if (!assignment_ok) {
  cat("\n  Worst mismatches (largest absolute difference):\n")
  cmp %>% filter(!agrees) %>%
    mutate(difference = prop_panel - prop_correct) %>%
    arrange(desc(abs(difference))) %>%
    left_join(eau_wmd %>% select(eau_id, wmd_id), by = "eau_id") %>%
    select(eau_id, wmd_id, prop_correct, prop_panel, difference) %>%
    head(10) %>% as.data.frame() %>% print(row.names = FALSE, digits = 5)
}


# ── 5. PERMUTATION TEST — relabelling, or genuine corruption? ─────────────────
hr("5. PERMUTATION TEST — what KIND of error is it?")

# If the two columns hold the SAME MULTISET of numbers but paired differently,
# the panel contains exactly the right values attached to the wrong parcels:
# a pure relabelling. Every value is still a real trajectory from a real parcel,
# which is what determines whether upstream computation can be trusted.
sorted_match <- isTRUE(all.equal(sort(cmp$prop_correct), sort(cmp$prop_panel),
                                 tolerance = 1e-8))

cat(sprintf("  %s the two columns hold the same set of values (sorted)\n",
            ok(sorted_match)))
cat(sprintf("       correct: mean %.6f | sd %.6f | min %.6f | max %.6f\n",
            mean(cmp$prop_correct), sd(cmp$prop_correct),
            min(cmp$prop_correct), max(cmp$prop_correct)))
cat(sprintf("       panel  : mean %.6f | sd %.6f | min %.6f | max %.6f\n",
            mean(cmp$prop_panel), sd(cmp$prop_panel),
            min(cmp$prop_panel), max(cmp$prop_panel)))

cat("\n  INTERPRETATION\n")
if (assignment_ok) {
  cat("    Values sit on the correct parcels. The positional join happens to be\n",
      "    valid here. No re-run needed on these grounds.\n", sep = "")
} else if (sorted_match) {
  cat("    PURE PERMUTATION. The panel holds exactly the right numbers on exactly\n",
      "    the wrong parcels. Each stored trajectory is a real, internally coherent\n",
      "    series from a real parcel — so distributions (including the trans_prob\n",
      "    distribution) are unaffected, and the low median hazard is a genuine\n",
      "    property of the FOREsce projections, not an artefact of this bug.\n",
      "    What IS wrong: habitat_share divides each parcel by a WMD sum assembled\n",
      "    from parcels outside that WMD, so scaled_abundance is genuinely\n",
      "    miscomputed, not merely relabelled.\n", sep = "")
} else {
  cat("    NOT a clean permutation — the value sets themselves differ. Something\n",
      "    beyond row order is involved (membership loss, a resample, or an issue\n",
      "    in script 02). Read section 3 before assuming the join is the whole story.\n",
      sep = "")
}


# ── 6. DISPLACEMENT — how far did each value travel? ──────────────────────────
hr("6. DISPLACEMENT — how far off is each value?")

# Reconstruct script 03's assumed mapping directly:
#   eau_row i  ->  eau_wmd row i        (its slice(seq_len(n)) + row_number())
#   eau_row i  ->  raster cell which(keep45)[i]   (the truth)
n_lu <- nrow(lu45)
if (n_lu <= nrow(eau_wmd)) {

  assumed_eau <- eau_wmd$eau_id[seq_len(n_lu)]          # what script 03 assigns
  true_cell   <- which(keep45)                           # where the data is from
  true_eau    <- eau_by_cell$eau_id[match(true_cell, eau_by_cell$cell)]

  perm <- tibble(assumed_eau_id = assumed_eau, true_eau_id = true_eau) %>%
    filter(!is.na(true_eau_id))

  coords <- eau_wmd %>% select(eau_id, x_coord, y_coord, wmd_id)

  perm <- perm %>%
    left_join(coords, by = c("assumed_eau_id" = "eau_id")) %>%
    rename(x_a = x_coord, y_a = y_coord, wmd_a = wmd_id) %>%
    left_join(coords, by = c("true_eau_id" = "eau_id")) %>%
    rename(x_t = x_coord, y_t = y_coord, wmd_t = wmd_id) %>%
    mutate(dist_km = sqrt((x_a - x_t)^2 + (y_a - y_t)^2) / 1000,
           same_eau = assumed_eau_id == true_eau_id,
           same_wmd = wmd_a == wmd_t)

  cat(sprintf("  Rows traced                     : %d\n", nrow(perm)))
  cat(sprintf("  Landed on the correct EAU       : %d (%.2f%%)\n",
              sum(perm$same_eau), 100 * mean(perm$same_eau)))
  cat(sprintf("  Landed in the correct WMD       : %d (%.2f%%)\n",
              sum(perm$same_wmd, na.rm = TRUE), 100 * mean(perm$same_wmd, na.rm = TRUE)))
  cat(sprintf("  Displacement (km): median %.1f | mean %.1f | 90th pct %.1f | max %.1f\n",
              median(perm$dist_km, na.rm = TRUE), mean(perm$dist_km, na.rm = TRUE),
              quantile(perm$dist_km, 0.90, na.rm = TRUE), max(perm$dist_km, na.rm = TRUE)))
  cat("\n  A median displacement of tens-to-hundreds of km means habitat values are\n",
      "  being drawn from a different part of the region entirely — which is what\n",
      "  breaks the natural correlation between habitat, hazard, and land cost.\n", sep = "")
} else {
  cat("  Skipped: Lu_prop rows exceed lookup rows; the assumed mapping cannot be\n",
      "  reconstructed as script 03 builds it.\n", sep = "")
}


# ── 7. VERDICT ────────────────────────────────────────────────────────────────
hr("7. VERDICT")

cat(sprintf("  %s membership  — the right EAUs are present\n", ok(membership_ok)))
cat(sprintf("  %s assignment  — each EAU holds its own value\n", ok(assignment_ok)))

if (assignment_ok && membership_ok) {
  cat("\n  No misalignment detected. Switching script 03 to a cell-number join is\n",
      "  still worth doing as a robustness measure, but no re-run is required.\n", sep = "")
} else {
  cat("\n  MISALIGNMENT CONFIRMED.\n\n",
      "  Fix (both pieces already exist on disk):\n",
      "    1. Script 01 — keep the original raster cell number as its own column\n",
      "       before eau_id is overwritten with the sequential index, and save it\n",
      "       into eau_wmd_lookup.\n",
      "    2. Script 03 — build the key from which(keep_rows_45) and join on cell\n",
      "       number instead of row position. Membership then handles itself: the\n",
      "       slice(seq_len(n_eau_panel)) and its warning both come out.\n",
      "    3. Promote section 4 of this script into a permanent logic check at the\n",
      "       end of 03, in the same PASS/FAIL style as the rest of the pipeline.\n\n",
      "  Scope of the re-run: scripts 03-06, then 08 and 10, then the results prose.\n",
      "  Cost is unaffected (script 06 georeferences via cellFromXY throughout).\n",
      sep = "")
}

cat("\n")
