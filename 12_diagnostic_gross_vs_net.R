# ══════════════════════════════════════════════════════════════════════════════
# 12_diagnostic_gross_vs_net.R  —  Does trans_prob measure what it should?
# ══════════════════════════════════════════════════════════════════════════════
#
# STANDALONE, READ-ONLY DIAGNOSTIC. Reads from input_data/ and derived_data/ but
# writes EXCLUSIVELY to its own folder, diag_gross_v_net/. Nothing this script
# produces can be mistaken for canonical derived or output data, and nothing the
# pipeline reads is touched. Safe to run at any time; safe to delete the whole
# output folder afterwards.
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
#   D. ANOMALY    — probe for the undiagnosed 2050 dip in per-period risk.
#
# ── INPUTS ────────────────────────────────────────────────────────────────────
#   input_data/prairie_potholes_gcam_ref_rcp45/*.tif
#   input_data/prairie_potholes_gcam_ref_rcp85/*.tif
#   input_data/prairie_potholes_gcam_ref_meanclim/..._meanclim_2020.tif
#   derived_data/foresce_eau_shared_mask.tif      (script 02)
#   derived_data/eau_wmd_lookup.csv               (script 01)
#   derived_data/data_panel.rds                   (script 06)  [for section A only]
#
# ── OUTPUTS ── all under diag_gross_v_net/ ────────────────────────────────────
#   diag_gross_v_net/diag_gross_vs_net.csv       long: eau_id x year x rcp
#   diag_gross_v_net/diag_summary.csv            the per-decade summary table
#   diag_gross_v_net/diag_gross_vs_net_scatter.png   (DRAW_FIGS)
#
# ── RUNTIME ───────────────────────────────────────────────────────────────────
# MEASURED baseline: ~1,230 s per year-pair for the naive 12-pass implementation,
# serial. Do NOT calibrate against script 02's "~4 min per scenario" comment --
# that is parallel wall-clock across all cores, and script 02 does far less work
# per unit.
#
# This version cuts each pair to 3 native passes (see the worker in section 1),
# taking temp I/O from ~56 GB to ~12 GB per pair, and dispatches pairs in parallel.
# Expect 2-3x from the I/O reduction on top of the parallel speedup: roughly
# 45-75 min at N_CORES = 4, against ~5.5 h for the naive version run serially.
# The run reports achieved speedup at the end, so the next estimate can be
# empirical rather than inferred.
#
# QUICK <- TRUE runs a 2-pair smoke test instead; it is OFF by default here
# because QUICK truncates to the 2020->2030 and 2030->2040 steps of RCP 4.5 only,
# which biases sections B and C toward the high-conversion head of the series and
# makes section D structurally unable to see the 2050 step.
# ══════════════════════════════════════════════════════════════════════════════

if (!exists(".SETUP_DONE")) source("00_setup.R")

if (!HAS_SPATIAL)
  stop("This diagnostic needs terra + sf. Install them, or run it on a machine ",
       "that has the spatial stack.")

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(tidyr); library(readr); library(stringr)
})


# ══ PARAMETERS ══════════════════════════════════════════════════════════════####

# Smoke test: process only the first two year-pairs of RCP 4.5 (2020->2030 and
# 2030->2040). Confirms the machinery and the section-A fidelity checks in ~2 min.
# Sections B/C/D are NOT interpretable under QUICK -- all 879 EAUs are present, but
# only 2 of 8 non-terminal periods and 1 of 2 RCPs, and the two retained periods are
# the highest-conversion ones. Outputs are suffixed "_QUICK" so a smoke-test artefact
# can never be mistaken for a real run.
QUICK <- FALSE

# Suitable-habitat classes. MUST match script 02 exactly, or section A will fail.
#   19 = Perennial Grass, 20 = Open Water, 26 = Grassland,
#   28 = Woody Wetland,   29 = Herbaceous Wetland
SUITABLE_CLASSES <- c(19, 20, 26, 28, 29)

# Tolerance for the fidelity checks in section A. prop_suitable is a ratio of
# pixel counts, so agreement should be near-exact; 1e-6 leaves room for float
# accumulation across ~3e5 pixels per EAU without masking a real discrepancy.
FIDELITY_TOL <- 1e-6

# ── Parallelism ───────────────────────────────────────────────────────────────
# Pairs are independent, so they are dispatched concurrently (script 02 does the
# same over files). Set to 1L to force the serial path.
#
# This job is DISK-bound, not CPU-bound: each pair writes several native-resolution
# intermediates to tempdir() and reads them back. Cores therefore contend for I/O
# and scaling is sublinear -- 7 cores will not give 7x. Start at 4; the run reports
# its achieved speedup at the end, so raise or lower it on evidence next time.
#
# Two things to check before raising it: free space in tempdir() (several GB per
# concurrent worker), and that a worker failure is loud -- it is, see section 3.
# ── Parallel backend ──────────────────────────────────────────────────────────
#   "psock"  — independent R processes. SAFE IN RSTUDIO. The default.
#   "fork"   — parallel::mclapply. Lower startup overhead, but forking inside
#              RStudio on macOS can hang or crash the IDE. Rscript only.
#   "serial" — one pair at a time. Debugging, or when N_CORES is 1.
#
# PSOCK pays a few seconds per worker at startup to launch a fresh R session, and
# workers cannot print to the RStudio console -- both handled below (the cost is
# trivial against ~20 min per pair, and progress goes to a log file instead).
PARALLEL_BACKEND <- "psock"

N_CORES <- 4L

# ── Scratch space ─────────────────────────────────────────────────────────────
# terra spills native-resolution intermediates to disk. By default that is R's
# tempdir(), which on macOS is under /var/folders/... on the BOOT volume. This run
# can ask for tens of GB across concurrent workers, so if the boot volume is tight,
# point this at a roomier disk (e.g. "/Volumes/Scratch/terra_tmp"). NULL = leave
# terra's default alone. The directory is created if it does not exist.
TERRA_TEMPDIR <- NULL

# ── Memory per worker ─────────────────────────────────────────────────────────
# terra's memfrac is the fraction of system RAM ONE terra process will use before
# chunking to disk, and forked workers each inherit it independently. Left at the
# default (~0.6), N concurrent workers can collectively commit N x 60% of RAM and
# drive the machine into swap -- which on macOS manifests as everything crawling
# rather than a clean error, so it is easy to misread as "the job is just slow".
# NULL = divide the default budget across workers automatically. Set a number to
# override (e.g. 0.1). Lower means more disk chunking and slightly slower workers,
# but you have plenty of scratch space and swapping is far more costly.
TERRA_MEMFRAC <- NULL

DRAW_FIGS <- TRUE     # write the gross-vs-net scatter to DIR_DIAG

# ── Output folder ─────────────────────────────────────────────────────────────
# Every artefact this script produces lands here and nowhere else. Deliberately
# NOT DIR_DERIVED or DIR_FIGS: diagnostic output must never sit alongside the
# canonical pipeline data, where a later reader could mistake it for a model input.
DIR_DIAG <- "diag_gross_v_net"
if (!dir.exists(DIR_DIAG)) dir.create(DIR_DIAG, recursive = TRUE)

# Filename tag. Empty on a full run; "_QUICK" on a smoke test, so partial output
# is self-identifying on disk and cannot silently overwrite a real run.
RUN_TAG <- if (QUICK) "_QUICK" else ""

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


# ── 0b. Scratch-space preflight ───────────────────────────────────────────####
#
# This run is disk-bound and the failure mode is ugly: workers die partway through
# with an out-of-space error, hours in. The requirement is computable from the grid
# rather than guessable, so compute it, compare against actual free space, and say
# so plainly BEFORE committing to a long run.
#
# Per pair the single-pass worker holds at most: the cropped 2-layer source stack
# and the 4-layer flow raster, all written as INT1U (1 byte/pixel, lossless for
# class codes and 0/1 flags). Allow 8 layer-equivalents for transients and
# compression scratch. This is ~5x smaller than the naive float implementation,
# which needed ~6 x 4-byte layers.

if (!is.null(TERRA_TEMPDIR)) {
  if (!dir.exists(TERRA_TEMPDIR)) dir.create(TERRA_TEMPDIR, recursive = TRUE)
  terra::terraOptions(tempdir = TERRA_TEMPDIR)
}

tmp_dir     <- terra::terraOptions(print = FALSE)$tempdir
native_px   <- as.numeric(n_cells) * scf
gb_per_rast <- native_px * 1 / 1024^3          # INT1U = 1 byte per pixel
gb_per_pair <- gb_per_rast * 8                 # 6 real layers + transients
gb_needed   <- gb_per_pair * max(1L, min(N_CORES, 16L))

# Free space on the volume holding tmp_dir. `df -Pk` is POSIX and works on macOS
# and Linux alike; if it is unavailable we simply skip the comparison rather than
# blocking the run on a diagnostic convenience.
gb_free <- tryCatch({
  df <- system2("df", c("-Pk", shQuote(tmp_dir)), stdout = TRUE, stderr = FALSE)
  as.numeric(strsplit(trimws(df[length(df)]), "\\s+")[[1]][4]) * 1024 / 1024^3
}, error = function(e) NA_real_, warning = function(w) NA_real_)

cat(sprintf("\n  Scratch dir      : %s\n", tmp_dir))
cat(sprintf("  Native pixels    : %.2e  (~%.1f GB per INT1U layer)\n",
            native_px, gb_per_rast))
cat(sprintf("  Peak per worker  : ~%.0f GB   x %d worker(s) = ~%.0f GB needed\n",
            gb_per_pair, min(N_CORES, 16L), gb_needed))
if (is.na(gb_free)) {
  cat("  Free space       : could not determine -- check manually before running.\n")
} else {
  cat(sprintf("  Free space       : ~%.0f GB\n", gb_free))
  if (gb_free < gb_needed)
    stop(sprintf(paste0("Not enough scratch space: ~%.0f GB free, ~%.0f GB needed ",
                        "at N_CORES = %d.\n  Options: lower N_CORES (need is linear ",
                        "in it, ~%.0f GB per worker), free space on that volume, or ",
                        "set TERRA_TEMPDIR to a roomier disk.\n  Stopping now rather ",
                        "than failing hours into the run."),
                 gb_free, gb_needed, N_CORES, gb_per_pair), call. = FALSE)
  if (gb_free < 1.5 * gb_needed)
    warning(sprintf("Scratch space is tight: ~%.0f GB free vs ~%.0f GB needed. ",
                    gb_free, gb_needed),
            "The estimate is an upper bound, so this may well be fine, but ",
            "consider lowering N_CORES.", call. = FALSE)
}


# ── 1. Per-year-pair flow decomposition ───────────────────────────────────####
#
# Everything happens at native FOREsce resolution BEFORE aggregation, which is the
# whole point -- script 02's resample-then-difference destroys the information that
# separates L from G.
#
#   st      = 1 where suitable at t,    0 otherwise,  NA outside coverage
#   st1     = 1 where suitable at t+1,  0 otherwise,  NA outside coverage
#   lost    = st * (1 - st1)     -> 1 exactly where suitable -> unsuitable
#   gained  = (1 - st) * st1     -> 1 exactly where unsuitable -> suitable
#
# All four are produced in a SINGLE pass over the native data (see the worker
# below) and summed to EAU resolution in a single multi-layer resample, then
# divided by scf. That gives the fraction of the EAU in each category -- the same
# units as prop_suitable, so the four are directly comparable and the stock
# identity ps_t1 = ps_t - lost + gained can be checked exactly.
#
# NOTE ON UNSUITABLE CELLS: script 02 maps non-suitable classes to NA and relies on
# resample-sum skipping them; we map them to 0 and let them contribute nothing to
# the sum. The counts of 1s are identical, so the two agree EXCEPT for an EAU with
# no suitable pixels at all, which yields NA under 02's encoding and 0 here.
# Section A detects any such case -- it is a candidate explanation for part of the
# 41-EAU coverage drop, and worth knowing about either way.

# ── Self-contained worker (PSOCK- and fork-safe) ──────────────────────────────
#
# The serial one-slot cache that used to live here has been REMOVED: it is
# inherently sequential and cannot coexist with parallel dispatch over pairs.
#
# CRITICAL for PSOCK: this function receives ONLY serialisable arguments -- file
# paths, a numeric vector, and numbers. It touches no object from the parent
# environment. terra objects are external C++ pointers that cannot cross a PSOCK
# connection at all (and are only accidentally safe across a fork), so the EAU mask
# is re-opened from its PATH inside the worker. Re-opening reads a header only.
#
# ── I/O BUDGET (this job is disk-bound, so passes over native data ARE the cost) ─
#
# The obvious implementation costs 12 full-raster passes per pair: 2 crop,
# 2 classify, 4 arithmetic (each complement and each product materialises its own
# native intermediate), 4 resample. At ~4 bytes/pixel once arithmetic promotes byte
# input to float, that is ~56 GB of temp I/O per pair.
#
# This version collapses it to 3 passes:
#   1. ONE crop over a 2-layer stack of the two source years.
#   2. ONE lapp over that stack emitting FOUR layers at once -- st, st1, lost,
#      gained. The complements and products never become separate rasters; they are
#      transient vectors inside the chunk callback. Written as INT1U, which is
#      lossless for 0/1 data and 4x smaller than the default float.
#   3. ONE resample of the resulting 4-layer stack (resample is multi-layer aware),
#      replacing four separate calls that each re-read the same geometry.
# ~56 GB -> ~12 GB per pair. Expect 2-3x wall-clock, not the 4.5x the I/O ratio
# implies, since CPU and compression overhead do not shrink proportionally.
#
# WHAT WAS DELIBERATELY NOT DONE: resample(method = "sum") is NOT replaced by
# aggregate(fact = ...), which would be far faster but would no longer reproduce
# script 02's aggregation exactly -- and check A2 (recomputed prop_suitable
# matching the panel) is the entire basis for comparing these results to the
# pipeline. ps_t1 is also still a real summed layer rather than being derived as
# ps_t - L + G, which would make check A1 true by construction.
#
# ON A1's STRENGTH: st and st1 now come from one lapp pass rather than two separate
# classify calls, so A1 is marginally less independent than before. It still tests
# that resample-sum is linear and consistent ACROSS LAYERS, which is the property
# that would break under a geometry or alignment fault, so it remains a real check
# -- but it is no longer an end-to-end test of two independently built rasters.

decompose_pair <- function(path_t, path_t1, mask_path, classes, scf, n_cells, label,
                           memfrac = NULL, tmpdir = NULL, progress_file = NULL) {
  t0 <- Sys.time()
  
  # Worker-local terra settings. Each PSOCK worker is a fresh R process with
  # default options, so these must be set HERE, not in the parent.
  if (!is.null(tmpdir))  terra::terraOptions(tempdir = tmpdir)
  if (!is.null(memfrac)) terra::terraOptions(memfrac = memfrac)
  
  aoi <- terra::rast(mask_path)      # rebuilt per worker, never inherited
  
  r_t  <- terra::rast(path_t)
  r_t1 <- terra::rast(path_t1)
  
  # Stacking requires identical geometry. FOREsce years share a grid, so this
  # should always hold -- but if it ever does not, failing here with both filenames
  # is far better than silently mis-registering loss against gain.
  if (!terra::compareGeom(r_t, r_t1, stopOnError = FALSE))
    stop("Source rasters do not share a geometry: ", basename(path_t), " vs ",
         basename(path_t1), ". They cannot be stacked for the single-pass ",
         "decomposition; the pair-level flows would be spatially mis-registered.")
  
  # PASS 1 -- one crop over both layers.
  s <- terra::crop(c(r_t, r_t1), aoi)
  
  # PASS 2 -- one arithmetic pass emitting all four native layers.
  # Returning a 4-column matrix from fun gives a 4-layer result.
  #
  # NA discipline: `%in%` yields FALSE (not NA) for NA input, so a no-coverage
  # pixel would silently become 0 and enter the sums as real unsuitable land. The
  # explicit na mask below restores NA across all four outputs, which is what keeps
  # the downstream `filter(!is.na(ps_t))` meaning "no coverage here" rather than
  # "covered, and none of it suitable".
  flow <- terra::lapp(
    s,
    fun = function(a, b) {
      st  <- as.integer(a %in% classes)
      st1 <- as.integer(b %in% classes)
      na  <- is.na(a) | is.na(b)
      lost   <- st * (1L - st1)      # suitable -> unsuitable
      gained <- (1L - st) * st1      # unsuitable -> suitable
      st[na] <- NA_integer_; st1[na]    <- NA_integer_
      lost[na] <- NA_integer_; gained[na] <- NA_integer_
      cbind(st, st1, lost, gained)
    },
    wopt = list(datatype = "INT1U")   # lossless for 0/1; NA flag is 255
  )
  names(flow) <- c("st", "st1", "lost", "gained")
  
  # PASS 3 -- one multi-layer resample, then scale and mask.
  props <- terra::mask(terra::resample(flow, aoi, method = "sum") / scf, aoi)
  
  v <- terra::values(props, na.rm = FALSE)   # n_cells x 4, layer order as named
  out <- data.frame(
    cell_id     = seq_len(n_cells),
    ps_t        = v[, 1],
    ps_t1       = v[, 2],
    prop_lost   = v[, 3],
    prop_gained = v[, 4]
  )
  
  # Values are extracted; drop the raster handles so terra can release the backing
  # temp files before the next pair starts. Concurrent workers would otherwise
  # leave tempdir() growing monotonically for the whole run.
  rm(v, props, flow, s, r_t, r_t1, aoi)
  invisible(gc(verbose = FALSE))
  
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  # PSOCK workers have no console: anything they cat() is invisible in RStudio,
  # so a long run would show no sign of life for hours. Append a line to a shared
  # log instead, which can be watched with `tail -f` or just opened periodically.
  if (!is.null(progress_file)) {
    try(cat(sprintf("[%s] done %-16s %7.1f s\n",
                    format(Sys.time(), "%H:%M:%S"), label, secs),
            file = progress_file, append = TRUE), silent = TRUE)
  }
  
  attr(out, "seconds") <- secs
  attr(out, "label")   <- label
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

n_pairs <- nrow(pair_tbl)
cat(sprintf("\n  Year-pairs to process: %d\n", n_pairs))


# ── 3. Run the decomposition (parallel over pairs) ────────────────────────####

hdr("Processing rasters (this is the slow part)")

n_use  <- max(1L, min(N_CORES, n_pairs))
waves  <- ceiling(n_pairs / n_use)

# memfrac is PER PROCESS, so divide the single-process budget across workers to
# stop N workers each sizing chunks against 60% of total RAM and swapping.
memfrac_use <- if (!is.null(TERRA_MEMFRAC)) TERRA_MEMFRAC else
  max(0.05, round(0.6 / n_use, 3))

progress_file <- file.path(DIR_DIAG, sprintf("progress%s.log", RUN_TAG))
cat(sprintf("run started %s | %d pairs | %d worker(s) | backend %s\n",
            format(Sys.time()), n_pairs, n_use, PARALLEL_BACKEND),
    file = progress_file)   # truncates any previous run's log

cat(sprintf("  Backend: %s | workers: %d of %d detected | pairs: %d | ~%d wave(s)\n",
            PARALLEL_BACKEND, n_use, parallel::detectCores(), n_pairs, waves))
cat(sprintf("  terra memfrac per worker: %.3f%s\n", memfrac_use,
            if (is.null(TERRA_MEMFRAC)) "  (auto: 0.6 / workers)" else "  (manual)"))
cat(sprintf("  Scratch: %s (checked at startup)\n", tmp_dir))
cat(sprintf("  Progress log: %s\n", progress_file))
if (identical(PARALLEL_BACKEND, "psock") && n_use > 1L)
  cat("  PSOCK workers are silent by design -- watch the progress log for per-pair\n  completions (`tail -f` in Terminal, or just reopen the file).\n")
cat("\n")

t_all <- Sys.time()

run_one <- function(k) {
  p <- pair_tbl[k, ]
  # try() so that a failed pair comes back as an ELEMENT rather than aborting the
  # whole dispatch. mclapply does this itself; parLapplyLB instead raises, which
  # would kill the run on the first failure and lose both the identity of the
  # failed pair and the 15 completed ones. Wrapping here makes both backends
  # behave identically and lets the check below name every pair that died.
  try(
    decompose_pair(
      path_t        = p$path_t,
      path_t1       = p$path_t1,
      mask_path     = SHARED_MASK,
      classes       = SUITABLE_CLASSES,
      scf           = scf,
      n_cells       = n_cells,
      label         = sprintf("rcp%s %d", p$rcp, p$year),
      memfrac       = memfrac_use,
      tmpdir        = TERRA_TEMPDIR,
      progress_file = progress_file
    ),
    silent = TRUE
  )
}

if (n_use == 1L || identical(PARALLEL_BACKEND, "serial")) {
  flows <- lapply(seq_len(n_pairs), run_one)
  
} else if (identical(PARALLEL_BACKEND, "psock")) {
  # PSOCK: independent R processes rather than forks. Slower to start (each worker
  # is a fresh R session) but safe inside RStudio on macOS, where forking can hang
  # or crash the IDE. Workers share no memory, so everything they need must be
  # exported explicitly -- which is only possible because decompose_pair() takes
  # nothing but paths, a matrix and numbers. A SpatRaster could not be sent at all.
  cl <- parallel::makeCluster(n_use, type = "PSOCK")
  
  # NB tryCatch/finally, not on.exit(): on.exit registers against a function frame,
  # and at top level in a sourced script there is none, so it would fire straight
  # away and stop the cluster before any work ran. finally guarantees teardown on
  # both success and error, so a failure part-way cannot leave orphan R processes.
  flows <- tryCatch({
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(terra)))
    parallel::clusterExport(
      cl,
      varlist = c("decompose_pair", "pair_tbl", "SHARED_MASK", "SUITABLE_CLASSES",
                  "scf", "n_cells", "memfrac_use", "TERRA_TEMPDIR",
                  "progress_file"),
      envir = environment()
    )
    # parLapplyLB load-balances: pairs are handed out as workers free up, which is
    # the PSOCK equivalent of mc.preschedule = FALSE and matters because the
    # meanclim-anchored pairs cost differently from the rest.
    parallel::parLapplyLB(cl, seq_len(n_pairs), run_one)
  }, finally = try(parallel::stopCluster(cl), silent = TRUE))
  
} else if (identical(PARALLEL_BACKEND, "fork")) {
  # Fork: lower overhead, but UNSAFE inside RStudio on macOS. Only for Rscript.
  flows <- parallel::mclapply(seq_len(n_pairs), run_one,
                              mc.cores = n_use, mc.preschedule = FALSE)
  
} else {
  stop("PARALLEL_BACKEND must be one of 'psock', 'fork', 'serial'.")
}

# Both backends can report a worker failure as a try-error element in the returned
# list rather than raising (mclapply always does; parLapplyLB does for some
# failures). Silently binding those would drop whole year-pairs and quietly shrink
# the panel. Fail loudly instead.
bad <- which(!vapply(flows, is.data.frame, logical(1)))
if (length(bad)) {
  cat("\n  Failed pairs:\n")
  for (k in bad) {
    cond <- attr(flows[[k]], "condition")
    cat(sprintf("    %s: %s\n",
                sprintf("rcp%s %d", pair_tbl$rcp[k], pair_tbl$year[k]),
                if (is.null(cond)) paste(utils::capture.output(print(flows[[k]])),
                                         collapse = " ") else conditionMessage(cond)))
  }
  stop(length(bad), " of ", n_pairs, " pairs failed. If these are memory or disk ",
       "errors, lower N_CORES and re-run; results from a partial set would be ",
       "silently incomplete.")
}

cat("  Per-pair timings:\n")
for (k in seq_len(n_pairs))
  cat(sprintf("    %-28s %7.1f s\n",
              attr(flows[[k]], "label"), attr(flows[[k]], "seconds")))

for (k in seq_len(n_pairs)) {
  flows[[k]]$rcp  <- pair_tbl$rcp[k]
  flows[[k]]$year <- pair_tbl$year[k]
}

cpu_s  <- sum(vapply(flows, function(d) attr(d, "seconds"), numeric(1)))
wall_s <- as.numeric(difftime(Sys.time(), t_all, units = "secs"))

flows <- bind_rows(flows) %>%
  inner_join(cell_to_eau, by = "cell_id")          # by-place join, as in script 03

# Drop rows lacking ANY flow component, not just ps_t. Filtering on ps_t alone
# lets through rows where ps_t is present but ps_t1 / prop_lost / prop_gained are
# NA, which then propagate silently into gross_loss and net_loss and turn every
# downstream sum() into NA. Report the two categories separately: rows missing
# ps_t are the expected no-coverage EAUs, whereas rows missing only a partner
# column are unexpected and worth knowing about.
n_pre        <- nrow(flows)
miss_ps_t    <- sum(is.na(flows$ps_t))
miss_partner <- sum(!is.na(flows$ps_t) &
                      (is.na(flows$ps_t1) | is.na(flows$prop_lost) |
                         is.na(flows$prop_gained)))

flows <- flows %>%
  filter(!is.na(ps_t), !is.na(ps_t1), !is.na(prop_lost), !is.na(prop_gained))

cat(sprintf("\n  Rows: %d joined -> %d retained\n", n_pre, nrow(flows)))
cat(sprintf("    dropped, no coverage (ps_t NA)      : %d  (expect ~41 per pair)\n",
            miss_ps_t))
cat(sprintf("    dropped, partial decomposition      : %d%s\n", miss_partner,
            if (miss_partner > 0) "  <-- UNEXPECTED, see note below" else ""))
if (miss_partner > 0)
  warning(miss_partner, " row(s) had ps_t but were missing another flow component. ",
          "These are excluded. If the count is more than a handful of edge EAUs, ",
          "reconcile it before citing sections B-D.", call. = FALSE)

cat(sprintf("\n  Wall clock: %.1f min | summed worker time: %.1f min | speedup %.1fx\n",
            wall_s / 60, cpu_s / 60, cpu_s / max(wall_s, 1e-9)))
cat(sprintf("  EAU-year rows: %d | distinct EAUs: %d\n",
            nrow(flows), n_distinct(flows$eau_id)))
if (n_use > 1L && cpu_s / max(wall_s, 1e-9) < 0.6 * n_use)
  cat(sprintf("  NOTE: speedup is well under %dx. This job is disk-bound, so extra\n        cores contend for I/O rather than adding throughput.\n", n_use))


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
#
# NOTE: this deliberately reports the NA count as well as the max residual. An
# earlier version used max(..., na.rm = TRUE) alone, which meant any row with a
# missing flow component was skipped by the very check meant to catch it -- A1
# reported PASS while broken rows sat in the data and blew up section B. With the
# filter above, n_resid_na should always be 0; if it is not, A1 FAILS.
n_resid_na <- sum(is.na(metrics$ident_resid))
max_resid  <- max(abs(metrics$ident_resid), na.rm = TRUE)
a1 <- is.finite(max_resid) && max_resid < FIDELITY_TOL && n_resid_na == 0
cat(sprintf("  [%s] stock identity  ps_t1 == ps_t - L + G   (max resid %.2e, NA rows %d)\n",
            ifelse(a1, "PASS", "FAIL"), max_resid, n_resid_na))

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
    # gross_minus_net is always defined and is the quantity of interest: how much
    # conversion the net metric drops in this decade, in loss-rate units.
    gross_minus_net = gross_mean - net_mean,
    # The ratio is only meaningful with a strictly positive denominator. Where mean
    # net loss is zero or negative -- restoration outweighing conversion on average,
    # which is exactly the regime this diagnostic exists to expose -- the ratio is
    # undefined or sign-flipped, so report NA rather than a guarded pseudo-value.
    ratio        = ifelse(net_mean > 0, gross_mean / net_mean, NA_real_),
    pct_censored = 100 * mean(censored),
    .groups = "drop"
  )

print(by_decade %>%
        mutate(across(c(gross_mean, gross_median, gross_max, gain_mean,
                        net_mean, net_median, gross_minus_net), ~ round(.x, 5)),
               ratio = round(ratio, 2), pct_censored = round(pct_censored, 1)) %>%
        as.data.frame(), row.names = FALSE)

if (any(is.na(by_decade$ratio)))
  cat(sprintf("\n  NOTE: ratio is NA in %d decade-RCP cell(s) where mean net loss <= 0.\n         Read gross_minus_net for those; the ratio has no finite value there.\n",
              sum(is.na(by_decade$ratio))))

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


# ══ D. THE 2050 ANOMALY ══════════════════════════════════════════════════════####
hdr("D. ANOMALY PROBE -- the undiagnosed 2050 dip")

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

f_long <- file.path(DIR_DIAG, sprintf("diag_gross_vs_net%s.csv", RUN_TAG))
f_summ <- file.path(DIR_DIAG, sprintf("diag_summary%s.csv",      RUN_TAG))

write_csv(out_long,  f_long)
write_csv(by_decade, f_summ)
cat(sprintf("  %s  (%d rows)\n", f_long, nrow(out_long)))
cat(sprintf("  %s\n", f_summ))


# ── Figures (optional) ────────────────────────────────────────────────────####
if (DRAW_FIGS && requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  
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
         title = "Every point above the dashed line understates conversion risk",
         subtitle = "Points left of the vertical line are floored to epsilon despite real habitat loss",
         # gross - net = G/ps_t >= 0 identically, so NO point can fall below the
         # diagonal. One that did would mean negative gain -- i.e. a broken
         # decomposition -- so the dashed line doubles as a visual check on A1.
         caption = "Gross - net = gain/ps_t >= 0, so points below the diagonal are impossible; the line is also a visual check on the decomposition.") +
    theme_minimal(base_size = 11) + theme(legend.position = "top")
  
  f_fig <- file.path(DIR_DIAG, sprintf("diag_gross_vs_net_scatter%s.png", RUN_TAG))
  ggsave(f_fig, p1, width = 9, height = 4.8, dpi = 200, bg = "white")
  cat(sprintf("  %s\n", f_fig))
}

hdr("Done -- no pipeline file was modified")
cat(sprintf("  All output written to %s/ only.\n", DIR_DIAG))
if (QUICK) cat("  *** QUICK mode: partial coverage, files tagged _QUICK.",
               "Set QUICK <- FALSE for the real numbers. ***\n")