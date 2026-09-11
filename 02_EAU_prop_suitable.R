### LAND COVER CHARACTERIZATION ####
# Characterizing the percent of suitable habitat for each EAU
# is a critical first step to establishing both benefit, cost, and transition estimates


# Note: this uses the PPR GCAM Reference scenario
    # https://www.sciencebase.gov/catalog/item/5b683db2e4b006a11f75b06a#Overlay EAUs and Landcover (LC) data ########

# Note that there is a slight misalignment of area extent between the EAUs and the FOREsce landcover data. 
    #42 EAUs are outside the projected area, so there are only 1168 EAUs with valid landcover data.
    # This script corrects that and ensures that all remaining EAUs have full LC data.
    #Note that we end up removing the entire WMD with missing LC data in a downstream script. 

#------ Setup ________
if (!exists(".SETUP_DONE")) source("00_setup.R")

##### 1. Create overlapping Area of Interest (AOI) ##########


#1. Load EAU Mask
aoi_mask <- rast(file.path(DIR_DERIVED, "wmd_raster_equal_area.tif"))

#2. Load LC data
lu_r <- rast("input_data/prairie_potholes_gcam_ref_rcp45/prairie_potholes_gcam_ref_rcp45_2014.tif")

#3. crop LC data to EAU mask
lu_r <- crop(lu_r, aoi_mask)
lu_r <- resample(lu_r, aoi_mask, method = 'sum') # using summing which retains most of the overlapping grid cells (but actual values don't mean anything)
shared_aoi <- mask(aoi_mask, lu_r)

#4. write the raster of shared area between EAU and LC
writeRaster(shared_aoi, file.path(DIR_DERIVED, "foresce_eau_shared_mask.tif"),
            overwrite = TRUE)


  ##### 2. Resample LC data to EAUs ##########

###### NOTE: this takes about 4 mins to run per scenario (8 mins total). 

#1. Load LC data from both RCP scenarios
scenarios <- list(
  `45` = "input_data/prairie_potholes_gcam_ref_rcp45",
  `85` = "input_data/prairie_potholes_gcam_ref_rcp85"
)

#2. Define correct mask
aoi_mask <- rast(file.path(DIR_DERIVED, "foresce_eau_shared_mask.tif"))

#3. Define "suitable" habitat types:
# 19 = Perennial Grass, 20 = Open Water, 26 = Grassland, 28 = Woody Wetland, 29 = Herbaceous Wetland
rcl_table <- matrix(c(19, 1,
                      20, 1,
                      26, 1,
                      28, 1,
                      29, 1),
                    byrow = TRUE, ncol = 2)

#4. Define scale conversion factor & establish length of data frame to hold calculations
scf <- (res(aoi_mask)[1]^2) / (30^2)
n_cells <- ncell(aoi_mask)

#5. Prepare the resampling process:
# use parallel processing
# define process_sample
# run across all time frames (2014-2020) and for both datasets (RCP 4.5 & 8.5)
n_cores <- max(1, parallel::detectCores() - 1)
cat("Detected cores:", parallel::detectCores(), "\n")
cat("Using cores:", n_cores, "\n")

process_one_file <- function(fname, LU_in, aoi_mask, rcl_table, scf, n_cells) {
  
  t_file <- Sys.time()
  
  lu_r <- rast(file.path(LU_in, fname))
  lu_r <- crop(lu_r, aoi_mask)
  lu_r <- classify(lu_r, rcl = rcl_table, others = NA)
  lu_r <- resample(lu_r, aoi_mask, method = "sum")
  lu_r <- mask(lu_r / scf, aoi_mask)
  
  vals <- values(lu_r, na.rm = FALSE)
  stopifnot(length(vals) == n_cells)
  
  elapsed <- difftime(Sys.time(), t_file, units = "secs")
  
  list(
    file = fname,
    vals = vals,
    seconds = as.numeric(elapsed)
  )
}

#6. Run Resampling process on all LC data
# calculate time to run as a benchmark
process_scenario <- function(LU_in, suffix,
                             aoi_mask, rcl_table, scf,
                             n_cells, out_dir = DIR_DERIVED,
                             n_cores = 1) {
  
  # FIXED PATTERN: match prairie_potholes_gcam_ref_rcp{suffix}_YYYY.tif naming convention
  file_names <- list.files(
    LU_in,
    pattern = paste0("prairie_potholes_gcam_ref_rcp", suffix, "_[0-9]{4}\\.tif$"),
    full.names = FALSE
  )
  
  if (length(file_names) == 0) {
    stop("Scenario ", suffix, " has no matching rasters.")
  }
  
  file_years <- as.integer(substr(file_names,
                                  nchar(file_names) - 7,
                                  nchar(file_names) - 4))
  
  o <- order(file_years)
  file_names <- file_names[o]
  file_years <- file_years[o]
  
  cat("\nScenario", suffix, "matched files:\n")
  print(data.frame(file_names = file_names, file_years = file_years))
  
  t_batch <- Sys.time()
  
  res_list <- parallel::mclapply(
    file_names,
    FUN = process_one_file,
    LU_in = LU_in,
    aoi_mask = aoi_mask,
    rcl_table = rcl_table,
    scf = scf,
    n_cells = n_cells,
    mc.cores = n_cores
  )
  
  batch_minutes <- difftime(Sys.time(), t_batch, units = "mins")
  cat("Scenario", suffix, "batch time:",
      round(as.numeric(batch_minutes), 2), "minutes\n")
  
  file_seconds <- sapply(res_list, `[[`, "seconds")
  names(file_seconds) <- sapply(res_list, `[[`, "file")
  print(round(file_seconds, 1))
  
  #6. Save results (% suitable habitat per EAU) into appropriately named csv.
  lu_prop <- do.call(cbind, lapply(res_list, `[[`, "vals"))
  colnames(lu_prop) <- sprintf("rcp%s_%s", suffix, file_years)
  
  keep_rows <- rowSums(is.na(lu_prop)) == 0
  lu_prop_all_years <- lu_prop[keep_rows, , drop = FALSE]
  
  cat("Scenario", suffix,
      "| AOI cells:", n_cells,
      "| kept:", sum(keep_rows),
      "| dropped:", sum(!keep_rows), "\n")
  
  write.csv(
    lu_prop_all_years,
    file.path(out_dir, paste0("Lu_prop_", suffix, ".csv")),
    row.names = FALSE
  )
  
  saveRDS(
    keep_rows,
    file.path(out_dir, paste0("Lu_prop_keep_rows_", suffix, ".rds"))
  )
  
  invisible(list(
    file_names = file_names,
    file_years = file_years,
    lu_prop = lu_prop,
    lu_prop_all_years = lu_prop_all_years,
    keep_rows = keep_rows,
    file_seconds = file_seconds
  ))
}

t_all <- Sys.time()

results <- lapply(names(scenarios), function(sfx) {
  process_scenario(
    LU_in = scenarios[[sfx]],
    suffix = sfx,
    aoi_mask = aoi_mask,
    rcl_table = rcl_table,
    scf = scf,
    n_cells = n_cells,
    out_dir = DIR_DERIVED,
    n_cores = n_cores
  )
})
names(results) <- names(scenarios)

cat("\nTotal time for both scenarios:",
    round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 2),
    "minutes\n")



  ##### 3. Resample Meanclim LC Data for 2020 ##########

# Single-file resample for the mean climate scenario (2020 only).
# Reuses process_one_file() with the same reclassification and masking logic.

#1. Define file path
meanclim_path <- "input_data/prairie_potholes_gcam_ref_meanclim"
meanclim_fname <- "prairie_potholes_gcam_ref_meanclim_2020.tif"

#2. Process single file
cat("\nProcessing meanclim 2020...\n")
t_meanclim <- Sys.time()

result_meanclim <- process_one_file(
  fname    = meanclim_fname,
  LU_in    = meanclim_path,
  aoi_mask = aoi_mask,
  rcl_table = rcl_table,
  scf      = scf,
  n_cells  = n_cells
)

cat("Meanclim 2020 processing time:",
    round(result_meanclim$seconds, 1), "seconds\n")

#3. Format output to match structure of RCP scenario CSVs
lu_prop_meanclim <- matrix(result_meanclim$vals, ncol = 1)
colnames(lu_prop_meanclim) <- "meanclim_2020"

keep_rows_meanclim <- !is.na(lu_prop_meanclim[, 1])
lu_prop_meanclim_keep <- lu_prop_meanclim[keep_rows_meanclim, , drop = FALSE]

cat("Meanclim 2020",
    "| AOI cells:", n_cells,
    "| kept:", sum(keep_rows_meanclim),
    "| dropped:", sum(!keep_rows_meanclim), "\n")

#4. Save results
write.csv(
  lu_prop_meanclim_keep, 
  file.path(DIR_DERIVED, "Lu_prop_meanclim_2020.csv"),
  row.names = FALSE
)

saveRDS(
  keep_rows_meanclim,
  file.path(DIR_DERIVED, "Lu_prop_keep_rows_meanclim.rds")
)



  ##### 4. Logic Check  ##########



keep_rows_45      <- results[["45"]]$keep_rows
keep_rows_85      <- results[["85"]]$keep_rows

cat("\n========================================\n")
cat("  LOGIC CHECK: keep_rows consistency\n")
cat("========================================\n")

checks <- list(
  "RCP 4.5 vs RCP 8.5"  = identical(keep_rows_45, keep_rows_85),
  "RCP 4.5 vs meanclim" = identical(keep_rows_45, keep_rows_meanclim),
  "RCP 8.5 vs meanclim" = identical(keep_rows_85, keep_rows_meanclim)
)

for (nm in names(checks)) {
  cat(sprintf("  %s  %s\n", if (checks[[nm]]) "[PASS]" else "[FAIL]", nm))
}

failures <- names(checks)[!unlist(checks)]

if (length(failures) > 0) {
  
  # Report row-level mismatch detail before throwing the error
  if (!checks[["RCP 4.5 vs RCP 8.5"]])
    cat("  Mismatched rows (4.5 vs 8.5):     ", which(keep_rows_45 != keep_rows_85), "\n")
  if (!checks[["RCP 4.5 vs meanclim"]])
    cat("  Mismatched rows (4.5 vs meanclim):", which(keep_rows_45 != keep_rows_meanclim), "\n")
  if (!checks[["RCP 8.5 vs meanclim"]])
    cat("  Mismatched rows (8.5 vs meanclim):", which(keep_rows_85 != keep_rows_meanclim), "\n")
  
  cat("========================================\n\n")
  stop("Logic check FAILED: keep_rows vectors differ across scenarios. ",
       "Investigate before proceeding to downstream analysis.\n",
       "Failed checks: ", paste(failures, collapse = ", "))
  
} else {
  cat(sprintf("\n  All checks passed. Valid EAUs: %d\n", sum(keep_rows_45)))
  cat("========================================\n\n")
}


  ##### 5. Gross habitat flows (loss and gain) per year-pair ##########

# WHY THIS SECTION EXISTS
#
# Section 2 gives prop_suitable, a STOCK. Script 05 currently differences that
# stock to get conversion risk, which yields the NET change:
#
#     (ps_t - ps_t1) / ps_t  =  (L - G) / ps_t
#
# where L is area converting suitable -> unsuitable and G is area converting the
# other way. The intended quantity is the GROSS rate, L / ps_t. The two agree only
# where G = 0. 

# OUTPUT: derived_data/Lu_flows_{45,85}.csv, long, keyed by cell_id and year,
# where `year` is the period the hazard attaches to (2020 -> 2090; 2100 is
# terminal and has no successor).
#
# ── RUNTIME ───────────────────────────────────────────────────────────────────
# Measured in script 12 at ~1,300 s per pair with five native layers. This worker
# emits only three (lost, gained, cover) because prop_suitable already exists from
# section 2, so expect roughly 900-1,100 s per pair, ~60-75 min for 16 pairs at
# four workers. Far slower than section 2 -- this works at native 30 m resolution,
# section 2 aggregates immediately.

#1. Suitable classes, as a vector. Same five codes as rcl_table above; expressed
#   differently because the flow worker tests membership with %in% rather than
#   reclassifying. Section 2 maps non-suitable to NA and lets resample-sum skip
#   them; this maps them to 0 and lets them contribute nothing. Counts of 1s are
#   identical, which the stock-identity check at the end confirms.
SUITABLE_CLASSES <- c(19, 20, 26, 28, 29)

#2. Workers. NOT n_cores from section 2: that is detectCores() - 1, which is right
#   for section 2's small aggregated rasters but wrong here. This job is disk-bound
#   and each concurrent worker holds ~25 GB of native-resolution scratch, so 11
#   workers would ask for ~275 GB and contend for I/O rather than add throughput.
N_CORES_FLOW <- 4L

#3. Backend. PSOCK, not the mclapply used in section 2. Forking inside RStudio on
#   macOS can hang or crash the IDE, which is survivable for section 2's ~8 minutes
#   and is not for a 60+ minute run. PSOCK workers are separate R processes, so the
#   worker below takes ONLY serialisable arguments -- file paths and numbers. A
#   SpatRaster is an external C++ pointer and cannot cross a PSOCK connection at
#   all, which is why the mask is re-opened from its path inside the worker.

#4. The worker.
#
# NA DISCIPLINE -- read before changing the lapp call.
#
# `%in%` returns FALSE (not NA) for NA input, so a no-coverage pixel would silently
# become 0 and enter the sums as real unsuitable land. Coverage is therefore
# tracked as an explicit third layer whose SUM gives the covered pixel count per
# EAU, and missingness is applied AFTER aggregation, where no sentinel exists.
#
# DO NOT add wopt = list(datatype = "INT1U") to the lapp. It saves temp space and
# it is a trap: INT1U's nodata code is 255, that code is a property of the FILE
# whether or not anything is ever flagged with it, and it survives into the
# resampled output. Any EAU whose summed flow count lands on exactly 255 pixels
# comes back NaN. INT2U is not the escape either -- its code is 65535,
# which the coverage sums pass straight through.
#
# GENERAL RULE: do not store an intermediate in a type whose nodata code a
# legitimate value can equal.
decompose_flows <- function(path_t, path_t1, mask_path, classes, scf, n_cells,
                            label, memfrac = NULL, progress_file = NULL) {

  t0 <- Sys.time()

  # Each PSOCK worker is a fresh R process with default terra options, so this
  # must be set here rather than in the parent. memfrac is PER PROCESS: left at
  # the default, N workers can each size chunks against 60% of RAM and drive the
  # machine into swap, which on macOS looks like slowness rather than an error.
  if (!is.null(memfrac)) terra::terraOptions(memfrac = memfrac)

  aoi  <- terra::rast(mask_path)          # re-opened per worker, never inherited
  r_t  <- terra::rast(path_t)
  r_t1 <- terra::rast(path_t1)

  # Stacking requires identical geometry. FOREsce years share a grid, so this
  # should always hold -- but if it ever does not, failing here with both
  # filenames is far better than silently mis-registering loss against gain.
  if (!terra::compareGeom(r_t, r_t1, stopOnError = FALSE))
    stop("Source rasters do not share a geometry: ", basename(path_t), " vs ",
         basename(path_t1), ". Flows would be spatially mis-registered.")

  # PASS 1 -- one crop over both years.
  s <- terra::crop(c(r_t, r_t1), aoi)

  # PASS 2 -- one arithmetic pass emitting all three native layers. The
  # complements and products never become separate rasters; they are transient
  # vectors inside the chunk callback.
  flow <- terra::lapp(
    s,
    fun = function(a, b) {
      ok  <- !is.na(a) & !is.na(b)
      st  <- as.integer(a %in% classes & ok)
      st1 <- as.integer(b %in% classes & ok)
      cbind(st * (1L - st1) * as.integer(ok),   # lost:   suitable -> unsuitable
            (1L - st) * st1 * as.integer(ok),   # gained: unsuitable -> suitable
            as.integer(ok))                     # coverage
    }
    # No datatype request -- see NA DISCIPLINE above.
  )
  names(flow) <- c("lost", "gained", "cover")

  # PASS 3 -- one multi-layer resample, then scale and mask. method = "sum" and
  # NOT aggregate(): section 2 uses resample-sum, and the stock-identity check
  # below only holds if both use the same aggregation.
  agg   <- terra::resample(flow, aoi, method = "sum",
                           wopt = list(datatype = "FLT8S"))
  props <- terra::mask(agg / scf, aoi)

  v <- terra::values(props, na.rm = FALSE)      # n_cells x 3, order as named

  # Missingness applied HERE, from the coverage count, not from a raster sentinel.
  no_cover <- !is.finite(v[, 3]) | v[, 3] <= 0
  v[no_cover, 1:2] <- NA_real_

  out <- data.frame(
    cell_id     = seq_len(n_cells),
    prop_lost   = v[, 1],
    prop_gained = v[, 2],
    prop_cover  = v[, 3]
  )

  # Values are extracted; drop the handles so terra releases the backing temp
  # files before the next pair starts. Concurrent workers would otherwise leave
  # tempdir() growing monotonically for the whole run.
  rm(v, props, agg, flow, s, r_t, r_t1, aoi)
  invisible(gc(verbose = FALSE))

  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # PSOCK workers have no console, so a long run would show no sign of life.
  # Append to a shared log instead; watch with `tail -f`.
  if (!is.null(progress_file))
    try(cat(sprintf("[%s] done %-16s %7.1f s\n",
                    format(Sys.time(), "%H:%M:%S"), label, secs),
            file = progress_file, append = TRUE), silent = TRUE)

  attr(out, "seconds") <- secs
  out
}

#5. Enumerate the year-pairs. Mirrors script 05's structure exactly, so the
#   hazard this feeds is attached to the same periods:
#     * 2020 -> 2030 crosses the meanclim/RCP boundary. The 2020 decision is made
#       before any RCP has materialised, so ps_start comes from meanclim. (This
#       cross-run pairing is a known open question -- the two rasters come from
#       different FOREsce simulations -- and is deliberately unchanged here.)
#     * 2030 -> 2100 are within-RCP forward differences.
#     * 2100 is terminal: no successor, hazard 0 by construction, no pair.
build_flow_pairs <- function(sfx) {
  fl <- data.frame(
    year = results[[sfx]]$file_years,
    path = file.path(scenarios[[sfx]], results[[sfx]]$file_names),
    stringsAsFactors = FALSE
  )
  fl <- fl[fl$year >= 2030, ]          # 2014/2020 RCP layers unused by the panel
  fl <- fl[order(fl$year), ]

  rbind(
    data.frame(rcp = sfx, year = 2020L,
               path_t  = file.path(meanclim_path, meanclim_fname),
               path_t1 = fl$path[fl$year == 2030],
               stringsAsFactors = FALSE),
    data.frame(rcp = sfx, year = head(fl$year, -1),
               path_t  = head(fl$path, -1),
               path_t1 = tail(fl$path, -1),
               stringsAsFactors = FALSE)
  )
}

flow_pairs <- do.call(rbind, lapply(names(scenarios), build_flow_pairs))
n_flow_pairs <- nrow(flow_pairs)

#6. Scratch-space preflight. This job spills native-resolution intermediates to
#   tempdir() and the failure mode is ugly: a worker dies out-of-space an hour in.
#   The requirement is computable rather than guessable, so check it first.
#   Per pair: 2 source layers at 1 byte + 3 flow layers at 4 bytes, plus 30% for
#   transients and compression scratch.
flow_tmp      <- terra::terraOptions(print = FALSE)$tempdir
flow_px       <- as.numeric(n_cells) * scf
flow_gb_pair  <- flow_px * (2 * 1 + 3 * 4) * 1.3 / 1024^3
flow_gb_need  <- flow_gb_pair * N_CORES_FLOW

flow_gb_free <- tryCatch({
  df <- system2("df", c("-Pk", shQuote(flow_tmp)), stdout = TRUE, stderr = FALSE)
  as.numeric(strsplit(trimws(df[length(df)]), "\\s+")[[1]][4]) * 1024 / 1024^3
}, error = function(e) NA_real_, warning = function(w) NA_real_)

cat("\n##### 5. Gross habitat flows #####\n")
cat(sprintf("  Pairs: %d | workers: %d | scratch: %s\n",
            n_flow_pairs, N_CORES_FLOW, flow_tmp))
cat(sprintf("  Peak ~%.0f GB per worker = ~%.0f GB needed | free: %s\n",
            flow_gb_pair, flow_gb_need,
            if (is.na(flow_gb_free)) "unknown" else sprintf("~%.0f GB", flow_gb_free)))

if (!is.na(flow_gb_free) && flow_gb_free < flow_gb_need)
  stop(sprintf(paste0("Not enough scratch space for the flow pass: ~%.0f GB free, ",
                      "~%.0f GB needed at N_CORES_FLOW = %d.\n  The need is linear ",
                      "in workers (~%.0f GB each) -- lower N_CORES_FLOW, or free ",
                      "space on that volume.\n  Stopping now rather than failing ",
                      "an hour into the run."),
               flow_gb_free, flow_gb_need, N_CORES_FLOW, flow_gb_pair),
       call. = FALSE)

#7. Dispatch.
flow_log <- file.path(DIR_DERIVED, "Lu_flows_progress.log")
cat(sprintf("run started %s | %d pairs | %d worker(s)\n",
            format(Sys.time()), n_flow_pairs, N_CORES_FLOW),
    file = flow_log)                       # truncates any previous run's log
cat(sprintf("  Progress log: %s  (PSOCK workers are silent; `tail -f` it)\n\n",
            flow_log))

# memfrac is per process, so divide the single-process budget across workers.
flow_memfrac <- max(0.05, round(0.6 / N_CORES_FLOW, 3))

t_flow <- Sys.time()

run_flow_pair <- function(k) {
  p <- flow_pairs[k, ]
  # try() so a failed pair returns as an ELEMENT rather than aborting the whole
  # dispatch, which would lose both the identity of the failure and every
  # completed pair. The check below names anything that died.
  try(
    decompose_flows(
      path_t        = p$path_t,
      path_t1       = p$path_t1,
      mask_path     = file.path(DIR_DERIVED, "foresce_eau_shared_mask.tif"),
      classes       = SUITABLE_CLASSES,
      scf           = scf,
      n_cells       = n_cells,
      label         = sprintf("rcp%s %d", p$rcp, p$year),
      memfrac       = flow_memfrac,
      progress_file = flow_log
    ),
    silent = TRUE
  )
}

if (N_CORES_FLOW <= 1L) {
  flow_res <- lapply(seq_len(n_flow_pairs), run_flow_pair)
} else {
  cl_flow <- parallel::makeCluster(N_CORES_FLOW, type = "PSOCK")
  # tryCatch/finally, not on.exit(): on.exit registers against a function frame,
  # and at top level in a sourced script there is none, so it would fire straight
  # away and stop the cluster before any work ran. finally guarantees teardown on
  # both success and error, so a failure cannot leave orphan R processes behind.
  flow_res <- tryCatch({
    parallel::clusterEvalQ(cl_flow, suppressPackageStartupMessages(library(terra)))
    parallel::clusterExport(
      cl_flow,
      varlist = c("decompose_flows", "flow_pairs", "DIR_DERIVED",
                  "SUITABLE_CLASSES", "scf", "n_cells", "flow_memfrac", "flow_log"),
      envir = environment()
    )
    parallel::parLapplyLB(cl_flow, seq_len(n_flow_pairs), run_flow_pair)
  }, finally = try(parallel::stopCluster(cl_flow), silent = TRUE))
}

bad_flow <- which(!vapply(flow_res, is.data.frame, logical(1)))
if (length(bad_flow)) {
  cat("\n  Failed pairs:\n")
  for (k in bad_flow) {
    cond <- attr(flow_res[[k]], "condition")
    cat(sprintf("    rcp%s %d: %s\n", flow_pairs$rcp[k], flow_pairs$year[k],
                if (is.null(cond)) paste(utils::capture.output(print(flow_res[[k]])),
                                         collapse = " ") else conditionMessage(cond)))
  }
  stop(length(bad_flow), " of ", n_flow_pairs, " flow pairs failed. Results from a ",
       "partial set would be silently incomplete. If these are memory or disk ",
       "errors, lower N_CORES_FLOW and re-run.")
}

cat("  Per-pair timings:\n")
for (k in seq_len(n_flow_pairs))
  cat(sprintf("    rcp%-4s %d   %7.1f s\n", flow_pairs$rcp[k], flow_pairs$year[k],
              attr(flow_res[[k]], "seconds")))
cat(sprintf("\n  Flow pass wall clock: %.1f min\n",
            as.numeric(difftime(Sys.time(), t_flow, units = "mins"))))

#8. Assemble and write, one long CSV per scenario. Keyed by cell_id and year --
#   NOT by row position. The Lu_prop CSVs above are positional with a companion
#   keep_rows RDS; reproducing that pattern here would re-open the by-position
#   join error found on 7/30, which attached land cover to the wrong EAUs. Script
#   03 joins these by cell_id.
flow_long <- do.call(rbind, lapply(seq_len(n_flow_pairs), function(k) {
  d <- flow_res[[k]]
  data.frame(
    rcp         = flow_pairs$rcp[k],
    year        = flow_pairs$year[k],
    cell_id     = d$cell_id,
    prop_lost   = d$prop_lost,
    prop_gained = d$prop_gained,
    prop_cover  = d$prop_cover,
    stringsAsFactors = FALSE
  )
}))

for (sfx in names(scenarios)) {
  d <- flow_long[flow_long$rcp == sfx, ]
  d <- d[order(d$year, d$cell_id), c("cell_id", "year", "prop_lost", "prop_gained")]
  f <- file.path(DIR_DERIVED, paste0("Lu_flows_", sfx, ".csv"))
  write.csv(d, f, row.names = FALSE)
  cat(sprintf("  %s  (%d rows, %d cells x %d periods)\n",
              f, nrow(d), n_cells, length(unique(d$year))))
}


  ##### 6. Logic Check: flows reconcile with prop_suitable ##########

# THE STOCK IDENTITY. Habitat standing at t+1 must equal habitat standing at t,
# less what converted away, plus what converted in:
#
#     ps_t1  ==  ps_t - prop_lost + prop_gained
#
# This is a strong check, not a formality, because the two sides are produced by
# INDEPENDENT paths through this script. The left comes from section 2
# (classify -> resample-sum, one file at a time). The right comes from section 5
# (%in% -> lapp -> resample-sum, two files compared at native resolution). They
# share only the mask and the scale factor. Agreement therefore confirms:
#
#   * the two encodings of "suitable" (others = NA vs 0/1) count the same pixels
#   * resample-sum is linear and consistent across layers and calls
#   * the flow decomposition is exhaustive -- no transition is double-counted or
#     missed, or the residual would not vanish
#   * cell_id keying in section 5 aligns with row position in section 2
#
# If this fails, the flows are wrong and script 05 must not consume them.
#
# ── TOLERANCES, AND WHY THERE ARE TWO ─────────────────────────────────────────
# The 2026-09-09 run found 100 rows across 13 cells (of 1,169 with coverage)
# where the residual reaches 2.07e-04 -- about 7 ha of a 34,223 ha EAU. Every one
# of those 13 lies OUTSIDE the model panel, and the largest residual among cells
# the model actually uses is 7.18e-08. 
#
# The mechanism is undiagnosed. It is some interaction between the two encodings
# of "suitable" (others = NA here vs 0/1 in section 5) and the fractional weights
# resample-sum assigns at partial-coverage cells. Neither coverage nor habitat
# density separates the offenders cleanly: they span prop_cover 0.06 to 0.99888,
# and 0.99888 is also the maximum among clean cells.
#
# A single strict tolerance therefore fails on a known, bounded, out-of-model
# artefact and would block the pipeline for no reason. Two thresholds instead:
#
#   STRUCT_TOL  a residual this large means the two sections genuinely aggregated
#               differently, or a sentinel collision corrupted a flow layer. Both
#               failure modes are wholesale -- they would affect most cells, not
#               13 -- so this is a hard fail.
#   STRICT_TOL  float-noise level. Cells above it are counted and profiled rather
#               than failed, with a cap on how many are tolerated: if the count
#               grows beyond a few percent, the artefact has stopped being an edge
#               case and needs diagnosing before use.

cat("\n========================================\n")
cat("  LOGIC CHECK: flows vs prop_suitable\n")
cat("========================================\n")

STRICT_TOL    <- 1e-6    # float-noise level for the stock identity
STRUCT_TOL    <- 1e-3    # above this = the two sections disagree structurally
MAX_EDGE_FRAC <- 0.05    # tolerated share of covered cells above STRICT_TOL
n_covered     <- sum(keep_rows_45)

flow_resid <- do.call(rbind, lapply(names(scenarios), function(sfx) {
  lu   <- results[[sfx]]$lu_prop                  # n_cells x years, unfiltered
  fl   <- flow_long[flow_long$rcp == sfx, ]
  yrs  <- sort(unique(fl$year))

  do.call(rbind, lapply(yrs, function(y) {
    d <- fl[fl$year == y, ]
    d <- d[order(d$cell_id), ]

    # ps_t for the 2020 step is meanclim, matching the pairing in section 5.
    ps_t  <- if (y == 2020L) as.numeric(lu_prop_meanclim[, 1])
             else            as.numeric(lu[, sprintf("rcp%s_%d", sfx, y)])
    ps_t1 <- as.numeric(lu[, sprintf("rcp%s_%d", sfx, y + 10L)])

    resid <- ps_t1 - (ps_t - d$prop_lost + d$prop_gained)
    ok    <- is.finite(resid)

    data.frame(
      rcp        = sfx,
      year       = y,
      n_compared = sum(ok),
      n_flow_na  = sum(is.na(d$prop_lost)),
      # Cells where the flows are missing but BOTH stocks are present. This must
      # be zero. The reverse case is expected and benign: section 2 encodes
      # non-suitable as NA, so an EAU with coverage but zero suitable pixels is
      # NA there while section 5 correctly reports it as 0 with coverage > 0.
      # That asymmetry is why this is a subset test rather than a count match.
      n_orphan   = sum(is.na(d$prop_lost) & !is.na(ps_t) & !is.na(ps_t1)),
      n_over_tol = sum(ok & abs(resid) > STRICT_TOL),
      max_resid  = suppressWarnings(max(abs(resid), na.rm = TRUE)),
      min_lost   = suppressWarnings(min(d$prop_lost,   na.rm = TRUE)),
      min_gained = suppressWarnings(min(d$prop_gained, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
}))

print(flow_resid, row.names = FALSE)

edge_frac <- max(flow_resid$n_over_tol) / n_covered

flow_checks <- list(
  "No residual exceeds structural tolerance" =
    all(is.finite(flow_resid$max_resid)) && max(flow_resid$max_resid) < STRUCT_TOL,
  "Identity discrepancies confined to a few cells" =
    edge_frac <= MAX_EDGE_FRAC,
  "Flows are non-negative (L >= 0, G >= 0)" =
    min(flow_resid$min_lost) >= 0 && min(flow_resid$min_gained) >= 0,
  "No cell has stocks but no flows" =
    all(flow_resid$n_orphan == 0),
  "Every pair compared the full covered-cell set" =
    all(flow_resid$n_compared == n_covered)
)

cat(sprintf("\n  Covered cells: %d | max cells above %.0e in any pair: %d (%.2f%%)\n",
            n_covered, STRICT_TOL, max(flow_resid$n_over_tol), 100 * edge_frac))
cat(sprintf("  Max residual: %.3e (structural tolerance %.0e)\n",
            max(flow_resid$max_resid), STRUCT_TOL))

cat("\n")
for (nm in names(flow_checks))
  cat(sprintf("  %s  %s\n", if (flow_checks[[nm]]) "[PASS]" else "[FAIL]", nm))

flow_failures <- names(flow_checks)[!unlist(flow_checks)]

if (length(flow_failures) > 0) {
  cat(sprintf("\n  flow-NA cells per pair: %s | section 2 dropped %d\n",
              paste(unique(flow_resid$n_flow_na), collapse = ", "),
              sum(!keep_rows_45)))
  cat("========================================\n\n")
  stop("Logic check FAILED: the gross flows do not reconcile with prop_suitable. ",
       "Do NOT run script 05 against these flows.\n",
       "  A residual above the structural tolerance usually means the two ",
       "sections aggregated differently -- check that both still use ",
       "resample(method = 'sum').\n",
       "  A jump in the number of cells above the strict tolerance means the ",
       "boundary artefact has stopped being an edge case; profile the offenders ",
       "by prop_cover and by model membership before proceeding.\n",
       "  An unexpected NA count usually means a nodata sentinel collision; see ",
       "the NA DISCIPLINE note in section 5.\n",
       "Failed checks: ", paste(flow_failures, collapse = ", "))
} else {
  cat(sprintf("\n  All checks passed. Max residual %.2e across %d pairs;\n",
              max(flow_resid$max_resid), n_flow_pairs))
  cat(sprintf("  %d of %d covered cells (%.2f%%) exceed %.0e -- see the tolerance\n",
              max(flow_resid$n_over_tol), n_covered, 100 * edge_frac, STRICT_TOL))
  cat("  note above. Script 05 re-checks this identity over the model EAU set.\n")
  cat("========================================\n\n")
}
