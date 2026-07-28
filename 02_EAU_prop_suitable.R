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
if (!isTRUE(.SETUP_DONE)) source("00_setup.R")

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
  
  
  
  
  
  
  
  
