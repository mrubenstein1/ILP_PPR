###ILP Formulation of the PPR Problem ###
  # March 2026

### ESTABLISH EAUs ####
  # spatial unit of analysis is a simulated parcel of uniform dimension (“equal area unit”, EAU)
  # EAU = 282km2 (70,000 acres). 
  # the number of EAUs per WMD varies based on WMD size, but there are on average 50 EAUs per WMD. 
  # There are ~1200 EAUs total across all 24 WMDs

#Credit Heini Kujala for original script; modified by Madeleine Rubenstein March 2026

#Load required libraries for all scripts
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(parallel)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)
library(forcats)

 
##### 1. Import Data ##########
 

#Load WMD shapefile
wmd <- st_read("input_data/Wetland_Management_Districts/FSMS_WMD.shp")

#Calculate target cell size
  #if we want an average of 50 EAUs per WMD, how large does each EAU need to be?
  #these calculations look at the size of WMDs and figure out the km2 distance
  #of each EAU to equal an average of 50 EAUs per WMD. 

total_area <- sum(wmd$area_km2) #total area of all WMDs (i.e., entire PPR region)
divideby <- 50 * length(wmd$WMD) 
cell_size <- sqrt(total_area/divideby) # in km
cell_size <- ceiling( cell_size * 1000 ) # convert cell size from km to m

#Create raster template (at target resolution)
wmd_nrows <- ceiling( (ext(wmd)[4]-ext(wmd)[3]) / cell_size ) #define rows of template
wmd_ncols <- ceiling( (abs( ext(wmd)[1] )+ext(wmd)[2]) / cell_size ) #define columns of template
wmd_r <- rast(wmd, nrow = wmd_nrows, ncol = wmd_ncols) #create the new raster at EAU resolution
wmd_r <- rast(wmd, resolution = cell_size)

#assign WMD ID to each EAU (each cell of this raster = 1 EAU)
wmd_split <- rasterize(wmd, wmd_r, field = "WMD", touches = FALSE) #touches=F means that EAUs are only assigned to the WMD if the center point falls within the polygon

#visualize
plot(wmd_split)

#save as geotiff
#writeRaster(wmd_split, "input_data/wmd_raster_equal_area.tif", filetype = 'GTiff') #this line is commented out because once the tif is saved, it doesn't need to be re-written

 
##### 2. Extract EAU x WMD lookup table ##########
 
# Reload the raster
wmd_r   <- rast("input_data/wmd_raster_equal_area.tif")

# ── 2. Build WMD reference table (numeric ID + area) ─────────────────────────
wmd_ref <- wmd %>%
  st_drop_geometry() %>%
  select(wmd_id = WMD, area_km2) %>%          # adjust if column names differ
  arrange(wmd_id) %>%
  mutate(wmd_id_num = row_number())            # numeric ID 1:N alphabetically

# ── 3. Extract EAU × WMD lookup table from raster ────────────────────────────
eau_wmd <- as.data.frame(wmd_r, xy = TRUE, cells = TRUE) %>%
  rename(
    eau_id  = cell,
    x_coord = x,
    y_coord = y,
    wmd_id  = "WMD"
  ) %>%
  filter(!is.na(wmd_id)) %>%
  left_join(wmd_ref, by = "wmd_id") %>%       # attach numeric ID and area
  arrange(wmd_id_num, eau_id) %>%
  mutate(eau_id = row_number()) %>%            # re-index EAUs 1:N sequentially
  select(eau_id, wmd_id_num, wmd_id, area_km2, x_coord, y_coord)

# ── 4. WMD summary table ──────────────────────────────────────────────────────
eau_summary <- eau_wmd %>%
  group_by(wmd_id_num, wmd_id, area_km2) %>%
  summarise(n_eaus = n(), .groups = "drop") %>%
  arrange(wmd_id_num)

# ── 5. Visualize ──────────────────────────────────────────────────────────────────
plot(wmd_split)
plot(as.polygons(wmd_split, dissolve = FALSE), add = TRUE, border = "white", lwd = 0.3)


# ── 6. Save ───────────────────────────────────────────────────────────────────
saveRDS(eau_wmd,     "input_data/eau_wmd_lookup.rds")
write.csv(eau_wmd,     "input_data/eau_wmd_lookup.csv",   row.names = FALSE)
write.csv(eau_summary, "input_data/wmd_summary.csv",       row.names = FALSE)

 
##### 3. Logic Check ##########
 

checks <- list(
  "No NA values in eau_id"      = !any(is.na(eau_wmd$eau_id)),
  "No NA values in wmd_id"      = !any(is.na(eau_wmd$wmd_id)),
  "No NA values in wmd_id_num"  = !any(is.na(eau_wmd$wmd_id_num)),
  "No NA values in area_km2"    = !any(is.na(eau_wmd$area_km2)),
  "EAU IDs are unique"          = !any(duplicated(eau_wmd$eau_id)),
  "EAU IDs are sequential 1:N"  = all(eau_wmd$eau_id == seq_len(nrow(eau_wmd))),
  "Expected 24 WMDs"            = nrow(eau_summary) == 24,
  "All WMDs have at least 1 EAU" = all(eau_summary$n_eaus > 0),
  "Total EAUs in expected range" = between(nrow(eau_wmd), 1000, 1500)
)

for (nm in names(checks)) {
  cat(sprintf("  %s  %s\n", if (checks[[nm]]) "[PASS]" else "[FAIL]", nm))
}

failures <- names(checks)[!unlist(checks)]

if (length(failures) > 0) {
  
  # Print detail on specific failures before halting
  if (!checks[["No NA values in eau_id"]])
    cat("  NA count in eau_id:     ", sum(is.na(eau_wmd$eau_id)), "\n")
  if (!checks[["No NA values in wmd_id"]])
    cat("  NA count in wmd_id:     ", sum(is.na(eau_wmd$wmd_id)), "\n")
  if (!checks[["No NA values in wmd_id_num"]])
    cat("  NA count in wmd_id_num: ", sum(is.na(eau_wmd$wmd_id_num)), "\n")
  if (!checks[["No NA values in area_km2"]])
    cat("  NA count in area_km2:   ", sum(is.na(eau_wmd$area_km2)), "\n")
  if (!checks[["EAU IDs are unique"]])
    cat("  Duplicate EAU IDs:      ", sum(duplicated(eau_wmd$eau_id)), "\n")
  if (!checks[["EAU IDs are sequential 1:N"]])
    cat("  EAU ID range:            ", min(eau_wmd$eau_id), "to", max(eau_wmd$eau_id),
        "(expected 1 to", nrow(eau_wmd), ")\n")
  if (!checks[["Expected 24 WMDs"]])
    cat("  WMDs found:             ", nrow(eau_summary), "(expected 24)\n")
  if (!checks[["All WMDs have at least 1 EAU"]])
    cat("  WMDs with 0 EAUs:       ", 
        paste(eau_summary$wmd_id[eau_summary$n_eaus == 0], collapse = ", "), "\n")
  if (!checks[["Total EAUs in expected range"]])
    cat("  Total EAUs found:       ", nrow(eau_wmd), "(expected 1000–1500)\n")
  
  cat("========================================\n\n")
  stop("Logic check FAILED: EAU × WMD lookup table has errors. ",
       "Investigate before proceeding to downstream analysis.\n",
       "Failed checks: ", paste(failures, collapse = ", "))
  
} else {
  cat(sprintf("\n  All checks passed. WMDs: %d | Total EAUs: %d\n",
              nrow(eau_summary), nrow(eau_wmd)))}
