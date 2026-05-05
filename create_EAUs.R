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

##########################
##### I. IMPORT ##########
##########################

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

##########################
##### II. EXTRACT EAU x WMD LOOKUP TABLE ##########
##########################
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

# ── 5. Check ──────────────────────────────────────────────────────────────────
cat("Total EAUs:", nrow(eau_wmd), "\n")
cat("Total WMDs:", nrow(eau_summary), "\n")
print(eau_summary)

#visualize the EAUs per WMD
plot(wmd_split)
plot(as.polygons(wmd_split, dissolve = FALSE), add = TRUE, border = "white", lwd = 0.3)


# ── 6. Save ───────────────────────────────────────────────────────────────────
saveRDS(eau_wmd,     "input_data/eau_wmd_lookup.rds")
write.csv(eau_wmd,     "input_data/eau_wmd_lookup.csv",   row.names = FALSE)
write.csv(eau_summary, "input_data/wmd_summary.csv",       row.names = FALSE)

