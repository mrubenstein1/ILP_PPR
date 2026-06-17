

# ============================================================
# Map % suitable habitat by EAU over time for both RCPs
# Uses the file/object names created in your scripts:
# - input_data/wmd_raster_equal_area.tif
# - input_data/foresce_eau_shared_mask.tif
# - input_data/Lu_prop_45.csv
# - input_data/eau_wmd_lookup.csv
# ============================================================

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)


####### I. MAP % SUITABLE HABITAT OVER TIME: RCP 45 ##########

# ------------------------------------------------------------
# 1. Load the EAU raster that matches the land-cover overlap
# ------------------------------------------------------------
# IMPORTANT:
# Lu_prop_45.csv only contains EAUs with full LC data across all years.
# Those rows correspond to the non-NA cells in the shared AOI mask,
# in raster cell order, because your script built Lu_prop_45.csv as:
#   lu_prop_all_years <- lu_prop[keep_rows, , drop = FALSE]
# where rows are the AOI raster cells.
aoi_mask <- rast("input_data/foresce_eau_shared_mask.tif")

# Optional: full EAU raster if you want WMD outlines or checks
wmd_r <- rast("input_data/wmd_raster_equal_area.tif")

# ------------------------------------------------------------
# 2. Read the suitable habitat proportions for RCP 45
# ------------------------------------------------------------
lu_45 <- read_csv("input_data/Lu_prop_45.csv", show_col_types = FALSE)

# Expect columns like:
# rcp45_2014, rcp45_2015, ..., etc.
year_cols <- grep("^rcp45_[0-9]{4}$", names(lu_45), value = TRUE)

if (length(year_cols) == 0) {
  stop("No columns matching the pattern '^rcp45_[0-9]{4}$' were found in Lu_prop_45.csv")
}

# ------------------------------------------------------------
# 3. Build spatial polygons for the valid EAUs only
# ------------------------------------------------------------
# ------------------------------------------------------------
# 3. Build spatial polygons for the EAUs retained in Lu_prop_45.csv
# ------------------------------------------------------------
aoi_df <- as.data.frame(aoi_mask, xy = TRUE, cells = TRUE, na.rm = FALSE)

# keep_rows is indexed to ALL cells in aoi_mask, not just non-NA WMD cells
keep_rows <- readRDS("input_data/Lu_prop_keep_rows_45.rds")

if (length(keep_rows) != nrow(aoi_df)) {
  stop(
    "Length mismatch:\n",
    "length(keep_rows) = ", length(keep_rows), "\n",
    "rows in full aoi_mask dataframe = ", nrow(aoi_df)
  )
}

# Apply keep_rows to the full raster-cell table first
aoi_keep <- aoi_df %>%
  mutate(keep_rows = keep_rows) %>%
  filter(keep_rows) %>%
  mutate(row_id = row_number())

# Optional check: all retained rows should correspond to actual EAUs
if (any(is.na(aoi_keep$WMD))) {
  stop("Some retained rows have NA WMD values; check alignment between keep_rows and aoi_mask.")
}

if (nrow(aoi_keep) != nrow(lu_45)) {
  stop(
    "Row mismatch after applying keep_rows:\n",
    "retained rows = ", nrow(aoi_keep), "\n",
    "rows in Lu_prop_45.csv = ", nrow(lu_45)
  )
}

# Convert all cells to polygons, then apply keep_rows in the same order
eau_polys <- st_as_sf(as.polygons(aoi_mask, dissolve = FALSE, na.rm = FALSE)) %>%
  mutate(keep_rows = keep_rows) %>%
  filter(keep_rows) %>%
  mutate(row_id = row_number())

# Join the map data
eau_map_wide <- eau_polys %>%
  left_join(
    bind_cols(
      aoi_keep %>% st_drop_geometry() %>% select(row_id, cell, x, y, WMD),
      lu_45
    ),
    by = "row_id"
  )






# ------------------------------------------------------------
# 4. Pivot to long format for faceted mapping
# ------------------------------------------------------------
eau_map_long <- eau_map_wide %>%
  pivot_longer(
    cols = all_of(year_cols),
    names_to = "time_period",
    values_to = "prop_suitable"
  ) %>%
  mutate(
    year = str_extract(time_period, "[0-9]{4}"),
    pct_suitable = prop_suitable * 100
  )

# ------------------------------------------------------------
# 5. Optional WMD outlines for context
# ------------------------------------------------------------
# This uses your original shapefile name from create_EAUs.R
wmd <- st_read(
  "input_data/Wetland_Management_Districts/FSMS_WMD.shp",
  quiet = TRUE
)

# Reproject if needed so outlines match polygons
if (st_crs(wmd) != st_crs(eau_map_long)) {
  wmd <- st_transform(wmd, st_crs(eau_map_long))
}

# ------------------------------------------------------------
# 6. Make the faceted map
# ------------------------------------------------------------
p_rcp45 <- ggplot() +
  geom_sf(
    data = eau_map_long,
    aes(fill = pct_suitable),
    color = NA
  ) +
  geom_sf(
    data = wmd,
    fill = NA,
    color = "grey25",
    linewidth = 0.2
  ) +
  facet_wrap(~ year, ncol = 3) +
  scale_fill_viridis_c(
    name = "% suitable\nhabitat",
    option = "C",
    limits = c(0, 100)
  ) +
  coord_sf(datum = NA) +
  labs(
    title = "Suitable habitat across EAUs over time (RCP 45)",
    subtitle = "Equal Area Units (EAUs); values show percent suitable habitat per EAU",
    caption = "Source objects: foresce_eau_shared_mask.tif and Lu_prop_45.csv"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_rcp45)

# ------------------------------------------------------------
# 7. Save output
# ------------------------------------------------------------
ggsave(
  filename = "output/EAU_pct_suitable_rcp45_faceted_map.png",
  plot = p_rcp45,
  width = 14,
  height = 10,
  dpi = 300
)

# ------------------------------------------------------------
# 8. Optional: map a single year
# ------------------------------------------------------------
map_one_year <- function(yr) {
  ggplot(filter(eau_map_long, year == yr)) +
    geom_sf(aes(fill = pct_suitable), color = NA) +
    geom_sf(data = wmd, fill = NA, color = "grey25", linewidth = 0.2) +
    scale_fill_viridis_c(
      name = "% suitable\nhabitat",
      option = "C",
      limits = c(0, 100)
    ) +
    coord_sf(datum = NA) +
    labs(
      title = paste("Suitable habitat across EAUs - RCP 45 -", yr),
      subtitle = "Percent suitable habitat per EAU"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()
    )
}

# Example:
# print(map_one_year("2014"))
# ggsave("output/EAU_pct_suitable_rcp45_2014.png", map_one_year("2014"),
#        width = 8, height = 6, dpi = 300)


####### II. MAP % SUITABLE HABITAT OVER TIME: RCP 85 ##########

# ------------------------------------------------------------
# 1. Load the EAU raster that matches the land-cover overlap
# ------------------------------------------------------------
# IMPORTANT:
# Lu_prop_85.csv only contains EAUs with full LC data across all years.
# Those rows correspond to the non-NA cells in the shared AOI mask,
# in raster cell order, because your script built Lu_prop_85.csv as:
#   lu_prop_all_years <- lu_prop[keep_rows, , drop = FALSE]
# where rows are the AOI raster cells.
aoi_mask <- rast("input_data/foresce_eau_shared_mask.tif")

# Optional: full EAU raster if you want WMD outlines or checks
wmd_r <- rast("input_data/wmd_raster_equal_area.tif")

# ------------------------------------------------------------
# 2. Read the suitable habitat proportions for RCP 85
# ------------------------------------------------------------
lu_85 <- read_csv("input_data/Lu_prop_85.csv", show_col_types = FALSE)

# Expect columns like:
# rcp85_2014, rcp85_2015, ..., etc.
year_cols <- grep("^rcp85_[0-9]{4}$", names(lu_85), value = TRUE)

if (length(year_cols) == 0) {
  stop("No columns matching the pattern '^rcp85_[0-9]{4}$' were found in Lu_prop_85.csv")
}

# ------------------------------------------------------------
# 3. Build spatial polygons for the valid EAUs only
# ------------------------------------------------------------
# ------------------------------------------------------------
# 3. Build spatial polygons for the EAUs retained in Lu_prop_85.csv
# ------------------------------------------------------------
aoi_df <- as.data.frame(aoi_mask, xy = TRUE, cells = TRUE, na.rm = FALSE)

# keep_rows is indexed to ALL cells in aoi_mask, not just non-NA WMD cells
keep_rows <- readRDS("input_data/Lu_prop_keep_rows_85.rds")

if (length(keep_rows) != nrow(aoi_df)) {
  stop(
    "Length mismatch:\n",
    "length(keep_rows) = ", length(keep_rows), "\n",
    "rows in full aoi_mask dataframe = ", nrow(aoi_df)
  )
}

# Apply keep_rows to the full raster-cell table first
aoi_keep <- aoi_df %>%
  mutate(keep_rows = keep_rows) %>%
  filter(keep_rows) %>%
  mutate(row_id = row_number())

# Optional check: all retained rows should correspond to actual EAUs
if (any(is.na(aoi_keep$WMD))) {
  stop("Some retained rows have NA WMD values; check alignment between keep_rows and aoi_mask.")
}

if (nrow(aoi_keep) != nrow(lu_85)) {
  stop(
    "Row mismatch after applying keep_rows:\n",
    "retained rows = ", nrow(aoi_keep), "\n",
    "rows in Lu_prop_85.csv = ", nrow(lu_85)
  )
}

# Convert all cells to polygons, then apply keep_rows in the same order
eau_polys <- st_as_sf(as.polygons(aoi_mask, dissolve = FALSE, na.rm = FALSE)) %>%
  mutate(keep_rows = keep_rows) %>%
  filter(keep_rows) %>%
  mutate(row_id = row_number())

# Join the map data
eau_map_wide <- eau_polys %>%
  left_join(
    bind_cols(
      aoi_keep %>% st_drop_geometry() %>% select(row_id, cell, x, y, WMD),
      lu_85
    ),
    by = "row_id"
  )






# ------------------------------------------------------------
# 4. Pivot to long format for faceted mapping
# ------------------------------------------------------------
eau_map_long <- eau_map_wide %>%
  pivot_longer(
    cols = all_of(year_cols),
    names_to = "time_period",
    values_to = "prop_suitable"
  ) %>%
  mutate(
    year = str_extract(time_period, "[0-9]{4}"),
    pct_suitable = prop_suitable * 100
  )

# ------------------------------------------------------------
# 5. Optional WMD outlines for context
# ------------------------------------------------------------
# This uses your original shapefile name from create_EAUs.R
wmd <- st_read(
  "input_data/Wetland_Management_Districts/FSMS_WMD.shp",
  quiet = TRUE
)

# Reproject if needed so outlines match polygons
if (st_crs(wmd) != st_crs(eau_map_long)) {
  wmd <- st_transform(wmd, st_crs(eau_map_long))
}

# ------------------------------------------------------------
# 6. Make the faceted map
# ------------------------------------------------------------
p_rcp85 <- ggplot() +
  geom_sf(
    data = eau_map_long,
    aes(fill = pct_suitable),
    color = NA
  ) +
  geom_sf(
    data = wmd,
    fill = NA,
    color = "grey25",
    linewidth = 0.2
  ) +
  facet_wrap(~ year, ncol = 3) +
  scale_fill_viridis_c(
    name = "% suitable\nhabitat",
    option = "C",
    limits = c(0, 100)
  ) +
  coord_sf(datum = NA) +
  labs(
    title = "Suitable habitat across EAUs over time (RCP 85)",
    subtitle = "Equal Area Units (EAUs); values show percent suitable habitat per EAU",
    caption = "Source objects: foresce_eau_shared_mask.tif and Lu_prop_85.csv"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_rcp85)

# ------------------------------------------------------------
# 7. Save output
# ------------------------------------------------------------
ggsave(
  filename = "output/EAU_pct_suitable_rcp85_faceted_map.png",
  plot = p_rcp85,
  width = 14,
  height = 10,
  dpi = 300
)

# ------------------------------------------------------------
# 8. Optional: map a single year
# ------------------------------------------------------------
map_one_year <- function(yr) {
  ggplot(filter(eau_map_long, year == yr)) +
    geom_sf(aes(fill = pct_suitable), color = NA) +
    geom_sf(data = wmd, fill = NA, color = "grey25", linewidth = 0.2) +
    scale_fill_viridis_c(
      name = "% suitable\nhabitat",
      option = "C",
      limits = c(0, 100)
    ) +
    coord_sf(datum = NA) +
    labs(
      title = paste("Suitable habitat across EAUs - RCP 85 -", yr),
      subtitle = "Percent suitable habitat per EAU"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()
    )
}

# Example:
# print(map_one_year("2014"))
# ggsave("output/EAU_pct_suitable_rcp85_2014.png", map_one_year("2014"),
#        width = 8, height = 6, dpi = 300)



#### III. Sanity check: look at absolute % change over time per EAU:85 ######
# ------------------------------------------------------------
# 9. Calculate change from first to last RCP 85 time step
# ------------------------------------------------------------

# Identify the RCP 85 columns using the naming convention from Lu_prop_85.csv
year_cols <- grep("^rcp85_[0-9]{4}$", names(eau_map_wide), value = TRUE)

if (length(year_cols) < 2) {
  stop("Need at least two RCP 85 time steps to calculate change.")
}

# Order columns by year
year_vals <- as.integer(sub("rcp85_", "", year_cols))
year_cols <- year_cols[order(year_vals)]
year_vals <- sort(year_vals)

first_col <- year_cols[1]
last_col  <- year_cols[length(year_cols)]

cat("Calculating delta as:", last_col, "-", first_col, "\n")

# Calculate delta in proportion and percentage points
eau_delta <- eau_map_wide %>%
  mutate(
    delta_prop = .data[[last_col]] - .data[[first_col]],
    delta_pct_point = delta_prop * 100
  )

# Optional quick summary
summary(eau_delta$delta_pct_point)

# ------------------------------------------------------------
# 10. Create delta map
# ------------------------------------------------------------

# Make the color scale symmetric around 0
max_abs_change <- max(abs(eau_delta$delta_pct_point), na.rm = TRUE)

p_delta_rcp85 <- ggplot() +
  geom_sf(
    data = eau_delta,
    aes(fill = delta_pct_point),
    color = NA
  ) +
  geom_sf(
    data = wmd,
    fill = NA,
    color = "grey25",
    linewidth = 0.2
  ) +
  scale_fill_gradient2(
    name = "Change in %\nsuitable habitat",
    low = "#b2182b",    # red = loss
    mid = "white",
    high = "#1a9850",   # green = gain
    midpoint = 0,
    limits = c(-max_abs_change, max_abs_change)
  ) +
  coord_sf(datum = NA) +
  labs(
    title = "Change in suitable habitat across EAUs (RCP 85)",
    subtitle = paste0(
      "Delta = ", last_col, " - ", first_col,
      " (percentage-point change)"
    ),
    caption = "Negative values indicate habitat loss; positive values indicate habitat gain"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "right"
  )

print(p_delta_rcp85)

# ------------------------------------------------------------
# 11. Save delta map
# ------------------------------------------------------------

ggsave(
  filename = "output/EAU_delta_pct_suitable_rcp85.png",
  plot = p_delta_rcp85,
  width = 9,
  height = 7,
  dpi = 300
)


#### IV. Sanity check: look at absolute % change over time per EAU:45 ######
# ------------------------------------------------------------
# 9. Calculate change from first to last RCP 45 time step
# ------------------------------------------------------------

# Identify the RCP 45 columns using the naming convention from Lu_prop_45.csv
year_cols <- grep("^rcp45_[0-9]{4}$", names(eau_map_wide), value = TRUE)

if (length(year_cols) < 2) {
  stop("Need at least two RCP 45 time steps to calculate change.")
}

# Order columns by year
year_vals <- as.integer(sub("rcp45_", "", year_cols))
year_cols <- year_cols[order(year_vals)]
year_vals <- sort(year_vals)

first_col <- year_cols[1]
last_col  <- year_cols[length(year_cols)]

cat("Calculating delta as:", last_col, "-", first_col, "\n")

# Calculate delta in proportion and percentage points
eau_delta <- eau_map_wide %>%
  mutate(
    delta_prop = .data[[last_col]] - .data[[first_col]],
    delta_pct_point = delta_prop * 100
  )

# Optional quick summary
summary(eau_delta$delta_pct_point)

# ------------------------------------------------------------
# 10. Create delta map
# ------------------------------------------------------------

# Make the color scale symmetric around 0
max_abs_change <- max(abs(eau_delta$delta_pct_point), na.rm = TRUE)

p_delta_rcp45 <- ggplot() +
  geom_sf(
    data = eau_delta,
    aes(fill = delta_pct_point),
    color = NA
  ) +
  geom_sf(
    data = wmd,
    fill = NA,
    color = "grey25",
    linewidth = 0.2
  ) +
  scale_fill_gradient2(
    name = "Change in %\nsuitable habitat",
    low = "#b2182b",    # red = loss
    mid = "white",
    high = "#1a9450",   # green = gain
    midpoint = 0,
    limits = c(-max_abs_change, max_abs_change)
  ) +
  coord_sf(datum = NA) +
  labs(
    title = "Change in suitable habitat across EAUs (RCP 45)",
    subtitle = paste0(
      "Delta = ", last_col, " - ", first_col,
      " (percentage-point change)"
    ),
    caption = "Negative values indicate habitat loss; positive values indicate habitat gain"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "right"
  )

print(p_delta_rcp45)

# ------------------------------------------------------------
# 11. Save delta map
# ------------------------------------------------------------

ggsave(
  filename = "output/EAU_delta_pct_suitable_rcp45.png",
  plot = p_delta_rcp45,
  width = 9,
  height = 7,
  dpi = 300
)

