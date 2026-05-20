library(terra)
library(sf)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load base raster, lookup, and WMD outlines
# ------------------------------------------------------------
wmd_r <- rast("input_data/wmd_raster_equal_area.tif")
eau_lookup <- readRDS("input_data/eau_wmd_lookup.rds")
wmd <- st_read("input_data/Wetland_Management_Districts/FSMS_WMD.shp", quiet = TRUE)

# ------------------------------------------------------------
# 2. Rebuild the raster-cell table in the same order used to make eau_lookup
# ------------------------------------------------------------
r_df <- as.data.frame(wmd_r, xy = TRUE, cells = TRUE) %>%
  rename(
    x_coord = x,
    y_coord = y,
    wmd_id  = WMD
  ) %>%
  filter(!is.na(wmd_id)) %>%
  left_join(
    eau_lookup %>% select(eau_id, wmd_id, x_coord, y_coord),
    by = c("wmd_id", "x_coord", "y_coord")
  )

# check: every raster cell should get an eau_id
sum(is.na(r_df$eau_id))

# ------------------------------------------------------------
# 3. Build polygons and attach the same raster cell index
# ------------------------------------------------------------
eau_grid <- st_as_sf(as.polygons(wmd_r, dissolve = FALSE, na.rm = FALSE)) %>%
  mutate(cell = 1:n()) %>%
  filter(!is.na(WMD)) %>%
  left_join(
    r_df %>% select(cell, eau_id, wmd_id),
    by = "cell"
  )

# check: every polygon should get an eau_id
sum(is.na(eau_grid$eau_id))

if (st_crs(wmd) != st_crs(eau_grid)) {
  wmd <- st_transform(wmd, st_crs(eau_grid))
}

# ------------------------------------------------------------
# 4. Choose scenario and join all decades
# ------------------------------------------------------------
map_rcp <- "85"
map_gcm <- "dry"

duck_panel_sf <- eau_grid %>%
  left_join(
    eau_panel_alloc %>%
      filter(rcp == map_rcp, gcm == map_gcm) %>%
      select(eau_id, year, scaled_density),
    by = "eau_id"
  ) %>%
  filter(!is.na(year))

# ------------------------------------------------------------
# 5. Diagnostic checks
# ------------------------------------------------------------
# should be 1 row per EAU per year
duck_panel_sf %>%
  st_drop_geometry() %>%
  count(year)

# check for missing scaled_density by year
duck_panel_sf %>%
  st_drop_geometry() %>%
  group_by(year) %>%
  summarise(
    n_total = n(),
    n_missing = sum(is.na(scaled_density)),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 6. Plot faceted panel with grid overlay
# ------------------------------------------------------------
p_duck_panel <- ggplot() +
  geom_sf(
    data = duck_panel_sf,
    aes(fill = scaled_density),
    color = NA
  ) +
  geom_sf(
    data = duck_panel_sf,
    fill = NA,
    color = "white",
    linewidth = 0.08
  ) +
  geom_sf(
    data = wmd,
    fill = NA,
    color = "grey20",
    linewidth = 0.25
  ) +
  facet_wrap(~ year, ncol = 3) +
  scale_fill_viridis_c(
    name = "Duck density\nper EAU",
    option = "C",
    na.value = "grey90"
  ) +
  coord_sf(datum = NA) +
  labs(
    title = paste("EAU-scale duck density across decades:", "RCP", map_rcp, map_gcm),
    subtitle = "Scaled density distributed across EAUs by suitable habitat share"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_duck_panel)


library(dplyr)
library(ggplot2)
library(sf)

# ------------------------------------------------------------
# 1. Choose scenario
# ------------------------------------------------------------
map_rcp <- "45"
map_gcm <- "dry"

# ------------------------------------------------------------
# 2. Compute percent change from previous time step
# ------------------------------------------------------------
eau_change <- eau_panel_alloc %>%
  group_by(eau_id, rcp, gcm) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    delta_pct = 100 * (scaled_density - lag(scaled_density)) / lag(scaled_density)
  ) %>%
  ungroup() %>%
  filter(rcp == map_rcp, gcm == map_gcm)

# optional: remove first year because it has no previous period
eau_change <- eau_change %>%
  filter(!is.na(delta_pct))

# ------------------------------------------------------------
# 3. Join change values onto EAU polygons
# ------------------------------------------------------------
eau_change_map <- eau_grid %>%
  left_join(
    eau_change %>% select(eau_id, year, rcp, gcm, delta_pct),
    by = "eau_id"
  ) %>%
  filter(!is.na(year))

# ------------------------------------------------------------
# 4. Make color scale symmetric around 0
# ------------------------------------------------------------
max_abs_change <- max(abs(eau_change_map$delta_pct), na.rm = TRUE)

# ------------------------------------------------------------
# 5. Faceted delta map
# ------------------------------------------------------------
p_delta_panel <- ggplot() +
  geom_sf(
    data = eau_change_map,
    aes(fill = delta_pct),
    color = NA
  ) +
  geom_sf(
    data = eau_change_map,
    fill = NA,
    color = "white",
    linewidth = 0.08
  ) +
  geom_sf(
    data = wmd,
    fill = NA,
    color = "grey25",
    linewidth = 0.2
  ) +
  facet_wrap(~ year, ncol = 3) +
  scale_fill_gradient2(
    name = "% change in\nscaled density",
    low = "#b2182b",
    mid = "white",
    high = "#1a9850",
    midpoint = 0,
    limits = c(-max_abs_change, max_abs_change)
  ) +
  coord_sf(datum = NA) +
  labs(
    title = paste("Change in EAU-scale duck density:", "RCP", map_rcp, map_gcm),
    subtitle = "Each panel shows percent change from the previous decade",
    caption = "Negative values indicate decline; positive values indicate increase"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_delta_panel)

# ------------------------------------------------------------
# 6. Save
# ------------------------------------------------------------
ggsave(
  filename = paste0("output/EAU_delta_scaled_density_panel_rcp", map_rcp, "_", map_gcm, ".png"),
  plot = p_delta_panel,
  width = 12,
  height = 10,
  dpi = 300
)