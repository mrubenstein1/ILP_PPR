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
    data_panel %>%
      filter(rcp == map_rcp, gcm == map_gcm) %>%
      select(eau_id, year, scaled_abundance),
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

# check for missing scaled_abundance by year
duck_panel_sf %>%
  st_drop_geometry() %>%
  group_by(year) %>%
  summarise(
    n_total = n(),
    n_missing = sum(is.na(scaled_abundance)),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 6. Plot faceted panel with grid overlay
# ------------------------------------------------------------
p_duck_panel <- ggplot() +
  geom_sf(
    data = duck_panel_sf,
    aes(fill = scaled_abundance),
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
eau_change <- data_panel %>%
  group_by(eau_id, rcp, gcm) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    delta_pct = 100 * (scaled_abundance - lag(scaled_abundance)) / lag(scaled_abundance)
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
  filename = paste0("output/EAU_delta_scaled_abundance_panel_rcp", map_rcp, "_", map_gcm, ".png"),
  plot = p_delta_panel,
  width = 12,
  height = 10,
  dpi = 300
)



library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

# ------------------------------------------------------------
# 1. Pick a random EAU (set seed for reproducibility)
# ------------------------------------------------------------
set.seed(42)
chosen_eau <- sample(unique(data_panel$eau_id), 1)
cat("Plotting EAU:", chosen_eau, "\n")

# ------------------------------------------------------------
# 2. Prepare data: build scenario label, subset, pivot long
# ------------------------------------------------------------
plot_data <- data_panel %>%
  filter(eau_id == chosen_eau) %>%
  mutate(
    scenario = case_when(
      rcp == "stationary" ~ "Stationary",
      TRUE ~ paste0("RCP ", ifelse(rcp == "45", "4.5", "8.5"),
                    " \u2013 ", str_to_title(gcm))
    )
  ) %>%
  select(year, scenario, prop_suitable, scaled_abundance) %>%
  pivot_longer(
    cols      = c(prop_suitable, scaled_abundance),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    prop_suitable    = "Proportion Suitable Habitat",
                    scaled_abundance = "Scaled Breeding Pair Abundance (EAU-level)")
  )

# ------------------------------------------------------------
# 3. Define shared color + linetype scales
# ------------------------------------------------------------
scenario_colors <- c(
  "RCP 4.5 \u2013 Dry" = "#2166ac",
  "RCP 4.5 \u2013 Wet" = "#74add1",
  "RCP 8.5 \u2013 Dry" = "#d73027",
  "RCP 8.5 \u2013 Wet" = "#f46d43",
  "Stationary"          = "#4d9221"
)

scenario_linetypes <- c(
  "RCP 4.5 \u2013 Dry" = "solid",
  "RCP 4.5 \u2013 Wet" = "dashed",
  "RCP 8.5 \u2013 Dry" = "solid",
  "RCP 8.5 \u2013 Wet" = "dashed",
  "Stationary"          = "dotted"
)

# ------------------------------------------------------------
# 4. Build plots per metric, combine with patchwork
# ------------------------------------------------------------
make_panel <- function(data, chosen_metric, show_x_label = FALSE) {
  
  sub_data <- data %>% filter(metric == chosen_metric)
  
  ggplot(sub_data,
         aes(x = year, y = value,
             color    = scenario,
             linetype = scenario)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_color_manual(values    = scenario_colors,    name = "Scenario") +
    scale_linetype_manual(values = scenario_linetypes, name = "Scenario") +
    scale_x_continuous(breaks = seq(2020, 2100, by = 10)) +
    labs(
      title = chosen_metric,
      x     = if (show_x_label) "Year" else NULL,
      y     = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 11),
      legend.position  = "none",
      axis.text.x      = if (show_x_label) element_text() else element_blank(),
      axis.ticks.x     = if (show_x_label) element_line() else element_blank(),
      panel.grid.minor = element_blank()
    )
}

p_suitable  <- make_panel(plot_data, "Proportion Suitable Habitat",                show_x_label = FALSE)
p_abundance <- make_panel(plot_data, "Scaled Breeding Pair Abundance (EAU-level)", show_x_label = TRUE)

# Combine panels with shared legend
combined <- (p_suitable / p_abundance) +
  plot_annotation(
    title    = paste0("Scenario Divergence Over Time \u2013 EAU ", chosen_eau),
    subtitle = paste0("WMD: ", unique(data_panel$wmd_id[data_panel$eau_id == chosen_eau])),
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  ) &
  theme(legend.position = "bottom") &
  guides(color    = guide_legend(nrow = 2),
         linetype = guide_legend(nrow = 2))

print(combined)

# ------------------------------------------------------------
# 5. Save
# ------------------------------------------------------------
ggsave(
  filename = paste0("output/eau_", chosen_eau, "_scenario_divergence.png"),
  plot     = combined,
  width    = 9,
  height   = 7,
  dpi      = 300
)

cat("Plot saved to output/eau_", chosen_eau, "_scenario_divergence.png\n", sep = "")




library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

# ============================================================
# SHARED SETUP
# ============================================================

# Build clean scenario labels
eau_plot_base <- data_panel %>%
  mutate(
    scenario = case_when(
      rcp == "stationary" ~ "Stationary",
      TRUE ~ paste0("RCP ", ifelse(rcp == "45", "4.5", "8.5"),
                    " \u2013 ", str_to_title(gcm))
    )
  )

scenario_colors <- c(
  "RCP 4.5 \u2013 Dry" = "#2166ac",
  "RCP 4.5 \u2013 Wet" = "#74add1",
  "RCP 8.5 \u2013 Dry" = "#d73027",
  "RCP 8.5 \u2013 Wet" = "#f46d43",
  "Stationary"          = "#4d9221"
)

# ============================================================
# PLOT 1: SPAGHETTI PLOT WITH SUMMARY RIBBON
# ============================================================

# Compute per-EAU trajectories and per-scenario summary stats
spaghetti_data <- eau_plot_base %>%
  select(eau_id, year, scenario, scaled_abundance)

ribbon_data <- spaghetti_data %>%
  group_by(scenario, year) %>%
  summarise(
    median_abd = median(scaled_abundance, na.rm = TRUE),
    q25_abd    = quantile(scaled_abundance, 0.25, na.rm = TRUE),
    q75_abd    = quantile(scaled_abundance, 0.75, na.rm = TRUE),
    .groups    = "drop"
  )

p_spaghetti <- ggplot() +
  # Individual EAU lines (faint)
  geom_line(
    data    = spaghetti_data,
    aes(x = year, y = scaled_abundance, group = eau_id, color = scenario),
    alpha   = 0.08,
    linewidth = 0.3
  ) +
  # IQR ribbon
  geom_ribbon(
    data = ribbon_data,
    aes(x = year, ymin = q25_abd, ymax = q75_abd, fill = scenario),
    alpha = 0.3
  ) +
  # Median line
  geom_line(
    data = ribbon_data,
    aes(x = year, y = median_abd, color = scenario),
    linewidth = 1.1
  ) +
  facet_wrap(~ scenario, nrow = 2) +
  scale_color_manual(values = scenario_colors, guide = "none") +
  scale_fill_manual(values  = scenario_colors, guide = "none") +
  scale_x_continuous(breaks = seq(2020, 2100, by = 20)) +
  labs(
    title    = "Breeding Pair Abundance Over Time by Scenario",
    subtitle = "Faint lines = individual EAUs; ribbon = IQR; bold line = median",
    x        = "Year",
    y        = "Scaled Breeding Pair Abundance"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold")
  )

print(p_spaghetti)

ggsave(
  "output/abundance_spaghetti_ribbon.png",
  plot   = p_spaghetti,
  width  = 10,
  height = 7,
  dpi    = 300
)

# ============================================================
# PLOT 2: CHANGE MAP (delta abundance 2100 vs 2020)
# ============================================================

# Join spatial coordinates from lookup table
# (x_coord, y_coord are in eau_wmd/eau_lookup)
eau_coords <- eau_wmd %>%
  select(eau_id, x_coord, y_coord)

# Compute delta: scaled_abundance_2100 - scaled_abundance_2020
# Exclude stationary (by definition it won't change)
delta_data <- eau_plot_base %>%
  filter(scenario != "Stationary") %>%
  filter(year %in% c(2020, 2100)) %>%
  select(eau_id, year, scenario, scaled_abundance) %>%
  pivot_wider(
    names_from  = year,
    values_from = scaled_abundance,
    names_prefix = "abd_"
  ) %>%
  mutate(delta = abd_2100 - abd_2020) %>%
  left_join(eau_coords, by = "eau_id")

# Symmetric color scale centered on 0
delta_limit <- max(abs(delta_data$delta), na.rm = TRUE)

p_map <- ggplot(delta_data,
                aes(x = x_coord, y = y_coord, color = delta)) +
  geom_point(size = 1.2, alpha = 0.85) +
  facet_wrap(~ scenario, nrow = 2) +
  scale_color_distiller(
    palette  = "RdBu",
    direction = 1,                      # blue = gain, red = loss
    limits   = c(-delta_limit, delta_limit),
    name     = "\u0394 Abundance\n(2100 \u2013 2020)"
  ) +
  coord_equal() +
  labs(
    title    = "Change in Breeding Pair Abundance: 2100 vs. 2020",
    subtitle = "Blue = gain; Red = loss",
    x        = NULL,
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right"
  )

print(p_map)

ggsave(
  "output/abundance_change_map.png",
  plot   = p_map,
  width  = 10,
  height = 7,
  dpi    = 300
)

cat("Plots saved to output/\n")