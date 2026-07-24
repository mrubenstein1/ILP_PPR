### Visualization: Cost Data (Exploratory) ####
# Four quick plots examining how acquisition cost varies spatially and temporally.
# Run after script 06; requires data_panel and eau_wmd in memory.

# ── Load if running standalone ────────────────────────────────────────────────
# data_panel <- readRDS("input_data/data_panel.rds")
# eau_wmd    <- read_csv("input_data/eau_wmd_lookup.csv", show_col_types = FALSE)


# ── Derive tile size from EAU raster ──────────────────────────────────────────
# geom_tile() needs explicit width/height to fill space without gaps.
# We read these directly from the raster so the map is geometrically correct
# regardless of what CRS Script 01 used.

eau_r  <- rast("input_data/wmd_raster_equal_area.tif")
tile_w <- res(eau_r)[1]
tile_h <- res(eau_r)[2]

cat(sprintf("CRS:       %s\n", crs(eau_r, proj = TRUE)))
cat(sprintf("Tile size: %.0f × %.0f CRS units\n", tile_w, tile_h))


# ── Prepare 2020 baseline cost table ──────────────────────────────────────────
# One row per EAU. Cost is scenario-invariant so any single rcp/gcm will do;
# the 2020 baseline row is the natural reference.

cost_2020 <- data_panel %>%
  filter(year == 2020, rcp == "baseline") %>%
  distinct(eau_id, cost) %>%
  left_join(eau_wmd %>% select(eau_id, x_coord, y_coord, wmd_id), by = "eau_id") %>%
  mutate(cost_m = cost / 1e6)   # work in millions USD throughout


# ── Shared theme ──────────────────────────────────────────────────────────────
theme_set(theme_bw(base_size = 11))


# ── Plot 1: Spatial map of 2020 cost ──────────────────────────────────────────
# Log scale because land values are roughly log-normally distributed.
# Axis text suppressed — raw CRS coordinates (likely metres) aren't informative.

p1 <- ggplot(cost_2020, aes(x = x_coord, y = y_coord, fill = cost_m)) +
  geom_tile(width = tile_w, height = tile_h) +
  scale_fill_viridis_c(
    option = "plasma",
    name   = "Cost\n(M USD)",
    trans  = "log10",
    labels = scales::label_dollar(suffix = "M", accuracy = 1)
  ) +
  coord_equal() +
  labs(
    title    = "EAU Acquisition Cost — 2020 Baseline",
    subtitle = "PLACES FMV: vacant, 2% inflation to 2020; log10 colour scale"
  ) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  )

print(p1)


# ── Plot 2: Distribution of 2020 cost ─────────────────────────────────────────
# Log-scale x-axis to match the log-normal shape of land values.
# Vertical lines mark median and mean.

med_cost <- median(cost_2020$cost_m)
mn_cost  <- mean(cost_2020$cost_m)

p2 <- ggplot(cost_2020, aes(x = cost_m)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white", linewidth = 0.3) +
  geom_vline(xintercept = med_cost, linetype = "solid",  color = "black", linewidth = 0.8) +
  geom_vline(xintercept = mn_cost,  linetype = "dashed", color = "black", linewidth = 0.8) +
  annotate("text", x = med_cost, y = Inf, vjust = 1.5, hjust = -0.1,
           label = sprintf("Median\n$%.0fM", med_cost), size = 3) +
  annotate("text", x = mn_cost,  y = Inf, vjust = 1.5, hjust =  1.1,
           label = sprintf("Mean\n$%.0fM", mn_cost),   size = 3) +
  scale_x_log10(labels = scales::label_dollar(suffix = "M", accuracy = 1)) +
  labs(
    title = "Distribution of EAU Acquisition Costs — 2020",
    x     = "Total acquisition cost (million USD, log scale)",
    y     = "Number of EAUs"
  )

print(p2)


# ── Plot 3: Cost by WMD, ordered by median ────────────────────────────────────
# Shows which districts are expensive vs cheap — useful for interpreting where
# the ILP is likely to face binding budget constraints.

p3 <- cost_2020 %>%
  mutate(wmd_id = forcats::fct_reorder(wmd_id, cost_m, median)) %>%
  ggplot(aes(x = wmd_id, y = cost_m)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.size = 0.8,
               outlier.alpha = 0.5) +
  scale_y_continuous(labels = scales::label_dollar(suffix = "M", accuracy = 1)) +
  coord_flip() +
  labs(
    title = "EAU Acquisition Cost by WMD — 2020",
    x     = NULL,
    y     = "Total acquisition cost (million USD)"
  )

print(p3)


# ── Plot 4: Cost trajectory 2020–2100 ─────────────────────────────────────────
# Since cost is scenario-invariant, we collapse to distinct(eau_id, year, cost).
# The shape of the distribution is preserved — only the scale shifts upward
# at 2% per year. Line = median across EAUs; ribbon = IQR.

cost_traj <- data_panel %>%
  distinct(eau_id, year, cost) %>%
  group_by(year) %>%
  summarise(
    median = median(cost) / 1e6,
    p25    = quantile(cost, 0.25) / 1e6,
    p75    = quantile(cost, 0.75) / 1e6,
    .groups = "drop"
  )

p4 <- ggplot(cost_traj, aes(x = year)) +
  geom_ribbon(aes(ymin = p25, ymax = p75), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = median), color = "steelblue", linewidth = 1) +
  scale_x_continuous(breaks = seq(2020, 2100, by = 10)) +
  scale_y_continuous(labels = scales::label_dollar(suffix = "M", accuracy = 1)) +
  labs(
    title    = "Projected EAU Acquisition Cost — 2020 to 2100",
    subtitle = "Median (line) and IQR (ribbon) across EAUs; 2% annual inflation from 2017 baseline",
    x        = "Decision year",
    y        = "Total acquisition cost (million USD)"
  )

print(p4)



library(patchwork)

# ── Shared spatial setup ───────────────────────────────────────────────────────
eau_r  <- rast("input_data/wmd_raster_equal_area.tif")
tile_w <- res(eau_r)[1]
tile_h <- res(eau_r)[2]


# ── Prepare 2020 baseline data ─────────────────────────────────────────────────
# One row per EAU. Both cost and abundance are scenario-invariant at 2020
# (cost is always scenario-invariant; abundance at 2020 is the observed baseline).

gut_check <- data_panel %>%
  filter(year == 2020, rcp == "baseline") %>%
  distinct(eau_id, scaled_abundance, cost) %>%
  left_join(eau_wmd %>% select(eau_id, x_coord, y_coord, wmd_id), by = "eau_id") %>%
  mutate(benefit_cost_ratio = scaled_abundance / cost)

cat(sprintf("EAUs in gut-check: %d\n", nrow(gut_check)))
cat(sprintf("Correlation (log cost ~ log abundance): r = %.3f\n",
            cor(log(gut_check$cost), log(gut_check$scaled_abundance + 1))))


# ── Shared map theme ───────────────────────────────────────────────────────────
map_theme <- theme_bw() +
  theme(axis.title = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank())


# ── Plot 1a: Cost map ──────────────────────────────────────────────────────────
map_cost <- ggplot(gut_check, aes(x = x_coord, y = y_coord, fill = cost / 1e6)) +
  geom_tile(width = tile_w, height = tile_h) +
  scale_fill_viridis_c(
    option = "plasma",
    name   = "Cost\n(M USD)",
    trans  = "log10",
    labels = scales::label_dollar(suffix = "M", accuracy = 1)
  ) +
  coord_equal() +
  labs(title = "Acquisition Cost (2020)") +
  map_theme


# ── Plot 1b: Duck density map ──────────────────────────────────────────────────
map_ducks <- ggplot(gut_check, aes(x = x_coord, y = y_coord, fill = scaled_abundance)) +
  geom_tile(width = tile_w, height = tile_h) +
  scale_fill_viridis_c(
    option = "viridis",
    name   = "Breeding\npairs",
    trans  = "log10",
    labels = scales::label_comma(accuracy = 1)
  ) +
  coord_equal() +
  labs(title = "Duck Breeding Pairs (2020)") +
  map_theme


# ── Print side by side ─────────────────────────────────────────────────────────
print(map_cost + map_ducks +
        plot_annotation(title    = "Cost vs Duck Density — 2020 Baseline",
                        subtitle = "Log10 colour scales"))


# ── Plot 2: Scatter — cost vs duck density ────────────────────────────────────
# Log-log scale: if the relationship is positive (expensive land = more ducks),
# it suggests conservation targets and high costs co-locate — a difficult budget
# environment. If negative or flat, cheap land may still have high duck value.

p_scatter <- ggplot(gut_check, aes(x = cost / 1e6, y = scaled_abundance)) +
  geom_point(aes(color = wmd_id), alpha = 0.5, size = 1.5, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  scale_x_log10(labels = scales::label_dollar(suffix = "M", accuracy = 1)) +
  scale_y_log10(labels = scales::label_comma(accuracy = 1)) +
  labs(
    title    = "Duck Breeding Pairs vs Acquisition Cost — 2020",
    subtitle = "Each point = one EAU; log-log scale; OLS trend line with 95% CI",
    x        = "Total acquisition cost (million USD, log scale)",
    y        = "Breeding pairs (EAU-level, log scale)"
  ) +
  theme_bw()

print(p_scatter)


# ── Plot 3: Benefit-to-cost ratio map ─────────────────────────────────────────
# Breeding pairs per dollar of acquisition cost.
# This is the spatial analogue of the ILP's efficiency criterion — EAUs in the
# upper tail of this distribution are where the model will want to invest first.
# EAUs that are cheap AND duck-rich should glow brightest here.

p_ratio <- ggplot(gut_check, aes(x = x_coord, y = y_coord,
                                 fill = benefit_cost_ratio)) +
  geom_tile(width = tile_w, height = tile_h) +
  scale_fill_viridis_c(
    option = "magma",
    name   = "Pairs\nper USD",
    trans  = "log10"
  ) +
  coord_equal() +
  labs(
    title    = "Benefit-to-Cost Ratio — 2020",
    subtitle = "Breeding pairs per acquisition dollar (log scale); bright = efficient"
  ) +
  map_theme

print(p_ratio)


# ── Quick tabular summary by WMD ──────────────────────────────────────────────
# Useful for spotting whether any WMDs are systematically expensive-but-duck-rich
# or cheap-but-duck-poor.

cat("\n--- WMD-level summary (2020 baseline, medians across EAUs) ---\n\n")

gut_check %>%
  group_by(wmd_id) %>%
  summarise(
    n_eaus             = n(),
    median_cost_m      = round(median(cost) / 1e6, 1),
    median_pairs       = round(median(scaled_abundance), 1),
    median_bcr         = round(median(benefit_cost_ratio), 8),
    .groups = "drop"
  ) %>%
  arrange(desc(median_bcr)) %>%
  print(n = Inf)


