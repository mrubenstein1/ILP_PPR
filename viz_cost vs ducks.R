### Visualization: Cost vs Duck Density Gut-Check ####
# Three plots examining whether high duck density coincides with expensive land.
# Run after script 06; requires data_panel and eau_wmd in memory.
# Requires: patchwork (install.packages("patchwork") if needed)

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