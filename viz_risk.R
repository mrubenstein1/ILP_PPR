### Transition Probability Visualizations ###
#
# Explores how trans_prob varies over time and space across EAUs.
#
# PREREQUISITES: Script 05 must have been run (trans_prob populated in data_panel).
# SPATIAL SECTION: requires EAU centroid coordinates — see Section 2 setup note.

library(dplyr)
library(ggplot2)
library(scales)
library(ggridges)


# ── Shared theme and palette ──────────────────────────────────────────────────

rcp_colors  <- c("45" = "#2171b5", "85" = "#cb181d")
rcp_labels  <- c("45" = "RCP 4.5", "85" = "RCP 8.5")

theme_thesis <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )


# ── Deduplicated trans_prob series (one row per EAU × year × RCP) ─────────────
# Drops GCM dimension since trans_prob is identical for wet/dry within an RCP.

tp_df <- data_panel %>%
  filter(rcp %in% c("45", "85"), year < 2100) %>%
  distinct(eau_id, wmd_id, year, rcp, trans_prob) %>%
  mutate(rcp_label = rcp_labels[rcp],
         rcp_label = factor(rcp_label, levels = rcp_labels))


# ════════════════════════════════════════════════════════════════════════════════
# SECTION 1: TEMPORAL PLOTS
# ════════════════════════════════════════════════════════════════════════════════

# ── Plot 1: Boxplots of trans_prob distribution over time by RCP ──────────────

p1 <- tp_df %>%
  ggplot(aes(x = factor(year), y = trans_prob, fill = rcp)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3, linewidth = 0.4) +
  scale_fill_manual(values = rcp_colors, labels = rcp_labels, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  facet_wrap(~ rcp_label, ncol = 1) +
  labs(
    title    = "Distribution of Transition Probability Across EAUs Over Time",
    subtitle = "Each box = distribution across all EAUs at that decade (terminal 2100 excluded)",
    x        = "Year",
    y        = "Transition Probability (Risk of Loss)"
  ) +
  theme_thesis +
  theme(legend.position = "none")

print(p1)


# ── Plot 2: Mean trajectory with ± 1 SD ribbon by RCP ────────────────────────

p2 <- tp_df %>%
  group_by(rcp_label, year) %>%
  summarise(mean_tp = mean(trans_prob),
            sd_tp   = sd(trans_prob),
            .groups = "drop") %>%
  ggplot(aes(x = year, y = mean_tp, color = rcp_label, fill = rcp_label)) +
  geom_ribbon(aes(ymin = pmax(0, mean_tp - sd_tp),
                  ymax = mean_tp + sd_tp),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_color_manual(values = rcp_colors, labels = rcp_labels, name = NULL) +
  scale_fill_manual(values  = rcp_colors, labels = rcp_labels, name = NULL) +
  scale_x_continuous(breaks = seq(2030, 2090, 10)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), limits = c(0, NA)) +
  labs(
    title    = "Mean Transition Probability Over Time by RCP",
    subtitle = "Ribbon = ± 1 SD across EAUs",
    x        = "Year",
    y        = "Mean Transition Probability"
  ) +
  theme_thesis

print(p2)


# ── Plot 3: Ridgeline — how the shape of the distribution shifts over time ────

p3 <- tp_df %>%
  ggplot(aes(x = trans_prob, y = factor(year), fill = rcp)) +
  geom_density_ridges(alpha = 0.6, scale = 0.9, rel_min_height = 0.01) +
  scale_fill_manual(values = rcp_colors, labels = rcp_labels, name = NULL) +
  scale_x_continuous(labels = percent_format(accuracy = 0.1), limits = c(0, NA)) +
  facet_wrap(~ rcp_label) +
  labs(
    title    = "Transition Probability Distribution by Decade",
    subtitle = "Ridgelines show how the cross-EAU distribution shifts shape over time",
    x        = "Transition Probability",
    y        = "Year"
  ) +
  theme_thesis +
  theme(legend.position = "none")

print(p3)


# ── Plot 4: High-risk EAU trajectories ───────────────────────────────────────
# Highlights the EAUs that ever exceed the 90th percentile of trans_prob.
# Useful for identifying the most urgently at-risk parcels.

p90_threshold <- quantile(tp_df$trans_prob, 0.90)

high_risk_eaus <- tp_df %>%
  filter(trans_prob >= p90_threshold) %>%
  distinct(eau_id)

p4 <- tp_df %>%
  mutate(high_risk = eau_id %in% high_risk_eaus$eau_id) %>%
  ggplot(aes(x = year, y = trans_prob,
             group = interaction(eau_id, rcp),
             color = high_risk, alpha = high_risk)) +
  geom_line(linewidth = 0.3) +
  geom_hline(yintercept = p90_threshold, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  annotate("text", x = 2031, y = p90_threshold,
           label = sprintf("90th pctile = %.1f%%", p90_threshold * 100),
           vjust = -0.5, hjust = 0, size = 3) +
  scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "#d94801"),
                     labels = c("FALSE" = "All other EAUs",
                                "TRUE"  = "Ever above 90th percentile"),
                     name = NULL) +
  scale_alpha_manual(values = c("FALSE" = 0.15, "TRUE" = 0.7), guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = seq(2030, 2090, 10)) +
  facet_wrap(~ rcp_label) +
  labs(
    title    = "High-Risk EAU Trajectories",
    subtitle = "Orange = EAUs that exceed the 90th percentile at any point; grey = all others",
    x        = "Year",
    y        = "Transition Probability"
  ) +
  theme_thesis

print(p4)


# ── Plot 5: Spaghetti — one line per EAU ─────────────────────────────────────
# Individual EAU trajectories with group mean overlaid in bold.
# Coloured by RCP; alpha kept low so the mean line reads clearly over the mass.

p5 <- tp_df %>%
  ggplot(aes(x = year, y = trans_prob,
             group = interaction(eau_id, rcp),
             color = rcp)) +
  geom_line(alpha = 0.08, linewidth = 0.25) +
  stat_summary(aes(group = rcp), fun = mean,
               geom = "line", linewidth = 1.2) +
  scale_color_manual(values = rcp_colors, labels = rcp_labels, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = seq(2030, 2090, 10)) +
  facet_wrap(~ rcp_label) +
  labs(
    title    = "Transition Probability Trajectories — All EAUs",
    subtitle = "Thin lines = individual EAUs; bold line = group mean",
    x        = "Year",
    y        = "Transition Probability"
  ) +
  theme_thesis +
  theme(legend.position = "none")

print(p5)


# ════════════════════════════════════════════════════════════════════════════════
# SECTION 2: SPATIAL MAPS
# ════════════════════════════════════════════════════════════════════════════════
#
# ── Setup: join EAU centroids ─────────────────────────────────────────────────
# These maps need EAU centroid coordinates from eau_lookup (Script 01).
#
# CHECK: run names(eau_lookup) and update the two lines below to match your
# actual coordinate column names before running this section.

coord_x_col <- "x_coord"
coord_y_col <- "y_coord"

eau_coords <- eau_wmd %>%
  select(eau_id, x = all_of(coord_x_col), y = all_of(coord_y_col))

tp_spatial <- tp_df %>%
  left_join(eau_coords, by = "eau_id")

# Confirm join worked before plotting
stopifnot("Coordinate join failed — check coord_x_col and coord_y_col above" =
            !all(is.na(tp_spatial$x)))


# ── Plot 6: Spatial map — mean trans_prob per EAU across all decades ──────────

p6 <- tp_spatial %>%
  group_by(eau_id, rcp_label, x, y) %>%
  summarise(mean_tp = mean(trans_prob), .groups = "drop") %>%
  ggplot(aes(x = x, y = y, color = mean_tp)) +
  geom_point(size = 1.2, alpha = 0.85) +
  scale_color_viridis_c(
    name   = "Mean\ntrans_prob",
    option = "plasma",
    labels = percent_format(accuracy = 0.1)
  ) +
  coord_equal() +
  facet_wrap(~ rcp_label) +
  labs(
    title    = "Spatial Pattern of Mean Transition Probability (2030–2090)",
    subtitle = "Each point = one EAU; colour = mean trans_prob averaged across all decades",
    x        = NULL, y = NULL
  ) +
  theme_thesis +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

print(p6)


# ── Plot 7: Spatial small multiples — trans_prob at each decade ───────────────
# Faceted by year × RCP. Shows how the spatial risk pattern evolves over time.

p7 <- tp_spatial %>%
  ggplot(aes(x = x, y = y, color = trans_prob)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_viridis_c(
    name   = "trans_prob",
    option = "plasma",
    labels = percent_format(accuracy = 0.1)
  ) +
  coord_equal() +
  facet_grid(rcp_label ~ year) +
  labs(
    title    = "Transition Probability Across Space and Time",
    subtitle = "Rows = RCP scenario; columns = decade",
    x        = NULL, y = NULL
  ) +
  theme_thesis +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 8))

print(p7)


# ── Plot 8: Spatial map — change in trans_prob from early to late century ─────
# Highlights EAUs where risk of loss accelerates vs. diminishes over time.

p8 <- tp_spatial %>%
  filter(year %in% c(2030, 2090)) %>%
  select(eau_id, rcp_label, x, y, year, trans_prob) %>%
  tidyr::pivot_wider(names_from = year, values_from = trans_prob,
                     names_prefix = "tp_") %>%
  mutate(tp_change = tp_2090 - tp_2030) %>%
  ggplot(aes(x = x, y = y, color = tp_change)) +
  geom_point(size = 1.2, alpha = 0.85) +
  scale_color_gradient2(
    low      = "#2171b5",
    mid      = "grey90",
    high     = "#cb181d",
    midpoint = 0,
    name     = "Change in\ntrans_prob\n(2090 − 2030)",
    labels   = percent_format(accuracy = 0.1)
  ) +
  coord_equal() +
  facet_wrap(~ rcp_label) +
  labs(
    title    = "Change in Transition Probability from 2030 to 2090",
    subtitle = "Red = risk increasing over century; blue = risk decreasing",
    x        = NULL, y = NULL
  ) +
  theme_thesis +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

print(p8)