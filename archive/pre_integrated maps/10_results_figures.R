# ══════════════════════════════════════════════════════════════════════════════
# 10_figures.R  —  Results table + figures for the PPR acquisition models
# ══════════════════════════════════════════════════════════════════════════════
#
# Reads the two model-output CSVs and writes a results table (flextable -> .docx/.png
# + a plain .csv) and six figures (.png at 300 dpi + vector .pdf) to OUT_DIR.
#
#   table_over_baseline   value each model creates over the do-nothing baseline
#                         (% improvement and Δpairs); all three models
#   table_gap_vs_rolling  greedy/myopic value-added gap vs. rolling (gap % and Δpairs)
#   fig1_gaps         value-added gaps, greedy vs. myopic, by scenario
#   fig2_delta_pairs  absolute pairs left on the table, by scenario
#   fig3_wetdry       wet vs. dry: stakes (A) vs. foresight premium (B)
#   fig4_J_vs_VA      same myopic gap on value-added vs. on total J
#   fig5_overtime_gap cumulative gap vs. rolling over the horizon (climate scenarios)
#   fig6_decline      landscape abundance over time (rolling)
#
# Everything you might restyle (paths, sizes, palette, font) is in the CONFIG block.
# Pure ggplot2/flextable — no model objects required, so it is fully reproducible
# from the CSVs alone.
# ══════════════════════════════════════════════════════════════════════════════

# ── CONFIG ────────────────────────────────────────────────────────────────────
IN_DIR   <- "output_data"     # folder holding the two CSVs (08's OUT_DIR)
OUT_DIR  <- "output_figs"     # where tables + figures are written
RESULTS  <- file.path(IN_DIR, "model_results.csv")
TRAJ     <- file.path(IN_DIR, "model_trajectories.csv")

DPI         <- 300
BASE_SIZE   <- 12
BASE_FAMILY <- ""             # "" = system default; set e.g. "Helvetica" if installed
SAVE_PDF    <- TRUE           # also write vector PDF alongside each PNG
TABLE_PNG   <- TRUE           # try to write a PNG of the table (needs webshot2 + magick)

# palette
COL_MODEL    <- c(greedy = "#B0413E", myopic = "#E2A23B", rolling = "#2E7D5B")
COL_MODEL_TXT<- c(greedy = "#B0413E", myopic = "#B8801E")          # darker amber for text
COL_MOIST    <- c(wet = "#3F7FA6", dry = "#C98A3A")
COL_METRIC   <- c(`value-added` = "#7E57C2", `total J` = "#B0BEC5")
GREEN_HEADER <- "#2E7D5B"
# ──────────────────────────────────────────────────────────────────────────────

## packages -------------------------------------------------------------------
.need <- c("readr", "dplyr", "tidyr", "ggplot2", "scales", "flextable", "patchwork")
.miss <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
if (length(.miss))
  stop("Install missing package(s): ", paste(.miss, collapse = ", "),
       "\n  install.packages(c(", paste(sprintf('\"%s\"', .miss), collapse = ", "), "))")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(scales); library(flextable); library(patchwork)
})

if (!file.exists(RESULTS) || !file.exists(TRAJ))
  stop("Could not find inputs in IN_DIR = '", IN_DIR, "'.\n  Expected: ",
       RESULTS, " and ", TRAJ, "\n  Edit IN_DIR at the top if they live elsewhere.")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

DELTA <- "\u0394"  # Δ, written portably

## labels / order -------------------------------------------------------------
SCEN_LEVELS <- c("rcp45_wet", "rcp45_dry", "rcp85_wet", "rcp85_dry", "stationary")
SCEN_LAB    <- c(rcp45_wet = "RCP 4.5\nwet", rcp45_dry = "RCP 4.5\ndry",
                 rcp85_wet = "RCP 8.5\nwet", rcp85_dry = "RCP 8.5\ndry",
                 stationary = "Stationary")
SCEN_LAB1   <- c(rcp45_wet = "RCP 4.5 wet", rcp45_dry = "RCP 4.5 dry",
                 rcp85_wet = "RCP 8.5 wet", rcp85_dry = "RCP 8.5 dry",
                 stationary = "Stationary")
CLIM <- SCEN_LEVELS[1:4]

## formatters -----------------------------------------------------------------
fmt_pct <- function(x) vapply(x, function(v) {
  if (is.na(v)) return(NA_character_)
  if (v < 1e-4) return("<0.0001%")
  if (v < 1) { s <- formatC(v, format = "f", digits = 4)
  s <- sub("0+$", "", s); s <- sub("\\.$", "", s); return(paste0(s, "%")) }
  paste0(formatC(v, format = "f", digits = 2), "%")
}, character(1))
fmt_pairs <- function(x) ifelse(abs(x) < 1, "<1", formatC(round(x), format = "d", big.mark = ","))

## theme ----------------------------------------------------------------------
theme_ppr <- function(base_size = BASE_SIZE, base_family = BASE_FAMILY) {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title    = element_text(face = "bold", hjust = 0, size = rel(1.12), margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0, colour = "grey35", margin = margin(b = 8)),
      plot.caption  = element_text(hjust = 0.5, colour = "grey45", size = rel(0.72), margin = margin(t = 8)),
      axis.title    = element_text(colour = "grey25"),
      axis.title.y  = element_text(angle = 90, margin = margin(r = 6)),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "top", legend.justification = "right",
      legend.title = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0, size = rel(0.95)),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(10, 14, 8, 10)
    )
}
theme_set(theme_ppr())

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = DPI, bg = "white")
  if (SAVE_PDF) {
    # Prefer cairo_pdf (nice font embedding) but fall back to the base 'pdf' device if
    # cairo is unavailable. On Windows builds without cairo this surfaces as a WARNING
    # ("failed to load cairo DLL"), not an error, so we catch warnings too and fall back
    # quietly; suppressWarnings on the fallback hides any glyph-substitution notices.
    ok <- tryCatch(
      { ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h,
               device = cairo_pdf, bg = "white"); TRUE },
      warning = function(w) FALSE,
      error   = function(e) FALSE)
    if (!ok)
      suppressWarnings(ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p,
                              width = w, height = h, bg = "white"))
  }
  print(p)            # also render to the active device (RStudio Plots pane)
  invisible(p)
}

## data -----------------------------------------------------------------------
res <- read_csv(RESULTS, show_col_types = FALSE) %>%
  mutate(scenario = factor(scenario, levels = SCEN_LEVELS))
trj <- read_csv(TRAJ, show_col_types = FALSE) %>%
  mutate(scenario = factor(scenario, levels = SCEN_LEVELS))

rollVA <- res %>% filter(model == "rolling") %>% select(scenario, rollVA = value_added)

gaps <- res %>%
  filter(model %in% c("greedy", "myopic")) %>%
  left_join(rollVA, by = "scenario") %>%
  mutate(dpairs = rollVA - value_added,
         model  = factor(model, levels = c("greedy", "myopic")))

# ══════════════════════════════════════════════════════════════════════════════
# TABLES  —  (1) value created over baseline, all models;  (2) gap vs. rolling
# Each written as flextable -> .docx (+ .png if webshot2/magick present) and a .csv.
# ══════════════════════════════════════════════════════════════════════════════

# percent formatter with 2 decimals (for the over-baseline %, which run ~0.3–0.7%)
fmt_pct2 <- function(x) ifelse(x < 1e-4, "<0.0001%", paste0(formatC(x, format = "f", digits = 2), "%"))
DP <- paste0(DELTA, "pairs")

# Build a styled flextable with a two-level header that groups paired columns under a
# model label. `df` has Scenario in column 1 and value columns thereafter; `groups`
# and `subs` are one entry per value column. Model text colours follow the palette.
style_grouped_table <- function(df, groups, subs, caption, footer_text) {
  vcols  <- names(df)[-1]
  rle_g  <- rle(groups)
  top_w  <- c(1, rle_g$lengths)
  top_l  <- c("Scenario", rle_g$values)
  seps   <- c(1, head(1 + cumsum(rle_g$lengths), -1))        # right-edge group separators
  txtcol <- c(COL_MODEL_TXT, rolling = unname(COL_MODEL[["rolling"]]))
  bd     <- officer::fp_border(color = "#dddddd", width = 1)
  
  ft <- flextable(df)
  ft <- set_header_labels(ft, values = setNames(as.list(c("", subs)), c("Scenario", vcols)))
  ft <- add_header_row(ft, top = TRUE, values = top_l, colwidths = top_w)
  ft <- merge_at(ft, i = 1:2, j = 1, part = "header")
  ft <- bg(ft, part = "header", bg = GREEN_HEADER)
  ft <- color(ft, part = "header", color = "white")
  ft <- bold(ft, part = "header")
  ft <- align(ft, align = "center", part = "all")
  ft <- align(ft, j = 1, align = "left", part = "all")
  ft <- valign(ft, j = 1, valign = "center", part = "header")
  ft <- bold(ft, j = 1, part = "body")
  for (g in unique(groups)) {
    gcols <- vcols[groups == g]; col <- txtcol[[tolower(g)]]
    if (!is.null(col) && length(gcols)) ft <- color(ft, j = gcols, color = col, part = "body")
  }
  ft <- bg(ft, i = seq(2, nrow(df), by = 2), bg = "#f5f5f2", part = "body")
  ft <- border_remove(ft)
  ft <- hline(ft, i = 1, part = "header", border = officer::fp_border(color = "white", width = 1))
  ft <- vline(ft, j = seps, border = bd, part = "header")
  ft <- vline(ft, j = seps, border = bd, part = "body")
  ft <- padding(ft, padding.top = 4, padding.bottom = 4, part = "all")
  ft <- add_footer_lines(ft, values = footer_text)
  ft <- fontsize(ft, size = 8.5, part = "footer")
  ft <- color(ft, color = "grey40", part = "footer")
  ft <- set_caption(ft, caption)
  autofit(ft)
}

save_table <- function(ft, df_csv, stem) {
  write_csv(df_csv, file.path(OUT_DIR, paste0(stem, ".csv")))
  save_as_docx(ft, path = file.path(OUT_DIR, paste0(stem, ".docx")))
  if (TABLE_PNG) {
    ok <- tryCatch({ save_as_image(ft, path = file.path(OUT_DIR, paste0(stem, ".png"))); TRUE },
                   error = function(e) FALSE)
    if (!ok) message("  ", stem, ".png skipped (needs webshot2 + magick); .docx and .csv written.")
  }
}

# ── TABLE 1 — value created over the do-nothing baseline (all three models) ────
base1 <- res %>%
  transmute(scenario, model,
            pct = value_added / J_baseline * 100,   # % improvement over baseline
            dp  = value_added) %>%                  # Δpairs over baseline (= value_added)
  pivot_wider(names_from = model, values_from = c(pct, dp)) %>%
  arrange(scenario)

tbl1 <- data.frame(
  Scenario    = unname(SCEN_LAB1[as.character(base1$scenario)]),
  greedy_pct  = fmt_pct2(base1$pct_greedy),  greedy_dp  = fmt_pairs(base1$dp_greedy),
  myopic_pct  = fmt_pct2(base1$pct_myopic),  myopic_dp  = fmt_pairs(base1$dp_myopic),
  rolling_pct = fmt_pct2(base1$pct_rolling), rolling_dp = fmt_pairs(base1$dp_rolling),
  check.names = FALSE, stringsAsFactors = FALSE
)
groups1 <- c("Greedy", "Greedy", "Myopic", "Myopic", "Rolling", "Rolling")
subs1   <- c("% over baseline", DP, "% over baseline", DP, "% over baseline", DP)
csv1    <- setNames(tbl1, c("Scenario", paste(groups1, subs1)))

ft1 <- style_grouped_table(
  tbl1, groups1, subs1,
  caption = "Conservation value created over the do-nothing baseline, by scenario",
  footer_text = paste0("% over baseline = value_added / J_baseline \u00d7 100 (improvement over the ",
                       "do-nothing landscape).  ", DP, " = value_added (discounted duck-pairs gained vs. baseline)."))
save_table(ft1, csv1, "table_over_baseline")

# ── TABLE 2 — value-added gap vs. rolling foresight (greedy, myopic) ───────────
gap2 <- gaps %>%
  transmute(scenario, model, gva = gap_value_added_pct, dp = dpairs) %>%
  pivot_wider(names_from = model, values_from = c(gva, dp)) %>%
  arrange(scenario)

tbl2 <- data.frame(
  Scenario   = unname(SCEN_LAB1[as.character(gap2$scenario)]),
  greedy_gap = fmt_pct(gap2$gva_greedy),  greedy_dp = fmt_pairs(gap2$dp_greedy),
  myopic_gap = fmt_pct(gap2$gva_myopic),  myopic_dp = fmt_pairs(gap2$dp_myopic),
  check.names = FALSE, stringsAsFactors = FALSE
)
groups2 <- c("Greedy", "Greedy", "Myopic", "Myopic")
subs2   <- c("gap %", DP, "gap %", DP)
csv2    <- setNames(tbl2, c("Scenario", paste(groups2, subs2)))

ft2 <- style_grouped_table(
  tbl2, groups2, subs2,
  caption = "Value-added gap vs. rolling (deterministic foresight), by scenario",
  footer_text = paste0("gap % = shortfall as % of rolling's value-added.  ", DP,
                       " = rolling value-added \u2212 model value-added (discounted duck-pairs)."))
save_table(ft2, csv2, "table_gap_vs_rolling")

# ══════════════════════════════════════════════════════════════════════════════
# FIG 1  —  value-added gaps, greedy vs. myopic
# ══════════════════════════════════════════════════════════════════════════════
dodge <- position_dodge(width = 0.78)
fig1 <- ggplot(gaps, aes(scenario, gap_value_added_pct, fill = model)) +
  geom_col(position = dodge, width = 0.72) +
  geom_text(aes(label = fmt_pct(gap_value_added_pct)), position = dodge,
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = COL_MODEL, labels = c("Greedy", "Myopic")) +
  scale_x_discrete(labels = SCEN_LAB) +
  scale_y_continuous("Gap vs. rolling (% of value-added)",
                     limits = c(0, 52), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "How far each policy falls short of deterministic foresight", x = NULL,
       caption = "Rolling = 0% by definition. Stationary null shown for reference (value-added \u2248 0, so % is noise on a near-zero base).")
save_fig(fig1, "fig1_gaps", 9, 5)

# ══════════════════════════════════════════════════════════════════════════════
# FIG 2  —  absolute pairs left on the table
# ══════════════════════════════════════════════════════════════════════════════
fig2 <- ggplot(gaps, aes(scenario, dpairs, fill = model)) +
  geom_col(position = dodge, width = 0.72) +
  geom_text(aes(label = fmt_pairs(dpairs)), position = dodge, vjust = -0.4, size = 3) +
  scale_fill_manual(values = COL_MODEL, labels = c("Greedy", "Myopic")) +
  scale_x_discrete(labels = SCEN_LAB) +
  scale_y_continuous("Foregone value vs. rolling\n(discounted duck-pairs)",
                     labels = label_comma(), expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Absolute cost of myopia and greed (duck-pairs left on the table)", x = NULL)
save_fig(fig2, "fig2_delta_pairs", 9, 5)

# ══════════════════════════════════════════════════════════════════════════════
# FIG 3  —  wet vs. dry: stakes (A) vs. premium (B)
# ══════════════════════════════════════════════════════════════════════════════
wd <- res %>%
  filter(scenario %in% CLIM) %>%
  mutate(rcp      = ifelse(grepl("45", scenario), "RCP 4.5", "RCP 8.5"),
         moisture = factor(ifelse(grepl("wet", scenario), "wet", "dry"), levels = c("wet", "dry")))

pA <- wd %>% filter(model == "rolling") %>%
  ggplot(aes(rcp, value_added, fill = moisture)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_text(aes(label = label_comma()(round(value_added))),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 2.9) +
  scale_fill_manual(values = COL_MOIST) +
  scale_y_continuous("Value at stake\n(rolling value-added, pairs)",
                     labels = label_comma(), expand = expansion(mult = c(0, 0.10))) +
  labs(title = "A  Stakes scale with abundance (wet > dry)", x = NULL)

pB <- wd %>% filter(model == "myopic") %>%
  ggplot(aes(rcp, gap_value_added_pct, fill = moisture)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_text(aes(label = sprintf("%.1f%%", gap_value_added_pct)),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 2.9) +
  scale_fill_manual(values = COL_MOIST) +
  scale_y_continuous("Foresight premium\n(myopic gap, % value-added)",
                     limits = c(0, 18), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "B  Premium is ~flat across wet/dry (risk is shared)", x = NULL)

fig3 <- (pA | pB) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Wet vs. dry: moisture sets the stakes, RCP sets the premium",
    caption = "Hazard (trans_prob) is identical across wet/dry within an RCP (GCM duplicates land cover); only WMD-level abundance differs.",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0, size = rel(1.15)),
                  plot.caption = element_text(hjust = 0.5, colour = "grey45", size = rel(0.72)),
                  legend.position = "top", legend.justification = "right")) &
  theme(legend.title = element_blank())
save_fig(fig3, "fig3_wetdry", 11, 5)

# ══════════════════════════════════════════════════════════════════════════════
# FIG 4  —  same myopic gap on value-added vs. on total J
# ══════════════════════════════════════════════════════════════════════════════
jva <- res %>%
  filter(model == "myopic", scenario %in% CLIM) %>%
  select(scenario, `value-added` = gap_value_added_pct, `total J` = gap_J_pct) %>%
  pivot_longer(-scenario, names_to = "metric", values_to = "gap") %>%
  mutate(metric = factor(metric, levels = c("value-added", "total J")))

fig4 <- ggplot(jva, aes(scenario, gap, fill = metric)) +
  geom_col(position = dodge, width = 0.72) +
  geom_text(aes(label = fmt_pct(gap)), position = dodge, vjust = -0.4, size = 3,
            colour = "grey25") +
  scale_fill_manual(values = COL_METRIC) +
  scale_x_discrete(labels = SCEN_LAB[CLIM]) +
  scale_y_continuous("Myopic gap (%)", limits = c(0, 18), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Why value-added is the right lens: the same gap is invisible on total J", x = NULL)
save_fig(fig4, "fig4_J_vs_VA", 9, 5)

# ══════════════════════════════════════════════════════════════════════════════
# FIG 5  —  cumulative gap vs. rolling over the horizon (climate scenarios)
# ══════════════════════════════════════════════════════════════════════════════
roll_disc <- trj %>% filter(model == "rolling") %>% select(scenario, year, roll = discounted)
cumg <- trj %>%
  filter(model %in% c("greedy", "myopic"), scenario %in% CLIM) %>%
  left_join(roll_disc, by = c("scenario", "year")) %>%
  arrange(scenario, model, year) %>%
  group_by(scenario, model) %>%
  mutate(cum_gap = cumsum(roll - discounted)) %>%
  ungroup() %>%
  mutate(model = factor(model, levels = c("greedy", "myopic")),
         scenario = factor(scenario, levels = CLIM))
ends <- cumg %>% group_by(scenario, model) %>% filter(year == max(year)) %>% ungroup()

fig5 <- ggplot(cumg, aes(year, cum_gap, colour = model)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  geom_text(data = ends, aes(label = fmt_pairs(cum_gap)), hjust = -0.12, size = 2.8, show.legend = FALSE) +
  facet_wrap(~scenario, ncol = 2, labeller = as_labeller(SCEN_LAB1)) +
  scale_colour_manual(values = COL_MODEL, labels = c("Greedy", "Myopic")) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  scale_y_continuous("Cumulative gap vs. rolling\n(discounted duck-pairs)", labels = label_comma()) +
  labs(title = "Where the foresight advantage accrues: cumulative gap builds over the horizon",
       x = "year",
       caption = "Each curve sums the year-by-year prevented-loss shortfall vs. rolling; the 2100 endpoint equals the value-added \u0394pairs.") +
  theme(panel.grid.major.x = element_line(colour = "grey92"))
save_fig(fig5, "fig5_overtime_gap", 10.5, 7.4)

# ══════════════════════════════════════════════════════════════════════════════
# FIG 6  —  landscape decline over time (rolling)
# ══════════════════════════════════════════════════════════════════════════════
SCEN_COL <- c(rcp45_wet = "#3F7FA6", rcp45_dry = "#9BB7C4",
              rcp85_wet = "#9C3848", rcp85_dry = "#D79A86", stationary = "#777777")
dec <- trj %>% filter(model == "rolling") %>%
  mutate(scenario = factor(scenario, levels = SCEN_LEVELS))

fig6 <- ggplot(dec, aes(year, ducks / 1e6, colour = scenario)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
  scale_colour_manual(values = SCEN_COL, labels = SCEN_LAB1) +
  scale_y_continuous("Landscape abundance\n(million duck-pairs)") +
  scale_x_continuous(breaks = seq(2020, 2100, 20)) +
  labs(title = "Landscape decline under each climate future (rolling policy)", x = "year",
       caption = "All scenarios share the 2020 anchor (5.11M). Policies are near-indistinguishable at landscape scale (~0.6% of total) \u2014 the separation lives in value-added.") +
  guides(colour = guide_legend(nrow = 2)) +
  theme(panel.grid.major.x = element_line(colour = "grey92"))
save_fig(fig6, "fig6_decline", 9, 5.2)

message("\nDone. Wrote 2 tables + 6 figures to '", OUT_DIR, "/'.")