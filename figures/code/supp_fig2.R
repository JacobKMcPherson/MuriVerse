# Supplementary Figure 2 (4 panels, 2 x 2): performance and label reproducibility.
#
# Row 1: a  IV vs pooled route-restriction delta-rho (lollipop/dumbbell)
#         b  Per-species performance: median fold error and Spearman rho, one facet each,
#            horizontal bars over the same three strata
# Row 2: c  OOF predicted vs observed (213-record IV set)
#         d  Between-source label agreement (cross-source pairs scatter)
#
# Data: MuriVerse/rodent/strict/ and figures/tabs/pub-tab__{2,3} for panels a and b (plus the
# frozen results/supp_fig2/ tables), and figures/tabs/table60_muriverse_halflife_cross_source_pairs_iv.csv
# for panel d. Panel c plots seed 0, like Figure 3c, and prints the HEADLINE metrics (mean of
# three seeds); panel b reports the median fold error the same way (mean of the three seeds' medians).
if (!exists("ROOT")) {
  source(file.path(system("git rev-parse --show-toplevel", intern = TRUE),
                   "figures/code/00_setup.R"))
}

DIR_REV  <- file.path(ROOT, "results", "supp_fig2")
DIR_SNAP <- file.path(ROOT, "MuriVerse", "rodent", "strict")
LIVE_TABS <- file.path(ROOT, "figures", "tabs")
DIR_TABS <- LIVE_TABS

# ---------- Load data ----------
meta <- read_csv(file.path(DIR_SNAP, "halflife_meta_with_split.csv"),
                 show_col_types = FALSE, guess_max = 1e5) |>
  mutate(species = str_to_title(species))

oof  <- read_csv(file.path(DIR_SNAP, "oof_predictions_iv.csv"),
                 show_col_types = FALSE) |>
  mutate(species = str_to_title(species))

rt   <- read_csv(file.path(DIR_TABS, "pub-tab__2_route_restriction.csv"),
                 show_col_types = FALSE)

perf <- read_csv(file.path(DIR_TABS, "pub-tab__3_performance_by_stratum.csv"),
                 show_col_types = FALSE)

s12  <- read_csv(file.path(DIR_REV, "S12_species_cell_contrasts_and_studies.csv"),
                 show_col_types = FALSE)

s01  <- read_csv(file.path(DIR_REV, "S01_final_IV_predictions_32_rows.csv"),
                 show_col_types = FALSE)

pal_sp <- c(Mouse = MV_MOUSE, Rat = MV_RAT, `Both species` = COL_ANALYSIS)

# Shared label-wrapping helper (mirrors the one in supp_fig1.R)
wrap_labs <- function(p, title_w = 28, sub_w = 36) {
  for (f in c("title", "subtitle")) {
    v <- p$labels[[f]]
    if (is.character(v) && length(v) == 1)
      p$labels[[f]] <- paste(strwrap(v, width = if (f == "title") title_w else sub_w),
                             collapse = "\n")
  }
  p + theme(plot.title.position = "plot",
            plot.title    = element_text(size = 9, face = "plain"),
            plot.subtitle = element_text(size = 7.5))
}

# =============================================================================
# Panel a: IV vs pooled Spearman rho (dumbbell / lollipop)
# old supp_fig4 panel C
# =============================================================================
# Route-restriction table has: Metric, Pooled-route fit, Intravenous fit
# Rows: Spearman rho, AUC, RMSE, Mouse rho

rat_iv_rho <- perf |> filter(Stratum == "Rat")   |> pull(`Spearman rho`)
rt_long <- rt |>
  filter(Metric %in% c("Spearman rho", "Mouse rho")) |>
  mutate(
    Stratum = case_when(
      Metric == "Spearman rho" ~ "Both species",
      Metric == "Mouse rho"    ~ "Mouse"
    ),
    Stratum = factor(Stratum, levels = c("Both species", "Mouse", "Rat"))
  ) |>
  bind_rows(
    tibble(Metric = "Rat rho", Stratum = factor("Rat", levels = c("Both species", "Mouse", "Rat")),
           `Pooled-route fit` = NA_real_, `Intravenous fit` = rat_iv_rho)
  ) |>
  pivot_longer(c(`Pooled-route fit`, `Intravenous fit`),
               names_to = "model", values_to = "rho") |>
  mutate(
    model   = factor(model, levels = c("Pooled-route fit", "Intravenous fit")),
    sp_col  = as.character(Stratum)
  )

p_a <- ggplot() +
  # segment connecting pooled to IV
  geom_segment(
    data = rt_long |>
      filter(!is.na(rho)) |>
      select(Stratum, model, rho) |>
      pivot_wider(names_from = model, values_from = rho) |>
      filter(!is.na(`Pooled-route fit`)),
    aes(x = `Pooled-route fit`, xend = `Intravenous fit`, y = Stratum, yend = Stratum),
    colour = "grey70", linewidth = 0.9
  ) +
  geom_point(
    data = rt_long |> filter(model == "Pooled-route fit", !is.na(rho)),
    aes(rho, Stratum), shape = 21, size = 3.5, fill = "grey75", colour = "grey50"
  ) +
  geom_point(
    data = rt_long |> filter(model == "Intravenous fit", !is.na(rho)),
    aes(rho, Stratum, fill = sp_col), shape = 23, size = 4,
    colour = "white"
  ) +
  scale_fill_manual(values = c("Both species" = COL_ANALYSIS,
                               "Mouse"        = MV_MOUSE,
                               "Rat"          = MV_RAT),
                    guide = "none") +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey70", linewidth = 0.5) +
  scale_x_continuous(limits = c(0.38, 0.70), breaks = seq(0.4, 0.7, 0.1)) +
  labs(x = "Spearman ρ",  y = NULL,
       title = "Route restriction: IV vs pooled ρ",
       subtitle = "Circle = pooled-route; diamond = IV-only") +
  theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))

# =============================================================================
# Panel b: per-species OOF performance (fold error + rho)
# old supp_fig4 panel D
# =============================================================================
SEED_COLS <- paste0("y_pred_cluster_PMID_seed", 0:2)
oof_iv <- oof |> filter(route_class == "iv")
# one median per seed, then the mean of the three: the convention every headline number uses.
# (Taking the median of the seed-AVERAGED prediction instead gave 2.13x and 2.29x, which agreed
# with neither the Discussion's 2.2x in rat nor the manuscript's headline convention.)
seed_fe <- oof_iv |>
  pivot_longer(all_of(SEED_COLS), names_to = "seed", values_to = "pred") |>
  mutate(fold_err = 10 ^ abs(pred - y_true))

perf_sp <- perf |>
  mutate(label = case_when(
    Stratum == "Both species" ~ "Both (n=213)",
    Stratum == "Rat"          ~ "Rat (n=179)",
    Stratum == "Mouse"        ~ "Mouse (n=34)"
  ) |> factor(levels = c("Mouse (n=34)", "Rat (n=179)", "Both (n=213)")))

fold_sp <- bind_rows(
  seed_fe |> group_by(seed) |> summarise(m = median(fold_err), .groups = "drop") |>
    summarise(med_fold = mean(m)) |>
    mutate(label = factor("Both (n=213)", levels = levels(perf_sp$label))),
  seed_fe |> group_by(seed, species) |> summarise(m = median(fold_err), .groups = "drop") |>
    group_by(species) |> summarise(med_fold = mean(m), .groups = "drop") |>
    mutate(label = factor(ifelse(species == "Rat",
                                 "Rat (n=179)", "Mouse (n=34)"),
                          levels = levels(perf_sp$label)))
)

col_sp <- c("Both (n=213)" = COL_ANALYSIS,
             "Rat (n=179)"  = MV_RAT,
             "Mouse (n=34)" = MV_MOUSE)

# b merges the old per-species and ranked-rho panels: the same three strata, two metrics, one
# facet each, horizontal bars. Fold error (left) is a median over the three seeds' medians; rho
# is the headline mean of three seeds.
perf_long <- bind_rows(
  fold_sp |> transmute(label, metric = "Median fold error", value = med_fold,
                       txt = sprintf("%.2f\u00d7", med_fold)),
  perf_sp |> transmute(label, metric = "Spearman \u03c1", value = `Spearman rho`,
                       txt = sprintf("%.3f", `Spearman rho`))
) |>
  mutate(metric = factor(metric, levels = c("Median fold error", "Spearman \u03c1")),
         label = factor(label, levels = rev(levels(perf_sp$label))))
perf_ref <- tibble(metric = factor(c("Median fold error", "Spearman \u03c1"),
                                   levels = levels(perf_long$metric)),
                   x = c(1, 0.5))

p_b <- ggplot(perf_long, aes(value, label, fill = label)) +
  geom_vline(data = perf_ref, aes(xintercept = x), linetype = "dashed", colour = "grey70",
             linewidth = 0.5, inherit.aes = FALSE) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_text(aes(label = txt), hjust = -0.12, size = 2.7) +
  facet_wrap(~metric, nrow = 1, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28)),
                     breaks = \(lim) if (lim[2] > 2) c(0, 2, 4) else c(0, 0.2, 0.4, 0.6)) +
  scale_fill_manual(values = col_sp, guide = "none") +
  labs(x = NULL, y = NULL, title = "Performance by species stratum",
       subtitle = "Fold error: median, mean of 3 seeds; \u03c1: headline mean of 3 seeds") +
  theme(strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_blank())

# =============================================================================
# Panel d: between-source label agreement (cross-source pairs)
# =============================================================================
# The pairs of the Methods' rule: every sequence-and-species combination that more than one source
# reports, each source reduced to its mean log10 half-life, every unordered pair of sources. The panel
# used to draw only the records flagged as repeat candidates (6 pairs, r = 0.94).
cs <- read_csv(file.path(LIVE_TABS, "table60_muriverse_halflife_cross_source_pairs_iv.csv"),
               show_col_types = FALSE) |>
  mutate(sp = str_to_title(species), combo = paste(sequence, sp))
dominant <- cs |> count(combo, sort = TRUE) |> slice(1) |> pull(combo)
cs <- cs |> mutate(dominant = combo == dominant)
r_cs   <- cor(cs$log10_h_a, cs$log10_h_b, method = "pearson")
rho_cs <- cor(cs$log10_h_a, cs$log10_h_b, method = "spearman")
lim_cs <- range(c(cs$log10_h_a, cs$log10_h_b), na.rm = TRUE) + c(-0.15, 0.15)

p_d <- ggplot(cs, aes(log10_h_a, log10_h_b, colour = sp, shape = dominant)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey70", linewidth = 0.6) +
  geom_point(size = 2.4, alpha = 0.75) +
  scale_colour_manual(values = c(Mouse = MV_MOUSE, Rat = MV_RAT), name = "Species") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), name = NULL,
                     labels = c(`FALSE` = "other peptides",
                                `TRUE`  = sprintf("one peptide, %d pairs", sum(cs$dominant)))) +
  coord_equal(xlim = lim_cs, ylim = lim_cs) +
  annotate("text", x = lim_cs[1] + 0.05, y = lim_cs[2] - 0.05,
           label = sprintf("n = %d pairs, %d peptide-species\nr = %.2f, ρ = %.2f\nmedian |Δ| = %.2f log₁₀",
                           nrow(cs), n_distinct(cs$combo), r_cs, rho_cs, median(cs$abs_diff_log10)),
           hjust = 0, vjust = 1, size = 2.5, colour = "grey30") +
  scale_x_loghl(breaks = c(1 / 60, 5 / 60, 0.25, 1, 6, 24)) + scale_y_loghl(breaks = c(1 / 60, 5 / 60, 0.25, 1, 6, 24)) +
  labs(x = "Half-life, source A",
       y = "Half-life, source B",
       title = "Between-source label agreement") +
  theme(legend.position  = "bottom",
        legend.box       = "vertical",
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(3, "mm"))

# =============================================================================
# Panel c: OOF predicted vs observed (213-record IV set)
# =============================================================================
# Points are seed 0, one real out-of-fold pass, as in Figure 3c and e (the convention of fig3.R). The numbers are the HEADLINE metrics, the mean of the three
# per-seed scores (analysis/metrics_iv.csv). The panel used to print rho 0.609, r 0.594 and RMSE 0.659,
# which are the scores of the seed-AVERAGED prediction - better on every metric, and not the published ones.
oof_iv <- oof |> filter(route_class == "iv") |> mutate(pred_seed2 = y_pred_cluster_PMID_seed0)  # name kept; the column is seed 0
met <- read_csv(file.path(DIR_REV, "metrics_iv.csv"), show_col_types = FALSE)
mval <- function(m) met$value[met$metric == m]

lim_oof <- range(c(oof_iv$pred_seed2, oof_iv$y_true)) + c(-0.15, 0.15)

p_c <- ggplot(oof_iv, aes(pred_seed2, y_true, colour = species)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey70", linewidth = 0.6) +
  geom_point(size = 1.3, alpha = 0.55) +
  scale_colour_manual(values = c(Mouse = MV_MOUSE, Rat = MV_RAT), name = "Species") +
  coord_equal(xlim = lim_oof, ylim = lim_oof) +
  annotate("text", x = lim_oof[1] + 0.05, y = lim_oof[2] - 0.05,
           label = sprintf("ρ = %.3f\nr = %.3f\nRMSE = %.3f log₁₀ h",
                           mval("Spearman rho"), mval("Pearson r"), mval("RMSE (log10 h)")),
           hjust = 0, vjust = 1, size = 2.5, colour = "grey30") +
  scale_x_loghl(breaks = c(1 / 60, 5 / 60, 0.25, 1, 6, 24)) + scale_y_loghl(breaks = c(1 / 60, 5 / 60, 0.25, 1, 6, 24)) +
  labs(x = "Predicted half-life",
       y = "Observed half-life",
       title = "Out-of-fold: predicted vs observed",
       subtitle = "Points: seed 0. Metrics: mean of 3 seeds") +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(3, "mm"))

# =============================================================================
# Assembly
# =============================================================================
LAB_S2 <- list(size = 8, face = "bold")

row1 <- ggarrange(wrap_labs(p_a), wrap_labs(p_b), ncol = 2, widths = c(0.8, 1.2),
                  labels = c("a", "b"), font.label = LAB_S2)
row2 <- ggarrange(wrap_labs(p_c), wrap_labs(p_d), ncol = 2,
                  labels = c("c", "d"), font.label = LAB_S2)

supp_fig2 <- ggarrange(row1, row2, nrow = 2, heights = c(1, 1.05))

save_fig(supp_fig2, "supp_fig2", width = 200, height = 215)
