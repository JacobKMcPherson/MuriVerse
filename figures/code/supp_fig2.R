# Supplementary Figure 2 (4 panels, 2 x 2): performance and label reproducibility.
#
# Row 1: a  IV vs pooled route-restriction delta-rho (lollipop/dumbbell)
#         b  Per-species performance: median fold error and Spearman rho, one facet each,
#            horizontal bars over the same three strata
# Row 2: c  OOF predicted vs observed (213-record IV set)
#         d  Between-source label agreement (cross-source pairs scatter)
#
# Data: MuriVerse/rodent/strict/ and figures/tabs/pub-tab__{2,3} for panels a and b (plus the
# frozen results/supp_fig2/ tables), and
# figures/tabs/table60_muriverse_halflife_cross_source_pairs_iv.csv
# for panel d. Panel c plots seed 0, like Figure 3c, and prints the HEADLINE metrics (mean of
# three seeds); panel b reports the median fold error the same way (mean of the three seeds' medians).
if (!exists("ROOT")) {
  source(file.path(system("git rev-parse --show-toplevel", intern = TRUE),
                   "figures/code/00_setup.R"))
}
# Floor for every piece of text in this figure, in pt and in geom_text's mm.
SZ_S2_MIN   <- 8
SZ_S2_ANNOT <- SZ_S2_MIN / ggplot2::.pt

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
            plot.subtitle = element_text(size = 8))
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
  geom_text(aes(label = txt), hjust = -0.12, size = SZ_S2_ANNOT) +
  facet_wrap(~metric, nrow = 1, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.4)),
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
oof_iv <- oof |> filter(route_class == "iv") |> mutate(pred_seed2 = y_pred_cluster_PMID_seed0)  # name kept; the column is seed 0
# One range for panels c and d, so their square plotting areas come out the same size.
lim_cd <- range(c(lim_cs, range(c(oof_iv$pred_seed2, oof_iv$y_true)) + c(-0.15, 0.15)))
CD_BREAKS_X <- c(1 / 60, 0.25, 6)
CD_BREAKS_Y <- c(1 / 60, 0.25, 1, 6)
# Panels c and d: the plot, then a fixed-width side column holding the legend over the statistics,
# so both square plotting areas get the same width whatever their legends contain.
side_panel <- function(p, stats) {
  leg <- ggpubr::as_ggplot(ggpubr::get_legend(p))
  txt <- ggpubr::as_ggplot(ggpubr::text_grob(stats, size = SZ_S2_MIN, color = "grey30", hjust = 0, x = 0.02,
                                             lineheight = 0.95))
  ggarrange(p + theme(legend.position = "none"), ggarrange(NULL, leg, txt, ncol = 1, heights = c(0.3, 1, 1)),
            ncol = 2, widths = c(1, 0.5))
}
SIDE_THEME <- theme(legend.position = "right", legend.justification = c(0, 1),
                    legend.box = "vertical", legend.box.just = "left",
                    legend.title = element_text(size = 8), legend.text = element_text(size = 8),
                    legend.key.size = unit(3, "mm"), legend.background = element_blank(),
                    axis.text.x = element_text(angle = 0, hjust = 0.5),
                    plot.margin = margin(14, 6, 4, 8))

p_d <- ggplot(cs, aes(log10_h_a, log10_h_b, colour = sp, shape = dominant)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey70", linewidth = 0.6) +
  geom_point(size = 2.4, alpha = 0.75) +
  scale_colour_manual(values = c(Mouse = MV_MOUSE, Rat = MV_RAT), name = "Species") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), name = NULL,
                     labels = c(`FALSE` = "other peptides",
                                `TRUE`  = sprintf("one peptide\n(%d pairs)", sum(cs$dominant)))) +
  coord_equal(xlim = lim_cd, ylim = lim_cd) +
  labs(x = "Half-life, source A",
       y = "Half-life, source B",
       title = "Between-source label agreement",
       subtitle = "Same peptide and species, two sources") +
  scale_x_loghl(breaks = CD_BREAKS_X) + scale_y_loghl(breaks = CD_BREAKS_Y) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 2)),
         shape = guide_legend(order = 2, ncol = 1, override.aes = list(colour = "grey35", size = 2))) +
  SIDE_THEME

# =============================================================================
# Panel c: OOF predicted vs observed (213-record IV set)
# =============================================================================
# Points are seed 0, one real out-of-fold pass, as in Figure 3c and e (the convention of fig3.R). The numbers are the HEADLINE metrics, the mean of the three
# per-seed scores (analysis/metrics_iv.csv). The panel used to print rho 0.609, r 0.594 and RMSE 0.659,
# which are the scores of the seed-AVERAGED prediction - better on every metric, and not the published ones.
met <- read_csv(file.path(DIR_REV, "metrics_iv.csv"), show_col_types = FALSE)
mval <- function(m) met$value[met$metric == m]

p_c <- ggplot(oof_iv, aes(pred_seed2, y_true, colour = species)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey70", linewidth = 0.6) +
  geom_point(size = 1.3, alpha = 0.55) +
  scale_colour_manual(values = c(Mouse = MV_MOUSE, Rat = MV_RAT), name = "Species") +
  coord_equal(xlim = lim_cd, ylim = lim_cd) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "Predicted half-life",
       y = "Observed half-life",
       title = "Out-of-fold: predicted vs observed",
       subtitle = "Points: seed 0; metrics: mean of 3 seeds") +
  scale_x_loghl(breaks = CD_BREAKS_X) + scale_y_loghl(breaks = CD_BREAKS_Y) +
  SIDE_THEME

# =============================================================================
# Panels e-h: label-certainty subsets and ensemble members (revision, round 1)
# =============================================================================
# e, f: the ensemble refitted on each restricted subset of the 213 intravenous records under the
# same grouped five-fold protocol (figures/tabs/supp_table8_cohort_sensitivity.csv), scored with
# Spearman rho (e) and mean absolute error (f), each beside the species-median baseline refitted in
# the same folds. The grey bars in e are the 5th-95th percentile of 30 random subsets of the
# body-verified size (the two bottom rows), so a drop that a smaller training set alone would
# produce is visible against the estimate.
# g, h: the shipped ensemble against its three members and the species-median constant on identical
# folds (supp_table7_ensemble_member_comparison.csv): the error with its 95% group-bootstrap interval
# (g) and the paired ensemble-minus-comparator difference in MAE (h).
SUB_LEVELS <- c("All intravenous records (shipped cohort)",
                "Exact, unambiguous half-life labels only",
                "Records whose source contradicts the label removed",
                "Contradicted records removed, three corrected values restored",
                "Verified in full primary-source text or abstract",
                "Verified in full primary-source text only",
                "Verified in full primary-source text, largest patent excluded")
SUB_SHORT <- c("All records", "Exact labels", "No contradiction", "+ 3 corrected",
               "Full text or abstract", "Full text only", "Full text, no top patent")
sens <- read_csv(file.path(LIVE_TABS, "supp_table8_cohort_sensitivity.csv"), show_col_types = FALSE) |>
  mutate(label = factor(sprintf("%s (n = %d)", SUB_SHORT[match(scenario, SUB_LEVELS)], n),
                        levels = rev(sprintf("%s (n = %d)", SUB_SHORT,
                                             sapply(SUB_LEVELS, \(x) first(n[scenario == x]))))),
         estimator = factor(estimator, levels = c("species median", "ensemble")))
size_ctl <- read_csv(file.path(ROOT, "results", "revision_analysis", "verification_size_control.csv"),
                     show_col_types = FALSE) |>
  group_by(matched_to) |>
  summarise(n = first(n), lo = quantile(ensemble_spearman, 0.05), hi = quantile(ensemble_spearman, 0.95),
            mean = mean(ensemble_spearman), .groups = "drop") |>
  mutate(label = factor(sprintf("%s (n = %d)",
                                SUB_SHORT[match(c("verified_body" = "Verified in full primary-source text only",
                                                  "verified_body_excluding_largest_patent" =
                                                    "Verified in full primary-source text, largest patent excluded")[matched_to],
                                                SUB_LEVELS)], n),
                        levels = levels(sens$label)))
SIZE_KEY <- "same-size random subsets, 5th-95th percentile"
# Grey is not used for any estimate: it would collide with the interval lines and the grid.
pal_est <- c(ensemble = "#332288", `species median` = "#CC6677", setNames("#DDCC77", SIZE_KEY))
DODGE_S <- position_dodge(width = 0.6)

p_e <- ggplot(sens, aes(spearman, label, colour = estimator)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_linerange(data = size_ctl, aes(xmin = lo, xmax = hi, y = label, colour = SIZE_KEY),
                 inherit.aes = FALSE, linewidth = 3.2) +
  geom_point(data = size_ctl, aes(x = mean, y = label), inherit.aes = FALSE, shape = 124,
             colour = "grey20", size = 3) +
  geom_linerange(aes(xmin = spearman_ci_low, xmax = spearman_ci_high), position = DODGE_S, linewidth = 0.5) +
  geom_point(position = DODGE_S, size = 1.9) +
  scale_colour_manual(values = pal_est, name = NULL, breaks = c("ensemble", "species median", SIZE_KEY)) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = c(0.5, 0.5, 3), shape = c(16, 16, NA)))) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
  labs(x = "Spearman \u03c1", y = NULL,
       title = "Label certainty",
       subtitle = "Lines: 95% CI; tick: random-subset mean") +
  theme(legend.position = "bottom", panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        axis.text.y = element_text(size = 8), axis.text.x = element_text(angle = 0, hjust = 0.5))

p_f <- ggplot(sens, aes(mae_log10, label, colour = estimator)) +
  geom_linerange(aes(xmin = mae_ci_low, xmax = mae_ci_high), position = DODGE_S, linewidth = 0.5) +
  geom_point(position = DODGE_S, size = 1.9) +
  scale_colour_manual(values = pal_est, name = NULL, breaks = c("ensemble", "species median")) +
  scale_x_continuous(breaks = c(0.4, 0.6, 0.8, 1.0)) +
  labs(x = "MAE, log\u2081\u2080 h", y = NULL, title = "Absolute error, same subsets",
       subtitle = "Species median refitted per fold") +
  theme(legend.position = "bottom", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

MEM_LEVELS <- c("ensemble", "ridge", "svr", "xgboost", "species_median")
MEM_LABELS <- c("Ensemble (shipped)", "Ridge regression", "SVR (RBF kernel)", "XGBoost", "Species median")
mem <- read_csv(file.path(LIVE_TABS, "supp_table7_ensemble_member_comparison.csv"), show_col_types = FALSE) |>
  mutate(label = factor(MEM_LABELS[match(estimator, MEM_LEVELS)], levels = rev(MEM_LABELS)))
p_g <- ggplot(mem, aes(mae_log10, label)) +
  geom_linerange(aes(xmin = mae_ci_low, xmax = mae_ci_high), colour = "#332288", linewidth = 0.5) +
  geom_point(colour = "#332288", size = 2) +
  geom_text(aes(x = mae_ci_high, label = sprintf("%.3f", mae_log10)), hjust = -0.15, size = SZ_S2_ANNOT) +
  scale_x_continuous(breaks = c(0.4, 0.6, 0.8, 1.0), expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = "MAE, log\u2081\u2080 h (95% CI)", y = NULL, title = "Ensemble versus its members",
       subtitle = "Identical folds; 213 records, 55 groups") +
  theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

p_h <- ggplot(filter(mem, estimator != "ensemble"), aes(ensemble_minus_comparator_mae, label)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
  geom_linerange(aes(xmin = difference_ci_low, xmax = difference_ci_high), colour = "#332288", linewidth = 0.5) +
  geom_point(colour = "#332288", size = 2) +
  scale_y_discrete(drop = FALSE, limits = levels(mem$label)) +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1)) +
  labs(x = "MAE difference (95% CI)", y = NULL, title = "Ensemble minus comparator",
       subtitle = "Left of zero: ensemble has lower error") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

# =============================================================================
# Assembly
# =============================================================================
LAB_S2 <- list(size = 10, face = "bold")

row1 <- ggarrange(wrap_labs(p_a), wrap_labs(p_b), ncol = 2, widths = c(0.8, 1.2),
                  labels = c("a", "b"), font.label = LAB_S2)
p_c_side <- side_panel(wrap_labs(p_c, 44, 50),
                       sprintf("\u03c1 = %.3f\nr = %.3f\nRMSE = %.3f\nlog\u2081\u2080 h", mval("Spearman rho"),
                               mval("Pearson r"), mval("RMSE (log10 h)")))
p_d_side <- side_panel(wrap_labs(p_d, 44, 50),
                       sprintf("%d pairs,\n%d peptide-species\nr = %.2f, \u03c1 = %.2f\nmedian |\u0394| = %.2f\nlog\u2081\u2080 h",
                               nrow(cs), n_distinct(cs$combo), r_cs, rho_cs, median(cs$abs_diff_log10)))
row2 <- ggarrange(p_c_side, p_d_side, ncol = 2,
                  labels = c("c", "d"), font.label = LAB_S2)
row3 <- ggarrange(wrap_labs(p_e, 60, 64), wrap_labs(p_f, 40, 40), ncol = 2, widths = c(1.3, 1), align = "h",
                  labels = c("e", "f"), font.label = LAB_S2, common.legend = TRUE,
                  legend = "bottom")
row4 <- ggarrange(wrap_labs(p_g, 60, 64), wrap_labs(p_h, 40, 40), ncol = 2, widths = c(1.45, 1), align = "h",
                  labels = c("g", "h"), font.label = LAB_S2)

supp_fig2 <- ggarrange(row1, row2, row3, row4, nrow = 4, heights = c(0.85, 1.05, 1.0, 0.62))

save_fig(supp_fig2, "supp_fig2", width = 180, height = 250)
