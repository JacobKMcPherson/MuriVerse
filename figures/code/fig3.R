# Figure 3 evaluates the shipped half-life model, which is fitted on the 213 INTRAVENOUS
# records of the 320-row cohort. It opens with the measurement that justifies that
# restriction - cross-species and cross-route transfer over the four species x route cells
# (panels a and b, the only panels drawn over all 320 rows) - then plots held-out predictions
# for every intravenous row (c-f) and compares all feature-block combinations (g). Panel h is
# the paired flip-flop evidence (all 320 rows) and panel i the saved split's label shift,
# which no performance number uses. Mouse and rat are descriptive strata of that one model;
# Figure 2 gives the input cohort.
#
# ggarrange order (the layout is the truth if a caption disagrees):
#   a p_e, b p_b_route | c p_f, d p_d, e p_g | f p_i, g p_j, h p_c, i p_b
# Built but NOT drawn: p_h (interval coverage), p_a (80/20 split counts).
if (!exists("ROOT")) {
  source(file.path(system("git rev-parse --show-toplevel", intern = TRUE), "figures/code/00_setup.R"))
}
if (!exists("hl_preds")) source(file.path(DIR_CODE, "05_data_halflife_model.R"))

# Mouse and rat take the species hues Figure 2 uses for the same two cohorts (00_setup.R).
# The both-species curve, which is the model's actual score, is the darker neutral the eye
# reads first; the two species curves are the decomposition of it. Labelled "Both species"
# rather than "Pooled" because "pooled" elsewhere in this figure means route-pooled.
pal_stratum <- c(`Both species` = "grey25",
                 Mouse  = MV_MOUSE,
                 Rat    = MV_RAT)

# Routes, ordered by how much of the cohort each carries (iv 213, sc 73, ip 22, im 12) so
# the key reads in the order the reader should weigh them. IV keeps the dark model colour
# used for the pooled model elsewhere in the figure because it supplies two thirds of the
# rows and effectively sets the headline score.
# Panels d and h both plot log10 half-life on both axes against an identity line, so they
# share one square range. That makes the two plot boxes render at the same size - coord_equal
# otherwise takes its aspect from each panel's own data - and lets the reader compare the
# spread in one against the other directly.
CD_LIM <- c(-1.9, 1.9)

ROUTE_LEVELS <- c("iv", "sc", "ip", "im")
pal_route <- c(iv = "#332288", sc = "#CC6677", ip = "#DDCC77", im = "#AA4499")

# The figure is drawn at the 180 mm journal column width (save_fig's default), so every panel
# is about 90 mm across and type has to be sized for that, not for a screen-scale canvas.
# Sizes are named once here rather than written into each panel's theme() so the whole figure
# retunes from one place and no panel can drift to an illegible 5 pt on its own. Floor is
# ~5.5 pt, which is the usual journal minimum at final size.
SZ_TITLE   <- 7.5
SZ_SUB     <- 6.2
SZ_AXTITLE <- 6.8
SZ_AXTEXT  <- 6.0
SZ_LEGEND  <- 5.8
SZ_STRIP   <- 6.2

# Applied to every panel, so the per-panel theme() calls below carry only what is genuinely
# specific to that panel.
theme_panel <- function(...) {
  theme(plot.title = element_text(size = SZ_TITLE, face = "plain"),
        plot.subtitle = element_text(size = SZ_SUB),
        axis.title = element_text(size = SZ_AXTITLE),
        axis.text = element_text(size = SZ_AXTEXT),
        strip.text = element_text(size = SZ_STRIP),
        legend.text = element_text(size = SZ_LEGEND),
        legend.title = element_text(size = SZ_LEGEND),
        plot.margin = margin(8, 4, 2, 4),
        ...)
}

# a: the split itself, in counts. Records, homology clusters and source PMIDs for the saved
# cluster-preserving 80/20 split, with the shared-group counts called out beside the clusters
# and sources rows. Zero shared clusters is the design property of the split: whole homology
# families remain together, so any peptide in validation shares no family with training. The
# four shared sources are the residual - a cluster can contain peptides from multiple PMIDs,
# so holding a cluster out does not fully disjoint sources. This panel sets up panel b, which
# reads out the label-range consequence of that cluster-preserving design.
split_stats <- bind_rows(
  tibble(stat = "records",
         split = "train", n = sum(hl_meta$split == "train")),
  tibble(stat = "records",
         split = "val",   n = sum(hl_meta$split == "val")),
  tibble(stat = "homology clusters",
         split = "train", n = length(unique(hl_meta$cluster[hl_meta$split == "train"]))),
  tibble(stat = "homology clusters",
         split = "val",   n = length(unique(hl_meta$cluster[hl_meta$split == "val"]))),
  tibble(stat = "source PMIDs",
         split = "train", n = length(unique(hl_meta$PMID[hl_meta$split == "train"]))),
  tibble(stat = "source PMIDs",
         split = "val",   n = length(unique(hl_meta$PMID[hl_meta$split == "val"])))
) |>
  # A discrete y scale draws its first level at the BOTTOM, so the levels run in reverse of
  # the intended reading order: records, source PMIDs, homology clusters top to bottom.
  # Likewise `split` is ordered val-then-train so the dodge puts training ABOVE validation
  # in each group; scale_fill_manual restores training-first in the key.
  mutate(stat = factor(stat, levels = c("homology clusters", "source PMIDs", "records")),
         split = factor(split, levels = c("val", "train")))

shared_clusters <- length(intersect(
  unique(hl_meta$cluster[hl_meta$split == "train"]),
  unique(hl_meta$cluster[hl_meta$split == "val"])
))
shared_sources <- length(intersect(
  unique(hl_meta$PMID[hl_meta$split == "train"]),
  unique(hl_meta$PMID[hl_meta$split == "val"])
))
overlap_labels <- tibble(
  stat = factor(c("homology clusters", "source PMIDs"), levels = levels(split_stats$stat)),
  label = c(sprintf("%d shared", shared_clusters), sprintf("%d shared", shared_sources))
)

p_a <- ggplot(split_stats, aes(x = n, y = stat, fill = split)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(aes(label = n), position = position_dodge(width = 0.72),
            hjust = 1.2, size = 2.0) +
  geom_text(data = overlap_labels, aes(x = Inf, y = stat, label = label),
            hjust = -0.05, size = 2.0, colour = "grey30", inherit.aes = FALSE) +
  scale_fill_manual(values = c(train = "#332288", val = "#DDCC77"), name = NULL,
                    breaks = c("train", "val"),
                    labels = c(train = "training", val = "validation")) +
  # Headroom for the count at the end of the longest bar ("257"); at a quarter of the page the
  # old 22% let it run off the facet edge.
  # Flipped: bars grow leftward from zero at the right edge, with the row labels beside zero.
  scale_x_reverse(expand = expansion(mult = c(0.4, 0))) +
  scale_y_discrete(position = "right") +
  labs(x = "Count", y = NULL) +
  # The key stacks in one column: side by side, "training" and "validation" together are wider
  # than this panel and the second entry was cut off.
  guides(fill = guide_legend(nrow = 1)) +
  theme_panel() +
  theme(legend.position = "bottom", legend.location = "plot", legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = 5.4),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"))

# d: the accuracy profile of the same predictions panel c plots - what fraction of held-out
# records land within a given fold of their observed half-life.
#
# This slot used to hold a quantile-binned calibration plot, and it was the weakest panel in
# the figure for a structural reason: binning at ~12 records per bin left RAT with ten points
# and MOUSE with two, so a panel whose job was to compare the two animals' absolute scale
# drew one of them with two dots. Twelve bin-means also hide the spread they are built from -
# the binned version showed mouse and rat biases of -0.09 and -0.04 and so looked like two
# animals on the same footing, while the records behind those means differ by nearly a factor
# of two in median fold error (mouse 4.1x, rat 2.2x). The curve below uses every record.
#
# Read it as: pick a fold on x, read off what share of predictions are at least that close.
# The predeclared reference is 3-fold, marked with the heavier line, because a dosing-interval
# decision is the use this model is being judged for and 2-fold is beyond it on this cohort.
# Each curve is the MEAN over the three seeds of that seed's empirical distribution, so the
# quoted percentages are the same quantity as the line they label and no seed caveat applies
# (unlike panels c, e and f, which draw one pass).
FOLD_GRID <- 10^seq(0, log10(32), length.out = 300)
FOLD_REF <- 3

fold_profile <- function(d) {
  m <- vapply(SEED_COLS, \(cl) {
    fe <- 10^abs(d$y_true - d[[cl]])
    vapply(FOLD_GRID, \(t) mean(fe <= t), numeric(1))
  }, numeric(length(FOLD_GRID)))
  tibble(fold = FOLD_GRID, share = rowMeans(m))
}

fold_strata <- bind_rows(
  mutate(fold_profile(hl_oof), stratum = "Both species"),
  mutate(fold_profile(filter(hl_oof, species == "mouse")), stratum = "Mouse"),
  mutate(fold_profile(filter(hl_oof, species == "rat")),   stratum = "Rat")
) |>
  mutate(stratum = factor(stratum, levels = c("Both species", "Mouse", "Rat")))

# The number each curve is labelled with: its own share at the reference fold, taken from the
# same averaged profile that draws the line.
fold_at_ref <- fold_strata |>
  group_by(stratum) |>
  summarise(share = approx(fold, share, xout = FOLD_REF)$y, .groups = "drop") |>
  mutate(lab = sprintf("%s %.0f%%", stratum, 100 * share))

# Median fold error per stratum, three-seed mean - the other half of the same picture, and
# the quantity the error-stratification panel (f) breaks down.
fold_median <- function(d) mean(vapply(SEED_COLS,
                                       \(cl) median(10^abs(d$y_true - d[[cl]])), numeric(1)))

p_d <- ggplot(fold_strata, aes(fold, share, colour = stratum)) +
  geom_vline(xintercept = 2, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = FOLD_REF, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = pal_stratum, name = NULL,
                      labels = setNames(fold_at_ref$lab, fold_at_ref$stratum)) +
  scale_x_log10(breaks = c(1, 2, 3, 5, 10, 30),
                labels = c("1\u00d7", "2\u00d7", "3\u00d7", "5\u00d7", "10\u00d7", "30\u00d7"),
                expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1), expand = expansion(mult = 0.02)) +
  labs(x = "Within this fold of observed", y = "Share of held-out records",
       title = sprintf("Accuracy: %.0f%% within %g-fold",
                       100 * fold_at_ref$share[fold_at_ref$stratum == "Both species"], FOLD_REF)) +
  # Bottom-right, in the light grey the other panels use for provenance: the median fold
  # error the curves imply, and the signed bias the binned calibration plot used to carry.
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("median %.1f\u00d7 mouse, %.1f\u00d7 rat\nbias %+.2f, %+.2f log10",
                           fold_median(filter(hl_oof, species == "mouse")),
                           fold_median(filter(hl_oof, species == "rat")),
                           mean(vapply(SEED_COLS, \(cl) mean(hl_oof[[cl]][hl_oof$species == "mouse"] -
                                       hl_oof$y_true[hl_oof$species == "mouse"]), numeric(1))),
                           mean(vapply(SEED_COLS, \(cl) mean(hl_oof[[cl]][hl_oof$species == "rat"] -
                                       hl_oof$y_true[hl_oof$species == "rat"]), numeric(1)))),
           hjust = 1.04, vjust = -0.3, size = 1.45, colour = "grey45", lineheight = 0.95) +
  guides(colour = guide_legend(ncol = 1)) +
  theme_panel() +
  theme(plot.title = element_text(size = 6.2),
        legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1), legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = 5.0), legend.background = element_blank(),
        legend.margin = margin(0, 0, 0, 0))

# f: observed vs predicted on the held-out clusters and sources. One model produced every
# point; the facets cut its validation set by species, so the two panels are the same
# regressor judged on the two animals it now has a column for - a mouse and a rat row are
# no longer the same input, though the offset it applies between them is small.
#
# The POINTS are seed 0 - one real out-of-fold pass, which is what a deployed model is - but
# the QUOTED rho is the three-seed mean, because that is the figure's headline quantity and a
# single seed is not it. The distinction is not cosmetic on this cohort: at seed 0 mouse
# scores 0.68 against rat's 0.65, which reads as "mouse is the better-ordered animal", while
# across the three seeds mouse runs 0.68/0.51/0.50 and rat 0.65/0.60/0.55 - mouse is BELOW
# rat in the other two seeds and the seed-0 reading overstates it.
# Mouse is the volatile stratum because it is 34 records; the seed range is printed beside
# the mean so the reader can size that volatility instead of inferring stability from one
# number. An earlier version quoted seed 0 alone with no indication that it did.
seed_rho <- function(d) vapply(SEED_COLS, \(cl) cor(d$y_true, d[[cl]], method = "spearman"),
                               numeric(1))
species_rho <- hl_oof |>
  mutate(species = str_to_title(species)) |>
  group_by(species) |>
  group_modify(\(d, k) { v <- seed_rho(d)
                         tibble(n = nrow(d), mean_rho = mean(v),
                                lo = min(v), hi = max(v), seed0 = v[[1]]) }) |>
  ungroup()

# Kept short: at a third of the page the longer "[... over 3 seeds]" wording was clipped by
# the facet edge, so what the bracket MEANS is stated once in the caption instead of twice
# in the strips.
# One p-value per facet, in the bottom-right corner. Reported as the WORST (largest) of the
# three seeds, so the number holds whichever pass a reader is looking at and needs no seed
# caveat of its own - unlike rho, where the mean is the quantity of interest.
species_p <- hl_oof |>
  mutate(species = str_to_title(species)) |>
  group_by(species) |>
  group_modify(\(d, k) tibble(p = max(vapply(SEED_COLS,
      \(cl) cor.test(d$y_true, d[[cl]], method = "spearman", exact = FALSE)$p.value,
      numeric(1))))) |>
  ungroup()

label_species <- species_rho |>
  mutate(lab = sprintf("%s (n = %d)\n\u03c1 = %.2f (%.2f\u2013%.2f)",
                       species, n, mean_rho, lo, hi)) |>
  select(species, lab) |> deframe()

p_f <- hl_preds |>
  mutate(species = factor(label_species[species], levels = unname(label_species))) |>
  ggplot(aes(y_true, y_pred, colour = species)) +
  geom_abline(linetype = "dashed", colour = "grey60") +
  geom_point(alpha = 0.85, size = 0.9) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.6) +
  scale_colour_manual(values = unname(pal_stratum[c("Mouse", "Rat")]), guide = "none") +
  geom_text(data = species_p |>
              mutate(species = factor(label_species[species], levels = unname(label_species)),
                     lab = if_else(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))),
            aes(x = Inf, y = -Inf, label = lab), inherit.aes = FALSE,
            hjust = 1.08, vjust = -0.5, size = 1.5, colour = "grey45") +
  # Mouse beside rat, each stratum's label above its facet.
  facet_wrap(~species, nrow = 1) +
  scale_x_loghl(breaks = c(1 / 60, 5 / 60, 1, 24)) + scale_y_loghl(breaks = c(1 / 60, 5 / 60, 1, 24)) +
  labs(x = "Observed half-life",
       y = "Predicted half-life",
       caption = "points: seed 0  \u00b7  \u03c1: 3-seed mean (range)") +
  theme_panel() +
  theme(strip.text = element_text(size = 5.0, lineheight = 0.95),
        plot.caption = element_text(size = 4.4, colour = "grey45", hjust = 0))

# g: the same held-out predictions read as a ranking problem at the predeclared 6 h target.
# The model is a regressor, so the label is dichotomised at HL_POS_THRESHOLD_H and the
# continuous prediction serves as the score.
#
# POOLED OVER SPECIES, deliberately. The earlier version split this curve by species, which
# the intravenous cohort cannot support: 213 IV records hold only 17 at or above 6 h, 8 mouse
# and 9 rat, so a per-species ROC would rest on single-digit positive counts and its AUC
# would be an artefact of which two or three records happened to land where. The positive
# count is printed on the panel rather than left to the caption, because it is the limiting
# quantity here and a reader should not have to look it up to size the claim.
#
# Note that `hl_roc`/`hl_clf_metrics` from 05_data_halflife_model.R are built at the 4 h
# SCREEN cutoff, so the 6 h curve is derived here from hl_preds rather than reused. AUC is
# computed from the same seed-0 points that draw the curve, so the legend number cannot
# drift from the line beside it.
roc_6h <- hl_preds |>
  mutate(positive_6h = half_life_h >= HL_POS_THRESHOLD_H)

roc_6h_curve <- roc_curve(roc_6h$y_pred, roc_6h$positive_6h)
n_pos_6h <- sum(roc_6h$positive_6h)
n_pos_by_species <- roc_6h |> filter(positive_6h) |> count(species)

# As in panel c, the CURVE is seed 0 and the quoted AUC is the three-seed mean, so the number
# on this panel is the same quantity the figure quotes elsewhere. Seed 0 is the most
# favourable of the three here (0.84 / 0.84 / 0.77), so quoting it alone would overstate the
# screen by about 0.03.
auc_6h_seeds <- vapply(SEED_COLS,
                       \(cl) rank_auc(hl_oof[[cl]], hl_oof$half_life_h >= HL_POS_THRESHOLD_H),
                       numeric(1))
roc_6h_auc <- mean(auc_6h_seeds)

p_g <- ggplot(roc_6h_curve, aes(fpr, tpr)) +
  geom_abline(linetype = "dashed", colour = "grey60") +
  geom_step(linewidth = 0.55, colour = "#332288") +
  coord_equal() +
  labs(x = "False positive rate", y = "True positive rate",
       title = sprintf("Ranking: half-life \u2265 %g h", HL_POS_THRESHOLD_H),
       subtitle = sprintf("AUC %.2f [%.2f\u2013%.2f over 3 seeds]",
                          roc_6h_auc, min(auc_6h_seeds), max(auc_6h_seeds))) +
  # The provenance and the positive count go to the empty bottom-right corner of a ROC plot
  # in light grey. They are what sizes the claim rather than the claim itself, so the
  # subtitle carries the AUC alone.
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("curve: seed 0\n%d of %d positive (mouse %d, rat %d)",
                           n_pos_6h, nrow(roc_6h),
                           sum(n_pos_by_species$n[n_pos_by_species$species == "Mouse"]),
                           sum(n_pos_by_species$n[n_pos_by_species$species == "Rat"])),
           hjust = 1.04, vjust = -0.3, size = 1.45, colour = "grey45", lineheight = 0.95) +
  theme_panel() +
  theme(plot.title = element_text(size = 6.2, face = "plain"),
        plot.subtitle = element_text(size = 5.0))

# h: leave-group-out residual conformal bands. For every connected cluster/source group, the
# interval radius is the quantile of absolute residuals from the other 74 groups. Thus the
# evaluated group's labels never set its own interval width. These are empirical uncertainty
# bands around the existing grouped out-of-fold predictions, not a substitute for a prospective
# PK study.
grouped_preds <- hl_preds |>
  mutate(cv_group = row_number())
repeat {
  before <- grouped_preds$cv_group
  grouped_preds <- grouped_preds |>
    group_by(cluster) |>
    mutate(cv_group = min(cv_group)) |>
    ungroup() |>
    group_by(PMID) |>
    mutate(cv_group = min(cv_group)) |>
    ungroup()
  if (identical(before, grouped_preds$cv_group)) break
}
grouped_preds <- grouped_preds |>
  mutate(abs_residual = abs(y_true - y_pred), cv_group = factor(cv_group))

interval_levels <- c(0.80, 0.95)
interval_radii <- crossing(cv_group = levels(grouped_preds$cv_group), level = interval_levels) |>
  rowwise() |>
  mutate(radius = quantile(grouped_preds$abs_residual[grouped_preds$cv_group != cv_group],
                           probs = level, names = FALSE, type = 8)) |>
  ungroup()
interval_eval <- grouped_preds |>
  crossing(level = interval_levels) |>
  left_join(interval_radii, by = c("cv_group", "level")) |>
  mutate(covered = abs_residual <= radius,
         width = 2 * radius)

# Group bootstrap respects the unit held out by the main evaluation.
coverage_ci <- function(data, level, n_boot = 500) {
  d <- filter(data, .data$level == level)
  members <- split(d$covered, d$cv_group)
  set.seed(as.integer(level * 1000))
  boot <- replicate(n_boot, {
    picked <- sample(members, length(members), replace = TRUE)
    mean(unlist(picked, use.names = FALSE))
  })
  c(lo = unname(quantile(boot, 0.025)), hi = unname(quantile(boot, 0.975)))
}
group_coverage <- function(data, level) {
  d <- filter(data, .data$level == level)
  mean(vapply(split(d$covered, d$cv_group), mean, numeric(1)))
}
interval_summary <- tibble(level = interval_levels) |>
  rowwise() |>
  mutate(coverage = group_coverage(interval_eval, level),
         ci = list(coverage_ci(interval_eval, level)),
         lo = ci["lo"], hi = ci["hi"],
         half_width_fold = 10^median(interval_eval$radius[interval_eval$level == level])) |>
  ungroup() |>
  mutate(label = sprintf("%.0f%% [%.0f, %.0f]", 100 * coverage, 100 * lo, 100 * hi),
         nominal = factor(sprintf("%d%%", 100 * level), levels = sprintf("%d%%", 100 * interval_levels)))

p_h <- ggplot(interval_summary, aes(nominal, coverage)) +
  geom_hline(aes(yintercept = level), linetype = "dashed", colour = "grey60") +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.55, colour = "#332288") +
  geom_point(size = 1.8, colour = "#332288") +
  geom_text(aes(label = label), vjust = -0.65, size = 1.7) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1.12),
                     breaks = c(0, 0.5, 0.8, 0.95, 1)) +
  scale_x_discrete(labels = setNames(sprintf("%d%%\n±%.1f-fold", 100 * interval_summary$level,
                                             interval_summary$half_width_fold),
                                      interval_summary$nominal)) +
  labs(x = "Nominal interval (median half-width)", y = "Group-weighted coverage") +
  theme_panel(axis.text.x = element_text(size = 5.0),
              axis.title.x = element_text(size = 6.0),
              axis.title.y = element_text(size = 6.0))

# i: where the absolute half-life error sits. These are descriptive strata of the same
# held-out predictions, not subgroup-specific models. Source is represented by the amount
# of evidence contributed by its source identifier, rather than a long, unstable list of
# individual PMIDs; length is a practical molecular-size proxy in this peptide-only cohort.
strata_data <- hl_preds |>
  select(id, y_true, y_pred, species, cluster, PMID) |>
  left_join(hl_meta |>
              select(id, route_class, assay_method, hl_engineering, n_residues), by = "id") |>
  mutate(fold_error = 10^abs(y_true - y_pred),
         half_life_range = cut(y_true,
                                breaks = quantile(y_true, c(0, 1/3, 2/3, 1)),
                                include.lowest = TRUE,
                                labels = c("short", "middle", "long")),
         length_tertile = cut(n_residues,
                              breaks = quantile(n_residues, c(0, 1/3, 2/3, 1)),
                              include.lowest = TRUE,
                              labels = c("short", "middle", "long")))
source_sizes <- strata_data |> count(PMID, name = "source_n")
strata_data <- strata_data |>
  left_join(source_sizes, by = "PMID") |>
  mutate(source_stratum = case_when(
    source_n == max(source_n) ~ "largest source",
    source_n >= 4 ~ "other sources, n ≥ 4",
    source_n >= 2 ~ "sources, n = 2–3",
    TRUE ~ "sources, n = 1"
  ))
# Route is NOT a section here any more: every record in the shipped cohort is intravenous,
# so the stratum would have one level carrying all 213 rows and say nothing.
error_strata <- bind_rows(
  transmute(strata_data, section = "Species", stratum = species, fold_error),
  transmute(strata_data, section = "Assay", stratum = assay_method, fold_error),
  transmute(strata_data, section = "Source contribution", stratum = source_stratum, fold_error),
  transmute(strata_data, section = "Observed half-life", stratum = half_life_range, fold_error),
  transmute(strata_data, section = "Peptide length", stratum = length_tertile, fold_error),
  transmute(strata_data, section = "Engineering", stratum = hl_engineering, fold_error)
) |>
  group_by(section, stratum) |>
  summarise(n = n(), median = median(fold_error), q25 = quantile(fold_error, 0.25),
            q75 = quantile(fold_error, 0.75), .groups = "drop") |>
  # A stratum of a handful of records has no stable interquartile range and invites the eye
  # to weigh it against strata a hundred times its size. The route restriction makes several
  # strata this small, so they are dropped rather than drawn.
  filter(n >= 5) |>
  mutate(label = sprintf("%s (n = %d)", stratum, n),
         section = factor(section, levels = c("Species", "Assay", "Source contribution",
                                               "Observed half-life", "Peptide length", "Engineering")))

p_i <- ggplot(error_strata, aes(median, reorder(label, median))) +
  geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.4) +
  geom_linerange(aes(xmin = q25, xmax = q75, colour = section), linewidth = 0.6) +
  geom_point(aes(colour = section), size = 1.5) +
  # Mirrored left-right, like the feature-block panel: error grows leftwards, and the row labels
  # and section names sit on the right.
  facet_grid(section ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(section = \(x) recode(x, `Source contribution` = "Source\ncontrib.",
                                                       `Observed half-life` = "Observed\nhalf-life",
                                                       `Peptide length` = "Peptide\nlength"))) +
  scale_x_continuous(transform = scales::transform_compose(scales::transform_log10(),
                                                           scales::transform_reverse()),
                     breaks = c(1, 10), labels = c("1×", "10×")) +
  # Labels wrapped to two lines and the section names set small and rotated, so the plotting
  # area keeps most of the column.
  scale_y_discrete(position = "right", labels = \(x) str_wrap(x, width = 26)) +
  scale_colour_manual(values = c("Species" = "#4169E1",
                                 "Assay" = "#CC6677", "Source contribution" = "#AA4499",
                                 "Observed half-life" = "#DDCC77", "Peptide length" = "#44AA99",
                                 "Engineering" = "#882255"), guide = "none") +
  # The IQR meaning of the bars is stated in the caption; spelling it out here overran the
  # half-page column width.
  labs(x = "Median\nfold error", y = NULL, caption = "seed 0") +
  theme_panel(axis.text.y = element_text(size = 4.4, lineheight = 0.85),
              axis.text.x = element_text(size = 5.4),
              axis.title.x = element_text(size = 5.8),
              strip.text.y.right = element_text(angle = -90, size = 4.6, lineheight = 0.85),
              strip.placement = "outside",
              strip.background = element_blank(),
              panel.spacing.y = unit(1.2, "pt")) +
  # Minimal plot margin so the panel's actual data region keeps as much of the half-page
  # column as possible once the stratum labels and section strips are drawn.
  theme(plot.margin = margin(6, 2, 2, 12),
        plot.caption = element_text(size = 4.4, colour = "grey45", hjust = 0))

# j: which feature blocks earn their place. Every combination of the three available blocks,
# scored by grouped 5-fold CV over every record (MuriVerse/evaluate_halflife_model.py), under
# two grouping schemes at once. This panel is the argument for the model's composition, and
# the composition is read off this measurement rather than asserted: `shipped` in
# 05_data_halflife_model.R is whichever combination scores best on held-out homology
# clusters, so the marker cannot drift from what the evaluator actually chose.
# Held-out clusters score LOWER than held-out sources because homology clusters bundle
# peptides sharing a half-life regime, so removing one removes a label range - visible in the
# baseline row, which is further below zero under cluster grouping.
pal_grouping <- c("held-out clusters and sources" = "#332288",
                  "held-out sequence clusters"    = COL_MV_TEAL,
                  "held-out sources"              = "#DDCC77")
GROUP_PRIMARY <- "held-out clusters and sources"

ablation_order <- hl_ablation |>
  filter(grouping == GROUP_PRIMARY) |>
  arrange(spearman) |>
  mutate(lab = if_else(shipped, paste0(features, "  \u25c0 shipped"), as.character(features)))

# The two rows the caption's central claim rests on: the shipped model reaches rho ~0.56, but
# the molecule covariates on their own reach ~0. Both were previously just two of sixteen
# undifferentiated rows, so the figure's strongest result was invisible in it. They are marked
# here from the measured table, not hard-coded, so the annotation cannot drift from the data.
MOLECULE_ROW <- "molecule annotation only"
# Row positions on the discrete axis, taken from the ordering above so the bands track the
# sort rather than being placed by eye.
band_rows <- tibble(
  features = c(ablation_order$features[ablation_order$shipped], MOLECULE_ROW),
  kind = c("shipped", "molecule")) |>
  mutate(idx = match(features, rev(ablation_order$features)))

# Rows run from the lowest rho at the TOP (molecule annotation only, red) to the highest at the
# BOTTOM (the shipped configuration, green). A discrete axis draws its first level at the
# bottom, so the levels are the descending order; band_rows above indexes into the same one.
p_j <- hl_ablation |>
  mutate(features = factor(features, levels = rev(ablation_order$features))) |>
  ggplot(aes(spearman, features)) +
  # The two rows the caption's argument turns on, banded so they are findable among sixteen:
  # what ships, and what the molecule alone is worth.
  geom_rect(data = band_rows, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = idx - 0.5, ymax = idx + 0.5, fill = kind)) +
  scale_fill_manual(values = c(shipped = "grey92", molecule = "#F6D5D3"), guide = "none") +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_linerange(data = \(d) filter(d, grouping == GROUP_PRIMARY),
                 aes(xmin = spearman_lo, xmax = spearman_hi), colour = "#332288",
                 linewidth = 0.35, alpha = 0.55) +
  geom_line(aes(group = features), colour = "grey70", linewidth = 0.5) +
  geom_point(aes(colour = grouping), size = 1.9) +
  scale_colour_manual(values = pal_grouping,
                      labels = c("held-out clusters and sources" = "clusters + sources",
                                 "held-out sequence clusters" = "clusters",
                                 "held-out sources" = "sources")) +
  # rho runs low to high from the left and the row labels sit on the right, so the long labels
  # have the free side of the panel instead of squeezing the data.
  scale_y_discrete(labels = setNames(ablation_order$lab, ablation_order$features),
                   position = "right") +
  scale_x_continuous(limits = c(-0.4, 0.8), breaks = c(0, 0.5), expand = expansion(mult = 0.02)) +
  labs(x = expression("Cross-validated Spearman "*rho), y = NULL, colour = NULL) +
  # labs(title = "Feature blocks, ranked",
  #      subtitle = local({
  #        a <- hl_ablation |> filter(grouping == GROUP_PRIMARY)
  #        best <- a[a$shipped, ]
  #        emb <- a[a$features == "sequence + SMILES", ]
  #        sprintf("%d compact columns vs %s embedding dims",
  #                best$dims, format(emb$dims, big.mark = ","))
  #      })) +
  # The molecule-only result, called out on the band itself: this is the number the caption's
  # central claim rests on, and as a bare row among sixteen it was invisible.
  # Short form: at half-page width the longer "molecule alone ..." wording ran back over the
  # row labels, and the banded row already names what is being quantified.
  # Both labels are written into the empty side of their row: the molecule-only row (rho ~0.09)
  # has its free space at the right, the shipped row (rho ~0.58) at the left.
  geom_text(data = \(d) filter(d, grouping == GROUP_PRIMARY, features == MOLECULE_ROW),
            aes(x = 0.78, label = sprintf("\u03c1 = %.2f", spearman)),
            hjust = 1, size = 1.8, colour = "grey25", inherit.aes = TRUE) +
  # The shipped model's headline rho, in green.
  geom_text(data = \(d) filter(d, grouping == GROUP_PRIMARY, shipped),
            aes(x = -0.38, label = sprintf("\u03c1 = %.2f", spearman)),
            hjust = 0, size = 2.0, fontface = "bold", colour = COL_MV_GREEN,
            inherit.aes = TRUE) +
  guides(colour = guide_legend(nrow = 1)) +
  theme_panel(axis.text.y = element_text(size = 5.8), plot.title.position = "plot") +
  theme(plot.title = element_text(size = SZ_TITLE, hjust = 0.05),
        plot.subtitle = element_text(size = SZ_SUB, hjust = 0.05),
        # Labels are on the right now, so the right edge needs no extra margin; a little at the
        # left keeps the figure letter clear of the axis.
        plot.margin = margin(8, 4, 2, 12)) +
  # A column to the right of the row labels, running down from the top right of the panel.
  # The row labels run to the right of the plot, and the top ones are short, so the key sits in
  # that gap rather than taking a column of its own; x is in panel widths, so it lies outside
  # the panel.
  # Along the bottom, in one row. The previous position was given in panel widths (x = 3.7)
  # to sit in the gap beside the old column-3 stack; after the route restriction changed that
  # column it landed on the neighbouring panel, and the panel's own corners are not actually
  # free - sixteen rows of points reach into all of them.
  theme(legend.position = "bottom", legend.location = "plot",
        legend.key = element_blank(), legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        legend.key.height = unit(6, "pt"), legend.key.width = unit(6, "pt"),
        legend.spacing.y = unit(0, "pt"), legend.spacing.x = unit(2, "pt"),
        legend.text = element_text(size = 5.2, margin = margin(l = 1, r = 4)))

# b: label shift across the saved cluster-preserving train/val split. The two densities
# share whole homology clusters between training and validation, so the mean label of each
# split is a property of the split itself, not noise - holding out whole peptide families
# removes a label range from the training set, which is why the cluster-grouped baseline in
# the feature-block panel (g) sits further below zero than the source-grouped one. Paired with panel a to make
# the split's label consequence legible alongside its group counts.
split_means <- hl_meta |> group_by(split) |> summarise(m = mean(log_label), .groups = "drop")
p_b <- ggplot(hl_meta, aes(log_label, fill = split)) +
  geom_density(alpha = 0.55, colour = NA) +
  geom_vline(data = split_means, aes(xintercept = m, colour = split),
             linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
  # Panel a carries the training/validation legend above the stacked pair, so suppress it here
  # to avoid showing the same key twice in the same column.
  scale_fill_manual(values = c(train = "#332288", val = "#DDCC77"), guide = "none") +
  scale_colour_manual(values = c(train = "#332288", val = "#DDCC77"), guide = "none") +
  scale_x_loghl(breaks = c(1 / 60, 5 / 60, 0.25, 1, 6, 24)) +
  labs(x = "Half-life", y = "Density") +
  theme_panel()

# c/d/e: is the pooled half-life label contaminated by flip-flop kinetics? After an
# intravenous dose the terminal half-life is an elimination half-life. After a subcutaneous
# dose it is only an elimination half-life if absorption is the faster step; when absorption
# is slower the terminal phase reports the ABSORPTION rate instead. The cohort pools both
# routes into one `half_life_h` column, so if flip-flop is present the target is a mixture of
# two different physical quantities. MuriVerse/evaluate_route_flipflop.py supplies all three
# panels; groups are connected components of (cluster, PMID) over the full cohort, so no
# cluster or source leaks between a training pool and the evaluated rows.
flip_pairs <- read_csv(file.path(DIR_TABS,
                                 "table32_muriverse_halflife_route_flipflop_pairs.csv"),
                       show_col_types = FALSE)
route_seeds <- read_csv(file.path(DIR_TABS,
                                  "table33b_muriverse_halflife_route_transfer_by_seed.csv"),
                        show_col_types = FALSE)
route_size_ctl <- read_csv(file.path(DIR_TABS,
                                     "table33c_muriverse_halflife_route_size_control.csv"),
                           show_col_types = FALSE)

# c: the only comparison here that isolates the route. Same sequence, same species, both
# routes, so the molecule and the animal are fixed and the route is all that changes. The
# marginal SC/IV median gap is ~4-fold, but most of that is composition - the two routes
# were used on different peptides - and the paired shift is the smaller, real effect.
flip_n <- nrow(flip_pairs)
flip_longer <- sum(flip_pairs$sc > flip_pairs$iv)
flip_p <- flip_pairs$wilcoxon_p[1]
flip_med <- flip_pairs$median_delta[1]

p_c <- ggplot(flip_pairs, aes(iv, sc, colour = species)) +
  geom_abline(linetype = "dashed", colour = "grey60") +
  geom_point(size = 1.1, alpha = 0.85) +
  scale_colour_manual(values = c(mouse = MV_MOUSE, rat = MV_RAT), name = NULL,
                      labels = c(mouse = "mouse", rat = "rat")) +
  coord_equal(xlim = CD_LIM, ylim = CD_LIM, expand = FALSE) +
  scale_x_loghl(breaks = c(1 / 60, 5 / 60, 0.25, 1, 6, 24)) + scale_y_loghl(breaks = c(1 / 60, 5 / 60, 0.25, 1, 6, 24)) +
  labs(x = "IV half-life", y = "SC half-life") +
  # Bottom-right, inside the panel, in the same grey as the calibration panel's pooled
  # figures: every
  # paired point sits on or above the identity line, so that corner is always empty.
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("%d/%d above identity\np = %.3f", flip_longer, flip_n, flip_p),
           hjust = 1.06, vjust = -0.35, size = 1.6, colour = "grey30") +
  guides(colour = guide_legend(ncol = 1)) +
  theme_panel() +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1), legend.key.size = unit(5, "pt"),
        legend.text = element_text(size = 4.8), legend.background = element_blank(),
        legend.margin = margin(0, 0, 0, 0))

# e: cross-species AND cross-route transfer, crossed. The earlier version of these two
# panels compared IV against SC pooled over both animals, which cannot separate a species
# effect from a route effect - and in this cohort the two are confounded, rat being 71%
# intravenous against mouse's 46%. So the cohort is cut into four species x route cells and
# each cell's held-out rows are scored four ways under one set of grouped folds. The training
# pools form a 2x2 of what they share with the target:
#
#   route/species-matched (evaluator key "same cell")  - the specialist
#   cross-species  other species, SAME route - does the animal transfer?
#   cross-route    same species, OTHER route - does the route transfer?
#   pooled         everything                - what the shipped model actually trains on
#
# Reading cross-species against cross-route within one cell separates the two axes directly.
# Routes are grouped IV against EV (extravascular: sc, ip, im) because only an intravenous
# dose has no absorption phase, and because ip (22 rows) and im (12) cannot stand alone.
# MuriVerse/evaluate_species_route_cells.py supplies both panels.
cell_seeds <- read_csv(file.path(DIR_TABS,
                                 "table36b_muriverse_halflife_species_route_cells_by_seed.csv"),
                       show_col_types = FALSE)

# Cells ordered by size, largest at the top, so the eye meets the cell that carries the
# cohort first - rat IV is 179 of the 320 rows and effectively sets the headline score.
CELL_LEVELS <- c("mouse EV", "mouse IV", "rat EV", "rat IV")
CELL_N <- c("rat IV" = 179, "rat EV" = 67, "mouse EV" = 40, "mouse IV" = 34)
cell_lab <- setNames(sprintf("%s\n(n = %d)", names(CELL_N), CELL_N), names(CELL_N))
# "same cell" is the key the evaluator, the tables and the number audit all use, so it is
# relabelled here at display time only. The label says what the pool shares with the target -
# both its route and its species - which "same cell" left for the reader to work out.
TRAIN_KEYS   <- c("same cell", "cross-species", "cross-route", "pooled")
TRAIN_LEVELS <- c("route/species-matched", "cross-species/route-fixed",
                  "cross-route/species-fixed", "pooled")
pal_train <- c("route/species-matched" = "#332288",
               "cross-species/route-fixed" = "#CC6677",
               "cross-route/species-fixed" = "#DDCC77",
               "pooled" = "#44AA99")

cell_tr <- cell_seeds |>
  mutate(cell = paste(species, admin),
         cell = factor(cell_lab[cell], levels = unname(cell_lab[CELL_LEVELS])),
         train = factor(TRAIN_LEVELS[match(train, TRAIN_KEYS)], levels = TRAIN_LEVELS))
cell_mean <- cell_tr |>
  group_by(cell, train) |>
  summarise(across(c(spearman, bias), mean), .groups = "drop")
DODGE_R <- position_dodge(width = 0.72)

cell_panel <- function(metric, xlab) {
  ggplot(cell_tr, aes(.data[[metric]], cell, colour = train)) +
    geom_point(position = DODGE_R, size = 0.7, alpha = 0.45) +
    geom_point(data = cell_mean, position = DODGE_R, size = 1.8) +
    scale_colour_manual(values = pal_train, name = NULL) +
    labs(x = xlab, y = NULL) +
    guides(colour = guide_legend(nrow = 2)) +
    theme_panel(axis.text.y = element_text(size = 4.8)) +
    theme(legend.position = "bottom", legend.key.size = unit(5, "pt"),
          legend.text = element_text(size = 4.6),
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"))
}

# e: one panel, two facets side by side - ranking on the left, signed bias on the right.
# They share the four cells, the colour key and one legend, and differ only in the metric on
# x, so they are facets rather than two panels: the reader compares a cell's ranking with its
# bias by moving across one row. Free x because rho and a log10 offset share no scale; the
# metric is named in the strip below each facet. The zero line matters in both - two cells
# rank below chance on the left, and on the right zero is the unbiased reference.
# No significance facet is shown: the rows of a grouped-CV pass are not independent, so a
# p-value on the Spearman rho would overstate the evidence and is omitted.
CELL_FACETS <- c(spearman = "Held-out Spearman \u03c1", bias = "Bias, log10 h")
long_cells <- function(d) {
  d |> pivot_longer(any_of(names(CELL_FACETS)), names_to = "metric", values_to = "value") |>
    mutate(metric = factor(CELL_FACETS[metric], levels = unname(CELL_FACETS)))
}
cell_seeds_long <- long_cells(cell_tr)
# Mean and standard deviation of the five seeds, per cell, training pool and metric. The SD
# is the spread of one model re-scored under five reshuffled fold assignments, so it measures
# how much the score depends on which groups land in which fold; it is not a confidence
# interval for the cohort, which five correlated passes over the same rows cannot give.
cell_mean_long <- cell_seeds_long |>
  group_by(cell, train, metric) |>
  summarise(mean = mean(value), sd = sd(value), .groups = "drop")

# Training pools take the four colours `pal_route` was defined with, because blue and green
# already mean mouse and rat elsewhere in the figure: matched navy, cross-species rose,
# cross-route sand, pooled magenta. (That palette no longer colours a panel of its own - the
# calibration panel splits by species now that every shipped row is intravenous - but it
# remains a well-separated four-colour set and is kept for this key.)
pal_train_bg <- setNames(unname(pal_route[c("iv", "sc", "ip", "im")]), TRAIN_LEVELS)

# Dodged y position of each dot row: position_dodge(0.72) puts the first training pool at the
# bottom of its cell, spaced 0.72 / 4 apart.
guide_rows <- expand_grid(cell = levels(cell_tr$cell), train = TRAIN_LEVELS) |>
  mutate(y = match(cell, levels(cell_tr$cell)) + (match(train, TRAIN_LEVELS) - 2.5) * 0.72 / 4)

# The five seeds are summarised as mean +/- 1 SD (a point with a horizontal bar), drawn first so
# the individual seeds sit on top of it and a seed that lands on the mean stays visible.
p_e <- ggplot(cell_seeds_long, aes(value, cell, colour = train)) +
  # A faint guide through each dot row (one per cell x training pool, at its dodged position),
  # drawn first so it sits below everything else.
  geom_hline(data = guide_rows, aes(yintercept = y), colour = "grey40", alpha = 0.2,
             linewidth = 0.3, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_pointrange(data = cell_mean_long, aes(x = mean, xmin = mean - sd, xmax = mean + sd),
                  position = DODGE_R, size = 0.35, linewidth = 0.5, fatten = 4) +
  geom_point(position = DODGE_R, size = 0.7, alpha = 0.45) +
  facet_wrap(~metric, nrow = 1, scales = "free_x", strip.position = "bottom") +
  # Rows are the TEST set (the held-out rows of one species x route cell); colours are four
  # separately fitted models, named for what they were TRAINED on. Both are said on the panel.
  scale_colour_manual(values = pal_train_bg, name = "Trained on:") +
  scale_x_continuous(breaks = scales::breaks_extended(n = 4),
                     expand = expansion(mult = c(0.12, 0.12))) +
  labs(x = NULL, y = "Tested on (held-out rows)") +
  guides(colour = guide_legend(nrow = 1)) +
  theme_panel(axis.text.y = element_text(size = 4.8)) +
  theme(strip.placement = "outside", strip.background = element_blank(),
        strip.text = element_text(size = 5.6),
        legend.position = "bottom", legend.key.size = unit(5, "pt"),
        legend.text = element_text(size = 4.6),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"))

# b: the contamination the restriction removes, and the reason the shipped model is fitted on
# intravenous records only. Panel a scores four training pools against each test cell; this
# panel reads the BIAS of two of them - the two that isolate one axis each - and shows that
# the axes are not comparable.
#
#   cross-species  other animal, SAME route   - does the ANIMAL transfer?
#   cross-route    same animal, OTHER route   - does the ROUTE transfer?
#
# Every cross-route transfer is offset by 0.45-0.65 log10, a factor of 3-4, and the sign is
# the one flip-flop kinetics predicts: after an extravascular dose the terminal phase reports
# the ABSORPTION rate whenever absorption is the slower step, so EV labels are systematically
# LONGER. An IV-trained model therefore UNDER-predicts extravascular half-lives, and an
# EV-trained model OVER-predicts intravenous ones. Every cross-species transfer, by contrast,
# stays inside +/-0.21. That the route offset reproduces independently in mouse and in rat, at
# near-identical magnitude and with opposite signs either way, is what makes it structural
# rather than sampling noise - it cannot be an artefact of species composition.
#
# The practical conclusion is asymmetric, and it is the whole argument for the cohort used
# from panel c onwards: the two species very nearly share an absolute half-life scale and the
# two routes do not, so a route-pooled target is a mixture of two different physical
# quantities. Species may be pooled; route may not.
FOLD_BREAKS <- c(-0.7, -0.3, 0, 0.3, 0.7)
fold_lab <- function(x) ifelse(abs(x) < 1e-9, "0",
                               sprintf("%+.1f\n(%.1f\u00d7)", x, 10^abs(x)))

bias_axes <- cell_seeds |>
  filter(train %in% c("cross-species", "cross-route")) |>
  mutate(cell = paste(species, admin),
         cell = factor(cell_lab[cell], levels = unname(cell_lab[CELL_LEVELS])),
         axis = factor(train, levels = c("cross-species", "cross-route"),
                       labels = c("other animal, same route", "same animal, other route")))
bias_mean <- bias_axes |>
  group_by(cell, axis) |>
  summarise(mean = mean(bias), sd = sd(bias), .groups = "drop")

# The species envelope, drawn as a band so the reader sees at a glance which points fall
# inside it. Taken from the measured cross-species spread rather than asserted.
SPECIES_ENV <- max(abs(bias_mean$mean[bias_mean$axis == "other animal, same route"]))

pal_axis <- c("other animal, same route" = MV_RAT, "same animal, other route" = "#CC6677")

p_b_route <- ggplot(bias_axes, aes(bias, cell, colour = axis)) +
  annotate("rect", xmin = -SPECIES_ENV, xmax = SPECIES_ENV, ymin = -Inf, ymax = Inf,
           fill = MV_RAT, alpha = 0.10) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.4) +
  geom_pointrange(data = bias_mean, aes(x = mean, xmin = mean - sd, xmax = mean + sd),
                  position = DODGE_R, size = 0.3, linewidth = 0.5, fatten = 4) +
  geom_point(position = DODGE_R, size = 0.6, alpha = 0.45) +
  scale_colour_manual(values = pal_axis, name = "Trained on:") +
  scale_x_continuous(breaks = FOLD_BREAKS, labels = fold_lab,
                     expand = expansion(mult = 0.08)) +
  # Title and subtitle are sized for a one-third-page column: the longer wording ran off the
  # panel's right edge at 180 mm.
  labs(x = "Bias, log10 h (fold)", y = "Tested on (held-out rows)",
       title = "Route does not transfer") +
  # The summary that used to be the subtitle, moved into the panel's empty bottom-right
  # corner in light grey: it is a readout of the points, so it belongs beside them rather
  # than in the heading, and the heading reads as a claim without it.
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("cross-route 3-4\u00d7\ncross-species within \u00b1%.2f", SPECIES_ENV),
           hjust = 1.04, vjust = -0.3, size = 1.45, colour = "grey45", lineheight = 0.95) +
  # The two arrows name the DIRECTION of the contamination, which is the half of the result a
  # bare magnitude cannot carry: the sign is what flip-flop predicts, and it flips with the
  # direction of transfer.
  annotate("text", x = -0.78, y = 4.52, label = "IV-trained\nunder-predicts EV",
           size = 1.4, colour = "grey30", hjust = 0, lineheight = 0.9) +
  annotate("text", x = 0.78, y = 4.52, label = "EV-trained\nover-predicts IV",
           size = 1.4, colour = "grey30", hjust = 1, lineheight = 0.9) +
  coord_cartesian(xlim = c(-0.8, 0.8), ylim = c(0.6, 4.75), clip = "off") +
  guides(colour = guide_legend(nrow = 2)) +
  theme_panel(axis.text.y = element_text(size = 4.8),
              axis.text.x = element_text(size = 4.6, lineheight = 0.85)) +
  theme(plot.title = element_text(size = 6.2),
        legend.position = "bottom", legend.key.size = unit(5, "pt"),
        legend.text = element_text(size = 4.6), legend.title = element_text(size = 4.8),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"))

# Layout: a three-column grid, lettered in reading order, using the historical object names in
# brackets. Rows 1 and 2 establish the cohort the model is fitted on; row 3 is that model.
#   Row 1: a species x route transfer (p_e) across columns 1-2 with two facets (ranking, bias),
#          b the route contamination it implies (p_b_route) in column 3
#   Row 2: c held-out prediction scatter (p_f; mouse beside rat), d calibration by species
#          (p_d), e ROC pooled over species (p_g)
#   Row 3: f error stratification (p_i), g feature blocks (p_j, widened for its long row labels),
#          h paired flip-flop evidence (p_c) over i label shift (p_b), stacked in column 3
# Row 2 is now the shipped model judged three ways on its own held-out rows - ordering,
# absolute scale, and ranking at the 6 h target - with the paired flip-flop evidence moved
# down to h, beside the other cohort-level diagnostics.
# Two panels are built above and deliberately left out of the layout: p_h (group-weighted
# interval coverage) and p_a (the saved 80/20 split in counts). p_a characterises a split that
# is not used for any performance number in the figure, and the route restriction cost the
# layout a cell; add either back here when it has a place.
LAB <- list(size = 9, face = "bold")
fig3 <- ggarrange(
  ggarrange(p_e, p_b_route, ncol = 2, widths = c(2, 1),
            labels = c("a", "b"), font.label = LAB),
  ggarrange(p_f, p_d, p_g, ncol = 3, widths = c(1, 1, 1),
            labels = c("c", "d", "e"), font.label = LAB),
  ggarrange(p_i, p_j, ggarrange(p_c, p_b, nrow = 2, labels = c("h", "i"), font.label = LAB),
            ncol = 3, widths = c(0.8, 1.4, 0.82),
            labels = c("f", "g", ""), font.label = LAB),
  nrow = 3, heights = c(62, 62, 93)
)

# 180 mm is save_fig's default and the single-page width the rest of the pipeline targets.
save_fig(fig3, "fig3", width = 180, height = 217)
