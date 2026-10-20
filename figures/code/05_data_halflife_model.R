# The single MuriVerse half-life regressor on the strict rodent cohort.
#
# Two design choices matter for what the figure can honestly show.
#
# First, there is ONE model, not separate sequence, SMILES and one-hot branches: fig3.R is built
# around a single concatenated model and a feature-block ablation panel.
#
# Second, the predictions are OUT-OF-FOLD, not from one held-out draw. The strict cohort is
# 402 rows over 101 homology clusters; a single 80/20 split leaves ~81 validation records,
# and evaluate_halflife_model.py records what that costs - a mouse Spearman of 0.00 on one
# split against ~0.5 under cross-validation over the same clusters. Every point plotted here
# was predicted by a model that saw neither the record nor any member of its homology
# cluster, which is the only version of "held out" this cohort is large enough to support.
if (!exists("ROOT")) {
  source(file.path(system("git rev-parse --show-toplevel", intern = TRUE), "figures/code/00_setup.R"))
}

DIR_MV <- file.path(ROOT, "MuriVerse", "rodent", "strict")

# ---------------------------------------------------------------------------------------
# THE COHORT THE SHIPPED MODEL IS FITTED ON: intravenous records only, both species.
#
# Figure 3's panel a measures why. The cohort splits into four species x route cells, and the
# two axes behave in opposite ways. Pooling ACROSS SPECIES is nearly free - every cross-species
# transfer holds its bias inside +/-0.21 log10 - while pooling ACROSS ROUTE is not: every
# cross-route transfer is offset by 0.45-0.65 log10 (a factor of 3-4), with the sign flip-flop
# kinetics predicts. After an intravenous dose the terminal half-life is an elimination
# half-life; after an extravascular dose it reports the ABSORPTION rate whenever absorption is
# the slower step. Pooling the routes therefore gives the regressor a target that is a mixture
# of two different physical quantities, and panel b shows the size of that contamination.
#
# Route is also a free lever in a route-pooled cohort: scoring every record by its own route's
# median half-life, with no molecule involved at all, already reaches Spearman rho = 0.42
# against the pooled label. So part of the old route-pooled headline was the model recovering
# which route a record came from. Restricting to IV removes both the mixed target and the
# lever. LIKE FOR LIKE - the same 213 IV rows, the same SHIPPED configuration, the same
# grouped protocol - it is better on every axis: rho 0.507 -> 0.574, AUC(>= 6 h) 0.785 ->
# 0.818, RMSE 0.699 -> 0.675, and the mouse stratum (54% extravascular) 0.409 -> 0.562.
#
# That rho gain is NOT on its own statistically resolvable: a paired bootstrap over the 55
# connected cluster/source groups puts it at +0.063 [-0.005, +0.157]. An earlier version of
# this comment read "0.514 -> 0.593 (+0.079 [+0.002, +0.194]), AUC 0.781 -> 0.826, RMSE 0.688
# -> 0.664, mouse 0.378 -> 0.634", which was wrong twice over: the right-hand numbers belong
# to `physchem + study context + termini`, the UNCONSTRAINED ablation winner that is not
# shipped (it cannot be applied - see REQUIRED_SHIPPED_BLOCKS), and 0.688/0.378 are the
# pooled model over all 320 rows rather than over the 213 being compared. Do not restore it.
#
# The case for the restriction is therefore the physics, the cross-route bias measured in
# panels a and b, and the removal of the route lever - not the headline correlation.
#
# HL_ROUTE is the switch. "iv" reads the artifacts written by
#   python MuriVerse/evaluate_halflife_model.py --route iv --nested --jobs 8
# and "" reads the route-pooled ones, which are left on disk unchanged for comparison.
if (!exists("HL_ROUTE")) HL_ROUTE <- "iv"
HL_TAG <- if (nzchar(HL_ROUTE)) paste0("_", HL_ROUTE) else ""
HL_ROUTES_KEPT <- if (nzchar(HL_ROUTE)) strsplit(HL_ROUTE, ",")[[1]] else NULL

# Keep the predeclared 6 h outcome target for headline reporting. Panels 3a-b also show a
# 4 h sensitivity cutoff so its classification results cannot overwrite that target.
HL_POS_THRESHOLD_H <- 6
# Figure 3a-b sensitivity view: compare a 4 h screen without changing the predeclared
# 6 h target reported in the overall outcome and stringency analysis.
HL_SCREEN_THRESHOLD_H <- 4

# `hl_meta_all` keeps every row of the strict cohort, because Figure 2 describes the whole
# input cohort and Figure 3's panels a and b are the measurement that JUSTIFIES the
# restriction - they need the extravascular rows they argue about. `hl_meta` is the cohort
# the shipped model is actually fitted and scored on, and every performance panel uses it.
hl_meta_all <- read_csv(file.path(DIR_MV, "halflife_meta_with_split.csv"),
                        show_col_types = FALSE, guess_max = 1e5) |>
  mutate(split = factor(split, levels = c("train", "val")))
hl_meta <- if (is.null(HL_ROUTES_KEPT)) hl_meta_all else
  filter(hl_meta_all, route_class %in% HL_ROUTES_KEPT)

# One row per record, carrying the prediction from the fold that held it out. Written by
# evaluate_halflife_model.py for whichever feature combination scores best on held-out
# homology clusters; `features` names that combination.
hl_oof <- read_csv(file.path(DIR_MV, paste0("oof_predictions", HL_TAG, ".csv")),
                   show_col_types = FALSE)
# The out-of-fold rows and the cohort rows must be the same set, or every panel below is
# silently describing two different cohorts.
stopifnot(setequal(hl_oof$id, hl_meta$id))

SHIPPED_FEATURES <- unique(hl_oof$features)
stopifnot(length(SHIPPED_FEATURES) == 1)

# `hl_preds` keeps the column names the figure uses. The cluster-grouped prediction is the
# headline one: held-out source IDs are the easier test on this cohort, because a source's
# peptides can have relatives in other sources while a homology cluster's do not.
# Each seed is a complete out-of-fold pass with its own fold shuffle. Every reported number
# is the MEAN OF THE PER-SEED SCORES, never the score of the seed-averaged prediction (which
# is a small bagging ensemble and reads higher than any single deployed model would). The
# curves and scatter are drawn from seed 0, one real out-of-fold pass.
SEEDS <- 0:2
# The per-seed prediction columns under the headline grouping. Figures that DRAW seed 0 but
# QUOTE the three-seed mean need these, so they live here beside SEEDS rather than being
# rebuilt in each panel.
SEED_COLS <- sprintf("y_pred_cluster_PMID_seed%d", SEEDS)
hl_preds <- hl_oof |>
  transmute(id, species = str_to_title(species), cluster, PMID,
            y_true, y_pred = y_pred_cluster_PMID_seed0,
            y_pred_cluster = y_pred_cluster_seed0, y_pred_pmid = y_pred_PMID_seed0,
            half_life_h, positive = half_life_h >= HL_SCREEN_THRESHOLD_H)

seed_preds <- function(grouping) {
  bind_rows(lapply(SEEDS, function(k) {
    hl_oof |>
      transmute(id, species = str_to_title(species), PMID, y_true,
                y_pred = .data[[sprintf("y_pred_%s_seed%d", grouping, k)]], seed = k,
                positive = half_life_h >= HL_SCREEN_THRESHOLD_H)
  }))
}

# Rank-based AUC (Mann-Whitney), identical in definition to sklearn's roc_auc_score.
rank_auc <- function(score, positive) {
  r <- rank(score); n1 <- sum(positive); n0 <- sum(!positive)
  (sum(r[positive]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Both species, then each species: the same model's predictions cut three ways, never three
# models. Labelled "Both species" rather than "Pooled": "pooled" elsewhere means route-pooled.
strata <- function(data) {
  bind_rows(mutate(data, stratum = "Both species"), mutate(data, stratum = species)) |>
    mutate(stratum = factor(stratum, levels = c("Both species", "Mouse", "Rat")))
}
hl_strata <- strata(hl_preds)

hl_seed_strata <- strata(seed_preds("cluster_PMID"))

hl_metrics <- hl_seed_strata |>
  group_by(stratum, seed) |>
  summarise(
    n_val    = n(),
    spearman = cor(y_true, y_pred, method = "spearman"),
    pearson  = cor(y_true, y_pred),
    auroc    = rank_auc(y_pred, positive),
    rmse     = sqrt(mean((y_true - y_pred)^2)),
    mae      = mean(abs(y_true - y_pred)),
    # negative wherever the model does not beat predicting the stratum's own mean
    r2       = 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2),
    .groups  = "drop"
  ) |>
  group_by(stratum) |>
  summarise(across(c(n_val, spearman, pearson, auroc, rmse, mae, r2), mean), .groups = "drop")

# fig3's panel c prints a held-out-split rho and a CV rho. Both are now the same quantity -
# every prediction is out of fold - so `hl_cv` is the pooled row, and the panel says so.
hl_cv <- hl_metrics |> filter(stratum == "Both species")

# --- the same predictions read as a ranking problem ------------------------------------
# A regressor, scored the way a screen would use it: can it put the long-lived peptides at
# the top? The label is dichotomised at HL_POS_THRESHOLD_H and the continuous prediction is
# the score. Computed directly rather than with pROC so the curves and the summary numbers
# come from one definition.
roc_curve <- function(score, positive) {
  ord <- order(score, decreasing = TRUE)
  pos <- cumsum(positive[ord]); neg <- cumsum(!positive[ord])
  tibble(tpr = c(0, pos / max(pos, 1)), fpr = c(0, neg / max(neg, 1)))
}

pr_curve <- function(score, positive) {
  # Apply each distinct score threshold as a unit. Ordering individual records within a
  # tied score invents a ranking the model did not produce (and can make PR jagged).
  n_positive <- max(sum(positive), 1)
  tibble(score = score, positive = positive) |>
    group_by(score) |>
    summarise(n = n(), tp = sum(positive), .groups = "drop") |>
    arrange(desc(score)) |>
    mutate(tp = cumsum(tp), seen = cumsum(n),
           recall = tp / n_positive, precision = tp / seen) |>
    transmute(recall, precision) |>
    add_row(recall = 0, precision = 1, .before = 1)
}

# Area by the trapezoid rule over the curve above; average precision is the step-wise sum
# that the precision-recall convention uses, not the trapezoid of the PR curve.
auc_of <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

average_precision <- function(score, positive) {
  # Match the threshold-wise curve: all examples tied at a score enter together.
  curve <- pr_curve(score, positive)
  sum(diff(curve$recall) * tail(curve$precision, -1))
}

hl_roc <- hl_strata |>
  group_by(stratum) |>
  reframe(roc_curve(y_pred, positive))

hl_pr <- hl_strata |>
  group_by(stratum) |>
  reframe(pr_curve(y_pred, positive))

hl_clf_metrics <- hl_seed_strata |>
  group_by(stratum, seed) |>
  summarise(
    n_val      = n(),
    n_positive = sum(positive),
    prevalence = mean(positive),
    auroc      = rank_auc(y_pred, positive),
    auprc      = average_precision(y_pred, positive),
    .groups    = "drop"
  ) |>
  group_by(stratum) |>
  summarise(across(c(n_val, n_positive, prevalence, auroc, auprc), mean), .groups = "drop") |>
  mutate(cutoff_h = HL_SCREEN_THRESHOLD_H)

# The displayed PR curves use seed 0, so their legend AP must describe those same scores.
hl_curve_ap <- hl_strata |>
  group_by(stratum) |>
  summarise(curve_ap = average_precision(y_pred, positive), .groups = "drop")

# --- feature-block ablation ------------------------------------------------------------
# The evidence for the model's composition, from evaluate_halflife_model.py. `shipped` is
# read off the data rather than hard-coded, so the flag cannot drift from what the evaluator
# actually chose.
GROUPING_LABELS <- c(cluster_PMID = "held-out clusters and sources",
                     cluster      = "held-out sequence clusters",
                     PMID         = "held-out sources")

hl_ablation <- read_csv(file.path(DIR_TABS, paste0("table31_muriverse_halflife_feature_ablation_strict", HL_TAG, ".csv")),
                        show_col_types = FALSE) |>
  filter(features != "baseline (train mean)") |>
  mutate(grouping = factor(GROUPING_LABELS[grouping], levels = unname(GROUPING_LABELS)),
         shipped  = features == SHIPPED_FEATURES)

stopifnot(any(hl_ablation$shipped))

# --- how the score depends on how strictly the test is held out ------------------------
# The compact-set row of the ablation under each grouping, the nested-selection estimate, and
# the WITHIN-source score: the same predictions with each source's mean removed from both
# the prediction and the label, so it measures only whether the model orders peptides inside
# one source identifier. The gap between the strict rows and the within-source row is the
# source offset the model cannot know for a new source.
within_study <- bind_rows(lapply(SEEDS, function(k) {
  hl_oof |>
    transmute(PMID, y_true, y_pred = .data[[sprintf("y_pred_cluster_PMID_seed%d", k)]],
              positive = half_life_h >= HL_POS_THRESHOLD_H, seed = k) |>
    group_by(seed, PMID) |> filter(n() > 1) |>
    mutate(y_true = y_true - mean(y_true), y_pred = y_pred - mean(y_pred)) |>
    ungroup()
})) |>
  group_by(seed) |>
  summarise(pearson = cor(y_true, y_pred), spearman = cor(y_true, y_pred, method = "spearman"),
            n = n(), .groups = "drop") |>
  summarise(pearson = mean(pearson), spearman = mean(spearman), auroc = NA_real_, n = first(n))

abl_all <- read_csv(file.path(DIR_TABS, paste0("table31_muriverse_halflife_feature_ablation_strict", HL_TAG, ".csv")),
                    show_col_types = FALSE)
nested <- read_csv(file.path(DIR_TABS, paste0("table31b_muriverse_halflife_nested_selection", HL_TAG, ".csv")),
                   show_col_types = FALSE)
pick <- function(g) abl_all |> filter(features == SHIPPED_FEATURES, grouping == g)
cols <- function(d, evaluation) {
  d |> transmute(evaluation = evaluation, n, pearson, pearson_lo, pearson_hi,
                 spearman, spearman_lo, spearman_hi, auc, auc_lo, auc_hi)
}
hl_stringency <- bind_rows(
  cols(pick("cluster_PMID"), "strict: clusters + sources held out"),
  cols(nested, "strict, configuration chosen inside folds"),
  cols(pick("cluster"), "clusters held out only"),
  cols(pick("PMID"), "sources held out only"),
  within_study |> transmute(evaluation = "within source (source mean removed)", n,
                            pearson, spearman, auc = NA_real_)
)
write_csv(hl_stringency, file.path(DIR_TABS, "table31c_muriverse_halflife_stringency.csv"))

# --- tables ----------------------------------------------------------------------------
write_csv(hl_metrics, file.path(DIR_TABS, "table10_muriverse_halflife_branch_metrics.csv"))
write_csv(hl_clf_metrics, file.path(DIR_TABS, "table10b_muriverse_halflife_ranking_metrics.csv"))

# The expanded peptide-scale cohort deliberately includes incomplete structures. Report
# model behaviour within each quality tier and species so the pooled score cannot hide a
# weak subgroup. These are descriptive slices of the same cluster-held-out predictions,
# not separate models or additional validation folds.
stopifnot(setequal(hl_preds$id, hl_meta$id),
          nrow(hl_preds) == nrow(hl_meta))
hl_tier_diagnostics <- hl_preds |>
  left_join(hl_meta |> select(id, structure_tier), by = "id", relationship = "one-to-one") |>
  group_by(species, structure_tier) |>
  summarise(
    records = n(),
    spearman = if (n() >= 3) cor(y_true, y_pred, method = "spearman") else NA_real_,
    rmse_log10_h = sqrt(mean((y_true - y_pred)^2)),
    .groups = "drop"
  ) |>
  arrange(species, desc(records))
write_csv(hl_tier_diagnostics,
          file.path(DIR_TABS, "table45_muriverse_performance_by_structure_tier.csv"))

write_csv(
  hl_meta |>
    count(smiles_status, name = "records") |>
    mutate(share = records / sum(records)) |>
    arrange(desc(records)),
  file.path(DIR_TABS, "table11_muriverse_halflife_smiles_construction.csv")
)

write_csv(
  hl_meta |>
    group_by(split) |>
    summarise(records = n(), sequences = n_distinct(sequence),
              hmmer_clusters = n_distinct(cluster), sources = n_distinct(PMID),
              median_half_life_h = median(half_life_h),
              mean_log10_h = mean(log_label), sd_log10_h = sd(log_label),
              .groups = "drop"),
  file.path(DIR_TABS, "table12_muriverse_halflife_split.csv")
)

COVARIATE_LABELS <- c(n_ter_mod = "N-terminal modification",
                      c_ter_mod = "C-terminal modification",
                      topology  = "topology",
                      chirality = "chirality")

hl_covariates <- hl_meta |>
  select(split, log_label, all_of(names(COVARIATE_LABELS))) |>
  mutate(chirality = recode(chirality, l = "all-L", d = "all-D", mix = "mixed L/D")) |>
  pivot_longer(all_of(names(COVARIATE_LABELS)),
               names_to = "covariate", values_to = "level") |>
  mutate(covariate = factor(COVARIATE_LABELS[covariate], levels = unname(COVARIATE_LABELS)))

write_csv(
  hl_covariates |>
    count(covariate, level, split) |>
    pivot_wider(names_from = split, values_from = n, values_fill = 0) |>
    arrange(covariate, desc(train + val)),
  file.path(DIR_TABS, "table13_muriverse_halflife_covariates.csv")
)
