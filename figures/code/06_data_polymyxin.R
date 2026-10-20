# Held-out application of the trained MuriVerse half-life model to colistin (polymyxin E)
# and polymyxin B, in mouse and in rat. Nothing is retrained here; the tables are written by
# MuriVerse/predict_polymyxins.py, which loads the models saved by
# evaluate_halflife_model.py after fitting them on every row of the strict set.
#
# The shipped configuration is whichever block combination scores best on held-out clusters
# and studies. All four polymyxin congeners are a lipidated, cyclic, mixed-chirality peptide
# dosed intravenously (or subcutaneously) into one species, so under that configuration they
# share few feature vectors and therefore few predictions. The STRUCTURAL configuration - the
# best-scoring one that includes PeptideCLM over the SMILES - can tell them apart, and the gap
# between the two is what the panels compare. The observed half-lives come from the studies
# flagged in_figure4 in MuriVerse/data/polymyxin_studies.csv.
if (!exists("ROOT")) {
  source(file.path(system("git rev-parse --show-toplevel", intern = TRUE), "figures/code/00_setup.R"))
}

pm_read <- function(file) read_csv(file.path(DIR_TABS, file), show_col_types = FALSE)

# WHICH COMPARISON IS DRAWN. The shipped model is fitted on intravenous records and predicts
# intravenous inputs, so the primary comparison uses intravenous observations only
# (PM_COMPARATOR_ROUTE = "IV", the default). Set PM_COMPARATOR_ROUTE <- "all" before sourcing
# this file for the secondary comparison that pools every route; fig4.R draws that one as a
# supplementary figure. The `_iv` tables are written by predict_polymyxins.py
# (compare_with_observed(..., iv_only = TRUE)). The summary tables 28 and 29 follow the
# primary comparison; the all-route versions carry a `_pooled` tag.
if (!exists("PM_COMPARATOR_ROUTE")) PM_COMPARATOR_ROUTE <- "IV"
stopifnot(PM_COMPARATOR_ROUTE %in% c("IV", "all"))
PM_IV <- PM_COMPARATOR_ROUTE == "IV"
PM_TAB_SUFFIX <- if (PM_IV) "_iv" else ""
PM_OUT_TAG <- if (PM_IV) "" else "_pooled"

# One configuration is drawn: the shipped model, which reads the molecule through its
# physicochemical descriptors and so returns a different number for each congener. The
# PeptideCLM-embedding variant is applied too but is NOT drawn - on molecules this far outside
# the training chemistry it extrapolates to ~1e-8 h, which would flatten every axis. Its raw
# predictions are written to table51_polymyxin_embedding_extrapolation.csv instead.
# The label names the feature configuration, so it is READ from the predictions table rather
# than written here. It used to be hard-coded, and when the shipped configuration changed
# with the move to the intravenous cohort the label silently kept describing the old one.
pm_pred_raw <- read_csv(file.path(DIR_TABS, "table20_polymyxin_predictions.csv"),
                        show_col_types = FALSE)
PM_SHIPPED_FEATURES <- pm_pred_raw |> filter(role == "shipped") |> pull(features) |> unique()
stopifnot(length(PM_SHIPPED_FEATURES) == 1)
PM_ROLE <- c(shipped = sprintf("shipped (%s)", PM_SHIPPED_FEATURES))
pal_role_pm <- setNames(COL_MV_RED, PM_ROLE[["shipped"]])

as_role <- function(x) factor(PM_ROLE[x], levels = unname(PM_ROLE))

# The strict cohort the two configurations were fitted on.
# The cohort the applied model was FITTED on, which since the move to the intravenous
# shipped model is the IV rows, not all 320. Panels that describe "the labels the model was
# trained on" read this, so the restriction has to be applied here or they would describe a
# training set the model never saw. PM_MODEL_ROUTES mirrors predict_polymyxins.MODEL_TAG.
PM_MODEL_ROUTES <- "iv"
pm_meta_all <- read_csv(file.path(ROOT, "MuriVerse", "rodent", "strict",
                                  "halflife_meta_with_split.csv"),
                        show_col_types = FALSE, guess_max = 1e5)
pm_meta <- filter(pm_meta_all, route_class %in% PM_MODEL_ROUTES)

# Short congener labels; the drug name is carried separately by facets
CONGENER_SHORT <- c("Polymyxin B1"              = "B1",
                    "Polymyxin B2"              = "B2",
                    "Colistin A (polymyxin E1)" = "E1 (colistin A)",
                    "Colistin B (polymyxin E2)" = "E2 (colistin B)")

PM_SPECIES <- c("Mouse", "Rat")

pm_pred_all <- pm_pred_raw
pm_pred <- pm_pred_all |>
  filter(role %in% names(PM_ROLE)) |>
  mutate(role = as_role(role), species = factor(str_to_title(species), levels = PM_SPECIES),
         congener = factor(CONGENER_SHORT[name], levels = unname(CONGENER_SHORT)),
         drug = factor(drug, levels = c("Colistin (polymyxin E)", "Polymyxin B")))

# Every retrieved value, with the study it came from. Only plasma or serum elimination
# half-lives of colistin or polymyxin B itself are drawn (`use_in_band`); the rest are kept in
# the table with the reason (tissue or ELF matrix, formed from the CMS prodrug, renal-failure
# substudy). The comparison value for a drug and species is the geometric mean of the dose-group
# values weighted by the animals dosed in each group.
# Who produced each observed value. The Tam lab (University of Houston: Abdelraouf 2012, He
# 2013, Manchandani 2016 and 2017) reported short polymyxin B half-lives from 1-compartment fits
# to sparse samples out to 6-7.5 h; Jiao 2023 (population PK, infected mice) reported about 4 h.
# The Monash group (R. L. Nation, J. Li, T. Velkov and colleagues; Monash University, and the
# University of South Australia for Li 2003) is the other lab with several studies. The grouping
# and each study's method are in data/polymyxin_study_groups.csv.
SOURCE_LEVELS <- c("Jiao 2023", "Tam lab", "Monash group", "other")

pm_obs_all <- pm_read("table21_polymyxin_observed_half_lives.csv") |>
  mutate(species = factor(str_to_title(species), levels = PM_SPECIES),
         drug = factor(drug, levels = levels(pm_pred$drug)),
         value_type = factor(case_when(kind == "mrt" ~ "mean residence time",
                                       str_detect(detail, "^derived") ~ "derived",
                                       TRUE ~ "reported"),
                             levels = c("reported", "derived", "mean residence time")),
         source_group = factor(source_group, levels = SOURCE_LEVELS),
         route_class = case_when(str_detect(route, "^IV") ~ "IV", str_detect(route, "^SC") ~ "SC",
                                 str_detect(route, "^IP") ~ "IP",
                                 TRUE ~ "not reported"))
# every usable value on any route: a study is "without a usable value" only if it has none here
pm_obs_any <- pm_obs_all |> filter(use_in_band, !is.na(drug))
pm_obs <- if (PM_IV) filter(pm_obs_any, route_class == "IV") else pm_obs_any
# A mean residence time is not an elimination half-life, so it is drawn but never averaged.
pm_mrt <- pm_obs_all |> filter(kind == "mrt")
if (PM_IV) pm_mrt <- filter(pm_mrt, route_class == "IV")

pm_obs_sum <- pm_read(sprintf("table49_polymyxin_observed_by_drug_species%s.csv", PM_TAB_SUFFIX)) |>
  mutate(species = factor(str_to_title(species), levels = PM_SPECIES),
         drug = factor(drug, levels = levels(pm_pred$drug)))

# One observed band per drug and species: the animal-weighted mean of the dose-group half-lives
# (line) and that mean +/- one animal-weighted SD (band; the study's own SD when a drug and
# species rest on a single dose group). Both are computed in predict_polymyxins.py.
pm_obs_band <- pm_obs_sum |>
  transmute(drug, species, gm = observed_log10_h, lo = band_low_log10_h, hi = band_high_log10_h,
            n_studies, n_groups, n_animals, band_basis)

pm_ctx <- pm_read("table26_polymyxin_training_label_context.csv") |>
  filter(role %in% names(PM_ROLE)) |>
  mutate(role = as_role(role), species = factor(str_to_title(species), levels = PM_SPECIES),
         train_median_log10_h = log10(train_median_half_life_h))

pm_err <- pm_read(sprintf("table27_polymyxin_fold_error_vs_baseline%s.csv", PM_TAB_SUFFIX)) |>
  filter(role %in% names(PM_ROLE)) |>
  mutate(role = as_role(role), species = factor(str_to_title(species), levels = PM_SPECIES),
         drug = factor(drug, levels = levels(pm_pred$drug)))

pm_degeneracy <- pm_read("table23_polymyxin_congener_degeneracy.csv") |>
  filter(role %in% names(PM_ROLE)) |>
  mutate(role = as_role(role), species = factor(str_to_title(species), levels = PM_SPECIES))

# The studies the observed half-lives come from (data/polymyxin_studies.csv, flagged in_figure4).
pm_studies <- pm_read("table50_polymyxin_studies_used.csv") |>
  mutate(
    # an organism tag in square brackets tells apart two papers of one author and year
    tag = str_extract(citation, "(?<=\\[)[^\\]]+(?=\\])"),
    base_label = str_squish(paste(str_extract(citation, "^.+?(?= et al\\.| and )"),
                                  str_extract(citation, "(19|20)\\d{2}"))),
    year = as.integer(str_extract(citation, "(19|20)\\d{2}")),
    source_group = factor(source_group, levels = SOURCE_LEVELS)) |>
  # the organism tag is shown only where it is needed to tell two studies apart
  group_by(base_label) |>
  mutate(label = if_else(n() > 1 & !is.na(tag), paste0(base_label, " (", tag, ")"), base_label)) |>
  ungroup()

# The observed value per drug, species and source, and every prediction against each of them.
pm_by_source <- pm_read(sprintf("table52_polymyxin_observed_by_source%s.csv", PM_TAB_SUFFIX)) |>
  mutate(species = factor(str_to_title(species), levels = PM_SPECIES),
         drug = factor(drug, levels = levels(pm_pred$drug)))
pm_fold_src <- pm_read(sprintf("table53_polymyxin_fold_error_by_source%s.csv", PM_TAB_SUFFIX)) |>
  filter(role %in% names(PM_ROLE)) |>
  mutate(species = factor(str_to_title(species), levels = PM_SPECIES),
         drug = factor(drug, levels = levels(pm_pred$drug)))

# The shipped model's prediction when a dose is supplied (the main predictions leave it missing).
pm_dose_fold <- pm_read("table54_polymyxin_prediction_vs_dose.csv") |>
  filter(!is.na(dose_mg_kg)) |>
  group_by(drug, dose_mg_kg) |>
  summarise(fold = median(fold_vs_dose_missing), .groups = "drop") |>
  mutate(drug = factor(drug, levels = levels(pm_pred$drug)))

pm_leak <- pm_read("table24_polymyxin_leakage_check.csv")
pm_homology <- pm_read("table25_polymyxin_surrogate_homology.csv")
pm_struct <- pm_read("table22_polymyxin_structures.csv")

# Constants the model has to beat: the training median of the same species (what a per-species
# lookup would give), of both species together (what a species-blind lookup would give) and of
# the few lipidated peptides in the training split (the class the polymyxins belong to).
pm_train <- pm_meta |> filter(split == "train")
PM_LIPID_N <- sum(pm_train$hl_engineering == "lipidated")
pm_baselines <- bind_rows(
  pm_train |> group_by(species = str_to_title(species)) |>
    summarise(h = 10^median(log_label), .groups = "drop") |>
    mutate(baseline = "median, same species"),
  tibble(species = PM_SPECIES, h = 10^median(pm_train$log_label),
         baseline = "median, both species"),
  tibble(species = PM_SPECIES,
         h = 10^median(pm_train$log_label[pm_train$hl_engineering == "lipidated"]),
         baseline = sprintf("median, lipidated (n = %d)", PM_LIPID_N))) |>
  mutate(species = factor(species, levels = PM_SPECIES))

# How far the pooled observed value moves under other reasonable choices (leave-one-study-out,
# each study once, reported values only, intravenous only).
pm_pooled_sens <- pm_read("table55_polymyxin_pooled_sensitivity.csv")

# The label distribution the applied model saw: the training split of the cohort it was
# fitted on, which under the intravenous shipped model is pm_meta (IV rows) and not all 320.
# predict_polymyxins.py takes its species-median baseline from the same restricted split, so
# the grey constants in panel b and the baseline_fold_error of table27 are one quantity.
pm_train_labels <- pm_meta |> filter(split == "train") |> transmute(log_label)

# --- the variant that is applied but not drawn -------------------------------------------
pm_embedding <- pm_pred_all |>
  filter(role == "structural_embedding") |>
  group_by(species, drug, name) |>
  summarise(median_pred_h = median(pred_half_life_h), min_pred_h = min(pred_half_life_h),
            max_pred_h = max(pred_half_life_h), .groups = "drop")
write_csv(pm_embedding, file.path(DIR_TABS, "table51_polymyxin_embedding_extrapolation.csv"))

# --- summary tables for the text -------------------------------------------------------

# Per role and drug: prediction spread across input variants, versus the observed span
pm_summary <- pm_pred |>
  group_by(role, species, drug) |>
  summarise(n_variants = n(),
            pred_min_h = min(pred_half_life_h), pred_max_h = max(pred_half_life_h),
            pred_median_h = median(pred_half_life_h),
            variant_fold_spread = max(pred_half_life_h) / min(pred_half_life_h),
            .groups = "drop") |>
  left_join(pm_obs_band |> transmute(species, drug, observed_band_low_h = 10^lo, observed_band_high_h = 10^hi,
                                     observed_geomean_h = 10^gm, observed_studies = n_studies,
                                     observed_animals = n_animals),
            by = c("species", "drug")) |>
  left_join(pm_err |> group_by(role, species, drug) |>
              summarise(model_fold_error = median(model_fold_error),
                        baseline_fold_error = first(baseline_fold_error),
                        baseline_half_life_h = first(baseline_half_life_h),
                        share_beating_baseline = mean(model_beats_baseline),
                        .groups = "drop"),
            by = c("role", "species", "drug"))

write_csv(pm_summary, file.path(DIR_TABS, sprintf("table28_polymyxin_role_drug_summary%s.csv", PM_OUT_TAG)))

# How far an arbitrary transcription choice moves the prediction, holding the molecule fixed,
# against how far the two DRUGS move it. The comparison only means anything because both are
# fold changes on the same prediction.
pm_encoding <- pm_pred |>
  group_by(role, species, name) |>
  summarise(fold = max(pred_half_life_h) / min(pred_half_life_h),
            n_variants = n(), .groups = "drop") |>
  group_by(role, species) |>
  summarise(n_variants = max(n_variants),
            max_encoding_fold_spread = max(fold),
            median_encoding_fold_spread = median(fold), .groups = "drop")

# Does either configuration separate the two drugs? Compared like with like: the same
# input-encoding variant, colistin against polymyxin B. Signed, because getting the ordering
# backwards is a different failure from not resolving them at all.
pm_discrimination <- pm_pred |>
  group_by(role, species, variant) |>
  filter(n_distinct(drug) == 2) |>
  summarise(colistin_h = median(pred_half_life_h[drug == "Colistin (polymyxin E)"]),
            polymyxin_b_h = median(pred_half_life_h[drug == "Polymyxin B"]),
            .groups = "drop") |>
  mutate(predicted_fold_gap = polymyxin_b_h / colistin_h) |>
  group_by(role, species) |>
  summarise(n_variants = n(),
            median_predicted_fold_gap = median(predicted_fold_gap),
            max_abs_predicted_fold_gap = max(pmax(predicted_fold_gap, 1 / predicted_fold_gap)),
            .groups = "drop") |>
  left_join(pm_encoding |> select(-n_variants), by = c("role", "species")) |>
  left_join(pm_obs_band |>
              select(species, drug, gm) |>
              pivot_wider(names_from = drug, values_from = gm) |>
              transmute(species, observed_fold_gap = 10^(`Polymyxin B` - `Colistin (polymyxin E)`)),
            by = "species") |>
  mutate(
    ordering = case_when(
      abs(log10(median_predicted_fold_gap)) < log10(1.05) ~ "not resolved",
      is.na(observed_fold_gap) ~ "no observed colistin value",
      sign(log10(median_predicted_fold_gap)) == sign(log10(observed_fold_gap)) ~
        "matches the observed direction",
      TRUE ~ "opposite to the observed direction"))

write_csv(pm_discrimination, file.path(DIR_TABS, sprintf("table29_polymyxin_drug_discrimination%s.csv", PM_OUT_TAG)))
