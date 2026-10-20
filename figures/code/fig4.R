# Figure 4: the trained MuriVerse half-life model applied to colistin (polymyxin E) and
# polymyxin B, two molecules that were never in the cohort, in mouse and in rat.
#
# The figure asks four things of the intravenous comparison. Panels a, b and d: are these
# predictions any better than constants that need no model, where do they sit relative to the
# labels the model was trained on, and which published study each observed value comes from.
# Panel c plots each congener's predicted half-life, in each species, against the intravenous
# observed value for its drug and species, so the departure from identity and the species
# ordering can be read directly. Everything in the figure is the shipped intravenous model
# against intravenous observations.
#
# One model is applied for the shipped predictions: the shipped configuration (annotation +
# study context + termini + physicochemical block read off the PubChem SMILES). The observed
# half-lives come only from the studies flagged in_figure4 in
# MuriVerse/data/polymyxin_studies.csv. Panel c reads
# figures/tabs/table20_polymyxin_predictions.csv (shipped role) and
# table49_polymyxin_observed_by_drug_species_iv.csv. The route-pooled route-shift diagnostic
# (table56_polymyxin_species_route_cells.csv) is still written but is no longer drawn.
if (!exists("ROOT")) {
  source(file.path(system("git rev-parse --show-toplevel", intern = TRUE), "figures/code/00_setup.R"))
}
if (!exists("PM_COMPARATOR_ROUTE")) PM_COMPARATOR_ROUTE <- "IV"   # "all" draws the supplementary pooled version
if (!exists("pm_pred")) source(file.path(DIR_CODE, "06_data_polymyxin.R"))

# Drawn at the 180 mm page width save_fig defaults to, like Figures 2 and 3. It used to be
# drawn at 330 x 430 mm, which print scaled to about 55%, so its 6-8 pt type read as 3-4 pt.
# Type sizes are named once here, as in Figure 3, so no panel drifts below the ~5.5 pt floor.
SZ_AXTITLE <- 6.8
SZ_AXTEXT  <- 6.0
SZ_LEGEND  <- 5.8
SZ_STRIP   <- 6.2
SZ_ANNOT   <- 2.0   # geom_text size in mm, about 5.7 pt

theme_panel <- function(...) {
  theme(axis.title = element_text(size = SZ_AXTITLE),
        axis.text = element_text(size = SZ_AXTEXT),
        strip.text = element_text(size = SZ_STRIP),
        legend.text = element_text(size = SZ_LEGEND),
        legend.title = element_text(size = SZ_LEGEND),
        legend.key.size = unit(6, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(8, 6, 2, 4)) +
    theme(...)
}

# One model is drawn, so a row of the y axis that repeats its long name adds nothing; the
# configuration is named in the caption instead.
NO_ROW_AXIS <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     axis.line.y = element_blank())

# log10-hour axis drawn in hours, so the reader never has to exponentiate
h_axis <- function(breaks = c(0.25, 1, 6, 24, 96), ...) {
  scale_x_continuous(breaks = log10(breaks), labels = hl_lab(breaks), ...)
}

# A fold error of 1.3 is not "1x". Keep a decimal while the number is small enough for one
# to matter, and drop it once the error is large enough that the digit is noise.
fmt_fold <- function(x) ifelse(x < 10, sprintf("%.1f", x), sprintf("%.0f", x))

STUDY_JIAO <- "Jiao2023"
# the observed reference the shaded band shows: intravenous observations (primary) or every route
PM_REF_LABEL <- if (PM_IV) "IV observed ± 1 SD" else "pooled observed ± 1 SD"

# One colour per source of the observed values (Jiao 2023, the Tam lab, everyone else); the
# model is always the red of the other panels.
pal_source <- c(`Jiao 2023` = "#009E73", `Tam lab` = "#332288", `Monash group` = "#E69F00",
                other = "grey25")
SOURCE_LABELS <- c(`Jiao 2023` = "Jiao 2023", `Tam lab` = "Tam lab",
                   `Monash group` = "Monash group", other = "other studies")
COL_MODEL <- pal_role_pm[[1]]
F_DRUG_ORDER <- c("Polymyxin B", "Colistin (polymyxin E)")
DRUG_STRIP <- labeller(drug = c("Colistin (polymyxin E)" = "Colistin", "Polymyxin B" = "Polymyxin B"))

# b: the control - does the model beat constants that need no model? ---------------------
# Every value is drawn once, in hours. Each congener is one mark - the median of its mouse
# and rat predictions, with a bar spanning them - rather than one per variant and species:
# the four congeners are two close pairs, the N-terminal encoding variants are bit-identical
# (literal vs normalised map to the same indicator column), and the species difference is
# about 3%, far inside the 1.46-fold spread across congeners. The
# species contrast has its own panel (d), where it is a difference rather than a position.
# The pooled
# observed values (bar: +/- 1 SD) and the three constants the training split alone would predict
# (the median of the same species, of both species, and of the few lipidated peptides) are each
# one mark too. Left of an observed value the model or constant is too short, right too long.
# The model's marks are the same in panels a and e: pink, one shape per congener, one size.
PM_SHAPES <- c(B1 = 21, B2 = 24, `E1 (colistin A)` = 22, `E2 (colistin B)` = 23)
PM_SPECIES_FILL <- c(Mouse = NA, Rat = COL_MODEL)
PM_MARK_SIZE <- 1.5
# Observed values wear the species colour of the predictions they are compared with, drawn 50%
# transparent (alpha 0.5), so the same hue means "this species" in every panel and the pale mark
# is the observation, the solid one the model.
OBS_ALPHA <- 0.5
# Drawn as the opaque tint an alpha-0.5 mark would show on white, so the mark does not let the grey
# bar or the training density show through it.
OBS_COL <- c(Mouse = mix_col(MV_MOUSE, "#FFFFFF", 1 - OBS_ALPHA), Rat = mix_col(MV_RAT, "#FFFFFF", 1 - OBS_ALPHA))
PM_MARK_ALPHA <- 0.85
B_OBS_LABEL <- \(d, sp) paste0("observed, ", ifelse(d == "Polymyxin B", "polymyxin B", "colistin"), ", ",
                               tolower(sp))
pm_b_const <- pm_baselines |>
  mutate(row = case_when(baseline == "median, same species" ~ paste0("median, ", tolower(species)),
                         TRUE ~ baseline)) |>
  distinct(row, h) |>
  transmute(group = "train median", row, x = h, lo = NA_real_, hi = NA_real_)
pm_b_dat <- bind_rows(
  pm_pred |>
    group_by(congener, species) |>
    summarise(x = median(pred_half_life_h), lo = min(pred_half_life_h), hi = max(pred_half_life_h),
              .groups = "drop") |>
    transmute(group = "model", row = as.character(congener), species = as.character(species), x, lo, hi),
  pm_obs_sum |>
    transmute(group = "observed", row = B_OBS_LABEL(drug, species),
              species = str_to_title(as.character(species)), x = observed_h,
              lo = 10^band_low_log10_h, hi = 10^band_high_log10_h),
  pm_b_const) |>
  mutate(group = factor(group, levels = c("model", "observed", "train median"))) |>
  # Within each block the rows run from the shortest half-life at the top to the longest at the
  # bottom; a discrete y axis starts at the bottom, so the levels are reversed.
  arrange(group, x) |>
  mutate(row = factor(row, levels = rev(unique(as.character(row)))))

# The three groups share one half-life axis and are drawn SUPERIMPOSED rather than in three
# stacked facets. Faceting gave each group its own rows and its own vertical span, so reading
# "is the model where the observed values are" meant carrying a position across two panel
# boundaries. Overlapped, that comparison is one horizontal glance: the question this panel
# exists to answer is whether the pink marks land on the dark ones, and whether either beats
# the grey constants.
#
# Each group gets a small vertical nudge so marks cannot hide one another; the y axis is
# otherwise meaningless and is drawn without labels or ticks. The nudges are deliberately
# smaller than the gap a facet would impose - the groups should read as one cloud compared
# along x, not as three rows.
B_NUDGE <- c(model = 0.30, observed = 0, `train median` = -0.30)
B_GROUP_LABEL <- c(model = "predicted", observed = "observed",
                   `train median` = "train median")

# Overlapping the three facets cost each mark its row label, and in a strip of twelve points
# the shapes alone cannot say WHICH congener, which drug x species, or which constant a mark
# is. So every point carries its own short label and the group names become the y axis.
# Labels are deliberately terse - the group is already given by the axis, so an observed row
# needs only its drug and animal, and a constant only what it is the median of.
b_point_label <- function(group, row) {
  r <- as.character(row)
  dplyr::case_when(
    group == "model" ~ stringr::str_replace(r, " \\(colistin [AB]\\)", ""),
    group == "observed" ~ r |>
      stringr::str_remove("^observed, ") |>
      stringr::str_replace("polymyxin B", "pmx B") |>
      stringr::str_replace(", ", " "),
    TRUE ~ r |>
      stringr::str_remove("^median, ") |>
      stringr::str_replace("both species", "both") |>
      stringr::str_replace(" \\(n = ", " (n="))
}

# The model marks are drawn per species (mouse and rat differ by about 3%), a hair apart in
# height, in the species colours; only the mouse mark carries the congener label, since the
# pair sits on one spot.
pm_b_dat <- pm_b_dat |>
  mutate(ypos = B_NUDGE[as.character(group)] +
           case_when(group == "model" & species == "Mouse" ~ 0.045,
                     group == "model" & species == "Rat" ~ -0.045, TRUE ~ 0),
         point_label = if_else(group == "model" & species == "Rat", "",
                               b_point_label(as.character(group), row)),
         col_key = case_when(group == "model" ~ paste("model", species),
                             group == "observed" ~ paste("observed", species),
                             TRUE ~ as.character(group)),
         lab_key = if_else(group == "observed", paste(col_key, "label"), col_key))

p_b <- ggplot(pm_b_dat, aes(x, ypos)) +
  # A single reference line through the middle, so the eye has something to run along.
  geom_hline(yintercept = 0, colour = "grey92", linewidth = 0.3) +
  geom_linerange(data = \(d) filter(d, group == "model"), aes(xmin = lo, xmax = hi, colour = col_key),
                 linewidth = 0.7, alpha = 0.5) +
  geom_linerange(data = \(d) filter(d, group == "observed"), aes(xmin = lo, xmax = hi),
                 colour = "grey25", linewidth = 0.8) +
  geom_point(data = \(d) filter(d, group == "observed"), aes(colour = col_key), shape = 16, size = 2.2) +
  geom_point(data = \(d) filter(d, group == "train median"), shape = 3, size = 2,
             stroke = 0.8, colour = "grey30") +
  # The model last, so its four congener marks sit on top of everything else.
  geom_point(data = \(d) filter(d, group == "model"), aes(shape = row, colour = col_key, fill = col_key),
             size = PM_MARK_SIZE + 0.6, alpha = PM_MARK_ALPHA) +
  scale_x_continuous(transform = "log10",
                     breaks = c(0.25, 0.5, 1, 2, 4),
                     labels = hl_lab(c(0.25, 0.5, 1, 2, 4)), expand = expansion(mult = 0.08)) +
  # The three group names are the y axis now, rather than text floating in the panel: that
  # frees the space above each mark for its own label and gives the strip a readable left
  # margin. Limits leave room for a label above the topmost band.
  scale_y_continuous(limits = c(-0.46, 0.50), breaks = unname(B_NUDGE),
                     labels = unname(B_GROUP_LABEL), expand = c(0, 0)) +
  scale_shape_manual(values = PM_SHAPES, name = "Congener", limits = names(PM_SHAPES),
                     labels = c("B1", "B2", "E1", "E2")) +
  # One label per point, above its mark. Repelled in x only, so a label never drifts onto a
  # neighbouring group's band and the reader can always tell which row a label belongs to;
  # the tight pairs (observed 1.13/1.21 and 1.90/2.00, and the two training medians at
  # 0.19/0.22) are exactly why this needs a repel rather than a fixed offset.
  ggrepel::geom_text_repel(aes(label = point_label, colour = lab_key),
                           size = 1.5, direction = "x", nudge_y = 0.085,
                           min.segment.length = 0.25, segment.size = 0.18,
                           segment.colour = "grey65", box.padding = 0.12,
                           point.padding = 0.1, max.overlaps = Inf,
                           seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = c(`model Mouse` = MV_MOUSE, `model Rat` = MV_RAT,
                                 `observed Mouse` = unname(OBS_COL["Mouse"]), `observed Rat` = unname(OBS_COL["Rat"]),
                                 `observed Mouse label` = MV_MOUSE, `observed Rat label` = MV_RAT,
                                 `train median` = "grey30"), guide = "none") +
  scale_fill_manual(values = c(`model Mouse` = MV_MOUSE, `model Rat` = MV_RAT), guide = "none") +
  labs(x = "Half-life", y = NULL) +
  theme_panel(legend.position = "none",
              axis.ticks.y = element_blank(),
              axis.text.y = element_text(size = 5.0, hjust = 1),
              panel.grid.major.y = element_blank())

# c: where the predictions sit against the labels the model was trained on ----------------
# One density, one configuration. The key for the three marks is written on the panel.
# One mark per congener, sitting on top of its own rug tick, purely so this panel can carry
# the congener key. The key lives here rather than in panel e's legend block because this
# panel has genuinely empty space above 10 h and panel e's block was four rows deep; the
# shapes are PM_SHAPES, so one key serves panels a, b and e.
pm_a_congener <- pm_pred |>
  group_by(congener, species) |>
  summarise(x = median(pred_log10_half_life_h), .groups = "drop") |>
  # mouse and rat predictions differ by about 3%, so the two species sit at slightly different
  # heights to stay visible; colour is the species colour of Figures 2 and 3
  mutate(y = if_else(species == "Mouse", 0.062, 0.030))

p_c <- ggplot(pm_train_labels, aes(log_label)) +
  geom_density(fill = "grey80", colour = NA, alpha = 0.8) +
  # The rug of all sixteen encoding variants is gone: the four congener marks below are the
  # same predictions, and panel b already carries the variant spread as the range bar, so the
  # ticks were a second drawing of something shown better twice over. The marks sit low in
  # the panel, where the ticks were, so they still read as positions on the half-life axis.
  geom_point(data = pm_a_congener, aes(x = x, y = y, shape = congener, colour = species, fill = species),
             inherit.aes = FALSE, size = PM_MARK_SIZE + 0.3, alpha = PM_MARK_ALPHA) +
  scale_colour_manual(values = c(Mouse = MV_MOUSE, Rat = MV_RAT,
                                 `obs Mouse` = unname(OBS_COL["Mouse"]), `obs Rat` = unname(OBS_COL["Rat"])),
                      breaks = c("Mouse", "Rat"), name = NULL) +
  scale_fill_manual(values = c(Mouse = MV_MOUSE, Rat = MV_RAT), name = NULL) +
  scale_shape_manual(values = PM_SHAPES, name = NULL, limits = names(PM_SHAPES),
                     labels = c("B1", "B2", "E1", "E2")) +
  guides(shape = guide_legend(nrow = 2, byrow = TRUE, order = 1,
                              override.aes = list(colour = "grey35", fill = "grey35")),
         colour = guide_legend(nrow = 1, order = 2, override.aes = list(shape = 16)),
         fill = "none") +
  geom_vline(data = pm_obs_band |> mutate(species = str_to_title(as.character(species))),
             aes(xintercept = gm, colour = paste("obs", species)), linetype = "dotted", linewidth = 0.7,
             show.legend = FALSE) +
  annotate("text", x = log10(1.6), y = 0.47, label = "observed",
           size = SZ_ANNOT, hjust = 0, colour = "grey25") +
  # Lifted clear of the congener marks, which now sit on top of the ticks at y = 0.075; at
  # y = 0.1 this label ran straight through them. "(ticks)" is dropped because the marks
  # above the ticks now say which prediction is which.
  annotate("text", x = log10(90), y = 0.14, label = "model predictions",
           size = SZ_ANNOT, hjust = 1, colour = "grey25") +
  annotate("text", x = log10(0.0035), y = 0.47, label = "training labels",
           size = SZ_ANNOT, hjust = 0, colour = "grey40") +
  h_axis(breaks = c(1 / 60, 5 / 60, 0.25, 1, 6, 24, 96)) +
  coord_cartesian(xlim = log10(c(0.003, 100))) +
  labs(x = "Half-life", y = "Training-label density") +
  # labs(title = "Where the predictions sit in the training label distribution",
  #      subtitle = "ticks = the polymyxin predictions; dotted = observed means") +
  theme_panel(legend.position = "inside", legend.position.inside = c(0.99, 0.99),
              legend.justification = c(1, 1), legend.key.size = unit(5, "pt"), legend.box = "vertical",
              legend.text = element_text(size = 4.8), legend.background = element_blank(),
              legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(0, "pt"))

# d: can the configuration tell the four molecules apart at all? --------------------------
# The count of DISTINCT predictions across the four congeners, holding the encoding variant
# fixed. One means the configuration is blind to the molecule.
p_d <- pm_degeneracy |>
  count(species, role, n_distinct_predictions, name = "variants") |>
  ggplot(aes(n_distinct_predictions, fct_rev(role), colour = role)) +
  geom_vline(xintercept = 4, colour = "grey30", linetype = "dashed", linewidth = 0.5) +
  geom_linerange(aes(xmin = 0, xmax = n_distinct_predictions), linewidth = 2,
                 alpha = MV_FILL_ALPHA) +
  geom_point(size = 2.6) +
  # labels to the RIGHT of the point rather than above it: stacked counts at the same y
  # collided with the lollipop head in the previous version.
  geom_text(aes(label = paste0(n_distinct_predictions, " of 4")), hjust = -0.35,
            size = SZ_ANNOT, colour = "grey20") +
  scale_x_continuous(breaks = 0:4, limits = c(0, 5.4)) +
  scale_colour_manual(values = pal_role_pm, guide = "none") +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  facet_wrap(~species, nrow = 1) +
  coord_cartesian(clip = "off") +
  labs(x = "Distinct predictions across the four congeners", y = NULL) +
  # labs(title = "Can it tell the molecules apart?",
  #      subtitle = sprintf("4 congeners = 4 predictions; %d encoding variants per species",
  #                         n_distinct(pm_degeneracy$variant))) +
  theme_panel(plot.margin = margin(6, 10, 2, 4)) + NO_ROW_AXIS

# e: what moves the prediction - the drug, or an arbitrary transcription choice? ----------
pm_sens <- pm_discrimination |>
  transmute(species, role,
            `input encoding choice` = max_encoding_fold_spread,
            `colistin vs polymyxin B`   = max_abs_predicted_fold_gap) |>
  pivot_longer(c(`input encoding choice`, `colistin vs polymyxin B`),
               names_to = "source", values_to = "fold") |>
  mutate(source = factor(source, levels = c("input encoding choice",
                                            "colistin vs polymyxin B")))

obs_gap <- pm_discrimination |> distinct(species, observed_fold_gap)

# Lollipops from 1 rather than bars from 0: these are fold changes, and 1 is the origin that
# means "this moves nothing". The dashed line is the observed polymyxin B / colistin ratio.
p_e <- ggplot(pm_sens, aes(y = fct_rev(role), colour = source)) +
  geom_vline(data = obs_gap, aes(xintercept = pmax(observed_fold_gap, 1 / observed_fold_gap)),
             colour = "grey25", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = 1, colour = "grey50", linewidth = 0.4) +
  geom_linerange(aes(xmin = 1, xmax = fold), position = position_dodge(width = 0.55),
                 linewidth = 2, alpha = 0.9) +
  geom_point(aes(x = fold), position = position_dodge(width = 0.55), size = 2) +
  scale_colour_manual(values = c("input encoding choice" = "grey60",
                                 "colistin vs polymyxin B" = "#332288"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  facet_wrap(~species, nrow = 1) +
  coord_cartesian(clip = "off") +
  labs(x = "Fold change in predicted half-life", y = NULL) +
  # labs(title = "What actually moves it",
  #      subtitle = "wording vs molecule; dashed = observed gap") +
  theme_panel(legend.position = "bottom", legend.box.margin = margin(-4, 0, 0, 0),
              plot.margin = margin(6, 10, 2, 4)) + NO_ROW_AXIS

# a: every observed value, by study ---------------------------
# One row per study and drug, grouped by who produced the value: Jiao 2023 first, then the Tam
# lab, then everyone else, polymyxin B above colistin. Colour is the source group; the text to
# the right says how that study calculated its half-life, because the polymyxin B values differ
# about as much by method as by animal. A study with no half-life we can use keeps its row.
pm_f_used <- pm_obs |> distinct(study_id, drug, species)
# A study is "without a usable value" when it has none on ANY route (pm_obs_any): under the
# intravenous comparison a study whose usable values are all subcutaneous is simply not drawn.
pm_f_unused <- pm_studies |>
  filter(!study_id %in% pm_obs_any$study_id, drug != "both") |>
  transmute(study_id, drug = factor(drug, levels = levels(pm_pred$drug)),
            species = factor(str_to_title(species), levels = levels(pm_obs$species)))

# The route and dose behind each row, so a reader can see where dose differs between studies
fmt_dose <- function(lo, hi) {
  f <- \(x) format(x, trim = TRUE, drop0trailing = TRUE)
  ifelse(lo == hi, f(lo), paste0(f(lo), "–", f(hi)))
}
pm_f_dose <- pm_obs |>
  group_by(study_id, drug) |>
  summarise(route_class = first(route_class),
            lo = if (all(is.na(dose_lo_mg_kg))) NA_real_ else min(dose_lo_mg_kg, na.rm = TRUE),
            hi = if (all(is.na(dose_hi_mg_kg))) NA_real_ else max(dose_hi_mg_kg, na.rm = TRUE),
            .groups = "drop") |>
  mutate(dose_lab = ifelse(is.na(lo), "dose not reported",
                           paste0(route_class, " ", fmt_dose(lo, hi), " mg/kg"))) |>
  select(study_id, drug, dose_lab)

pm_f_mean <- pm_obs |>
  group_by(study_id, drug) |>
  summarise(mean_log10 = mean(log10_half_life_h), .groups = "drop")

pm_f_rows <- bind_rows(pm_f_used, pm_f_unused) |>
  left_join(pm_studies |> select(study_id, label, year, source_group, half_life_type, half_life_analysis,
           animal_state, matrix_basis, sampling_window_h),
            by = "study_id") |>
  left_join(pm_f_dose, by = c("study_id", "drug")) |>
  left_join(pm_f_mean, by = c("study_id", "drug")) |>
  mutate(unused = study_id %in% pm_f_unused$study_id,
         row_label = paste0(label, if_else(unused, " (no usable value)",
                                           paste0(" · ", dose_lab))),
         key = paste(drug, species, row_label),
         # ascending by the study's mean half-life, lowest at the top; a study with no usable
         # value goes last
         key = fct_reorder(key, -replace_na(mean_log10, 99)),
         half_life_type = replace_na(half_life_type, ""),
         half_life_analysis = replace_na(half_life_analysis, ""),
         # what else differs between studies and could move a half-life: the animals' state, the
         # matrix and how long they were sampled (blank where the text we had did not say)
         # a missing part says so, instead of blanking the whole label
         window_lab = case_when(is.na(sampling_window_h) ~ "not read",
                                sampling_window_h == "not stated" ~ "not stated",
                                TRUE ~ paste0("to ", sampling_window_h, " h")),
         context = paste(animal_state, matrix_basis, window_lab, sep = " · "),
         # the type of half-life, then the analysis that gave it
         half_life_label = if_else(half_life_analysis == "", paste(half_life_type, "not read", sep = " · "),
                                   paste(half_life_type, half_life_analysis, sep = " · ")),
         drug = factor(as.character(drug), levels = F_DRUG_ORDER),
         # one facet per drug and species, so the rows are grouped by animal
         panel = factor(paste(drug, species, sep = " · "),
                        levels = paste(rep(F_DRUG_ORDER, each = 2), c("Mouse", "Rat"), sep = " · ")))
pm_f_obs <- bind_rows(pm_obs, pm_mrt) |>
  mutate(drug = factor(as.character(drug), levels = F_DRUG_ORDER)) |>
  left_join(pm_f_rows |> select(study_id, drug, species, key, panel),
            by = c("study_id", "drug", "species"))
F_METHOD_X <- log10(15) + 0.12
F_CONTEXT_X <- log10(15) + 1.65
F_XLIM <- c(log10(0.17), log10(15) + 2.8)

# The pooled observed value that panels a, b and f compare the model with (animal-weighted mean of the
# dose groups, +/- 1 SD), and the same with every study counted once.
pm_f_ref <- pm_obs_sum |>
  mutate(drug = factor(as.character(drug), levels = F_DRUG_ORDER),
         panel = factor(paste(drug, species, sep = " · "),
                        levels = paste(rep(F_DRUG_ORDER, each = 2), c("Mouse", "Rat"), sep = " · ")),
         lo = band_low_log10_h, hi = band_high_log10_h, mid = observed_log10_h)
pm_f_once <- pm_by_source |> filter(source_group == "pooled: each study once") |>
  mutate(drug = factor(as.character(drug), levels = F_DRUG_ORDER),
         panel = factor(paste(drug, species, sep = " · "),
                        levels = paste(rep(F_DRUG_ORDER, each = 2), c("Mouse", "Rat"), sep = " · ")),
         mid = log10(observed_h))

pm_f_pred <- pm_pred |>
  mutate(drug = factor(as.character(drug), levels = F_DRUG_ORDER),
         species = factor(as.character(species), levels = c("Mouse", "Rat"))) |>
  # Per drug AND species. Species is a model input, so the two facets of a drug no longer
  # share one band; drawing the pooled band in both would hide the (small) difference the
  # model does apply.
  group_by(drug, species) |>
  summarise(lo = min(pred_log10_half_life_h), hi = max(pred_log10_half_life_h),
            med = median(pred_log10_half_life_h), .groups = "drop") |>
  mutate(panel = factor(paste(drug, species, sep = " · "),
                        levels = paste(rep(F_DRUG_ORDER, each = 2), c("Mouse", "Rat"), sep = " · ")))

# Each prediction (every variant of every congener) as a point on a "model" row at the top of
# its block, so the individual values sit under the band that summarises them. Each species'
# block carries its own predictions, since the two are no longer the same numbers.
MODEL_KEYS <- paste(rep(F_DRUG_ORDER, each = 2), c("Mouse", "Rat"), "model")
KEY_LEVELS <- c(levels(pm_f_rows$key), MODEL_KEYS)
pm_f_rows <- pm_f_rows |> mutate(key = factor(as.character(key), levels = KEY_LEVELS))
pm_f_obs <- pm_f_obs |> mutate(key = factor(as.character(key), levels = KEY_LEVELS))
pm_f_model <- pm_pred |>
  mutate(drug = factor(as.character(drug), levels = F_DRUG_ORDER),
         species = factor(str_to_title(as.character(species)), levels = c("Mouse", "Rat")),
         panel = factor(paste(drug, species, sep = " · "),
                        levels = paste(rep(F_DRUG_ORDER, each = 2), c("Mouse", "Rat"), sep = " · ")),
         key = factor(paste(drug, species, "model"), levels = KEY_LEVELS))

p_a <- ggplot(pm_f_rows, aes(y = key)) +
  # first, so every row is trained in its table order whatever the highlight layers hold
  geom_blank() +
  # the model's predictions (all variants), shaded from the lowest to the highest, with the
  # median as a line. One band per drug AND species, since species is a model input.
  geom_rect(data = pm_f_ref, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf,
                                 fill = PM_REF_LABEL), alpha = 0.22, inherit.aes = FALSE) +
  geom_vline(data = pm_f_ref, aes(xintercept = mid), colour = "grey35", linewidth = 0.5) +
  # the same pooled value with each study counted once
  geom_vline(data = pm_f_once, aes(xintercept = mid), colour = "grey35", linewidth = 0.5,
             linetype = "dotted") +
  geom_rect(data = pm_f_pred, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf,
                                  fill = paste0("model, ", tolower(as.character(species)))),
            alpha = 0.15, inherit.aes = FALSE) +
  geom_vline(data = pm_f_pred, aes(xintercept = med, colour = paste("model", tolower(as.character(species)))),
             linewidth = 0.6) +
  geom_point(data = pm_f_model, aes(pred_log10_half_life_h, key, shape = congener,
                                    colour = paste("model", tolower(as.character(species))),
                                    fill = paste0("model, ", tolower(as.character(species)))),
             size = PM_MARK_SIZE, alpha = PM_MARK_ALPHA,
             inherit.aes = FALSE, position = position_jitter(height = 0.18, width = 0, seed = 1986)) +
  # Jiao 2023, whose terminal half-lives are about 4 h, in its source colour
  geom_tile(data = pm_f_rows |> filter(study_id == STUDY_JIAO),
            aes(x = mean(F_XLIM), y = key), width = diff(F_XLIM), height = 0.92,
            fill = pal_source[["Jiao 2023"]], alpha = 0.16, inherit.aes = FALSE) +
  geom_text(data = pm_f_obs |> filter(study_id == STUDY_JIAO) |> group_by(drug, species, panel) |>
              summarise(x = min(log10_half_life_h), key = first(key),
                        label = paste0("t½ ", paste(format(sort(half_life_h), digits = 3),
                                                   collapse = " and "), " h"),
                        .groups = "drop"),
            aes(x - 0.08, key, label = label), colour = pal_source[["Jiao 2023"]],
            fontface = "bold", hjust = 1, size = SZ_ANNOT) +
  geom_text(aes(x = F_CONTEXT_X, label = context), hjust = 0, size = SZ_ANNOT,
            colour = "grey45") +
  # headings of the two text columns, above the first block only
  geom_text(data = pm_f_rows |> filter(panel == levels(panel)[1]) |> slice(1),
            aes(x = F_METHOD_X, y = Inf, label = "half-life · analysis"), hjust = 0, vjust = -0.8,
            size = SZ_ANNOT, colour = "grey25", fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = pm_f_rows |> filter(panel == levels(panel)[1]) |> slice(1),
            aes(x = F_CONTEXT_X, y = Inf, label = "animals · matrix · sampled"), hjust = 0, vjust = -0.8,
            size = SZ_ANNOT, colour = "grey25", fontface = "bold", inherit.aes = FALSE) +
  geom_text(aes(x = F_METHOD_X, label = half_life_label), hjust = 0, size = SZ_ANNOT,
            colour = "grey25") +
  # MRTs are not half-lives and are not drawn in panel e; panel g keeps them (see pm_h_mrt)
  geom_point(data = pm_f_obs |> filter(value_type != "mean residence time"),
             aes(log10_half_life_h, key, colour = source_group, shape = species,
                 alpha = value_type),
             stroke = 0.8, size = 1.8) +
  # Jiao 2023 is labelled on its own row (shaded green), so it is left out of the key
  scale_colour_manual(values = c(pal_source, `model mouse` = MV_MOUSE, `model rat` = MV_RAT), name = NULL,
                      breaks = setdiff(names(pal_source), "Jiao 2023"),
                      labels = SOURCE_LABELS[setdiff(names(pal_source), "Jiao 2023")]) +
  scale_fill_manual(values = c(`model, mouse` = MV_MOUSE, `model, rat` = MV_RAT,
                               setNames("grey45", PM_REF_LABEL)),
                    breaks = c(PM_REF_LABEL, "model, mouse", "model, rat"),
                    name = "Reference") +
  scale_shape_manual(values = c(Mouse = 4, Rat = 1, PM_SHAPES), name = NULL,
                     breaks = names(PM_SHAPES)) +
  # Derived values (we computed ln2 x V / CL) are drawn paler: they are not reported numbers.
  scale_alpha_manual(values = c(reported = 1, derived = 0.5, `mean residence time` = 1),
                     name = NULL, breaks = c("reported", "derived"),
                     labels = c(reported = "reported", derived = "derived by us")) +
  scale_y_discrete(labels = \(k) ifelse(str_detect(k, " model$"), "model",
                                         str_remove(k, "^(Colistin \\(polymyxin E\\)|Polymyxin B) (Mouse|Rat) ")),
                   drop = TRUE, expand = expansion(add = 1.0)) +
  scale_x_continuous(breaks = log10(c(0.25, 0.5, 1, 2, 4, 6)),
                     labels = hl_lab(c(0.25, 0.5, 1, 2, 4, 6)),
                     expand = expansion(mult = c(0.03, 0))) +
  # drug outside, species inside
  facet_grid(drug + species ~ ., scales = "free_y", space = "free_y", labeller = DRUG_STRIP) +
  coord_cartesian(xlim = F_XLIM, clip = "off") +
  labs(x = "Observed half-life", y = NULL) +
  # labs(title = "Which study each value comes from",
  #      subtitle = "plasma or serum t1/2; polymyxin B studies above, colistin below") +
  # The congener key moved to panel a's top-right whitespace, so it is suppressed here: the
  # shapes are the same PM_SHAPES in both panels and one key is enough for the figure.
  guides(colour = guide_legend(order = 1, nrow = 1, title = "Observed", override.aes = list(shape = 4)),
         shape = "none",
         fill = guide_legend(order = 3, nrow = 1, title = "Reference",
                             override.aes = list(alpha = 0.3, colour = NA)),
         alpha = guide_legend(order = 4, nrow = 1, title = "Value",
                              override.aes = list(shape = 4, colour = "grey25"))) +
  theme_panel(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
              legend.justification = "left", legend.spacing.y = unit(0, "pt"),
              legend.title = element_text(face = "bold"),
              strip.text.y = element_text(size = 5.2),
              plot.margin = margin(16, 6, 2, 4),
              legend.box.margin = margin(-4, 0, 0, 0),
              panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))

# g: the model against each source, as a signed fold ----------------------------------------
# The same 16 predictions per drug and species, each divided by the observed value of one
# source. Right of 1 the model is too long, left of 1 too short. Dots are the median variant;
# the thick line spans the variants and the thin line is the observed value's own +/- 1 SD, so
# the two kinds of uncertainty are not confused. Values the studies report are kept apart from
# values we derived (paler) because the two are not the same quantity. A mean residence time is
# not a half-life and is not compared with the model here (it is drawn in panel a). The last
# four rows are the pooled value of panels a, b and f and three other ways to pool.
G_SOURCES <- c("Jiao 2023", "Tam lab", "Monash group", "other")
G_LAB <- c(`Jiao 2023` = "Jiao 2023", `Tam lab` = "Tam lab", `Monash group` = "Monash group",
           other = "other studies")
G_POOL <- c(pooled_primary = "all studies pooled",
            `pooled: each study once` = "pooled, each study once",
            `pooled: reported only` = "pooled, reported only",
            `pooled: IV only` = "pooled, IV only")
G_ROW_LEVELS <- c(unlist(lapply(G_SOURCES[-1], \(g) paste0(G_LAB[[g]], ", ", c("reported", "derived")))),
                  unname(G_POOL))                              # Jiao 2023 is added at the top
G_ROW_LEVELS <- c("Jiao 2023", G_ROW_LEVELS)
pm_g_primary <- pm_pred |>
  left_join(pm_obs_sum |> transmute(drug, species, observed_h, n_studies, n_animals,
                                    sd_log10_h = band_high_log10_h - observed_log10_h),
            by = c("drug", "species")) |>
  transmute(drug, species, source_group = "pooled_primary", value_basis = "",
            n_studies, n_animals, sd_log10_h, pred_over_observed = pred_half_life_h / observed_h)
pm_g_dat <- bind_rows(
  pm_g_primary,
  pm_fold_src |> filter(quantity == "elimination half-life") |>
    transmute(drug, species, source_group = as.character(source_group),
              value_basis = replace_na(value_basis, ""), n_studies, n_animals, sd_log10_h,
              pred_over_observed)) |>
  mutate(label = case_when(source_group %in% names(G_POOL) ~ unname(G_POOL[source_group]),
                           source_group == "Jiao 2023" ~ "Jiao 2023",
                           TRUE ~ paste0(G_LAB[source_group], ", ", value_basis)),
         label = factor(label, levels = G_ROW_LEVELS)) |>
  mutate(value_basis = if_else(value_basis == "", "reported", value_basis)) |>
  group_by(drug, species, label, source_group, value_basis) |>
  summarise(med = median(pred_over_observed), lo = min(pred_over_observed),
            hi = max(pred_over_observed), sd = first(sd_log10_h), n_an = first(n_animals),
            n_st = first(n_studies), .groups = "drop") |>
  mutate(sd_lo = if_else(is.na(sd), NA_real_, med / 10^sd),
         sd_hi = if_else(is.na(sd), NA_real_, med * 10^sd),
         right = pmax(hi, sd_hi, na.rm = TRUE),
         lab = paste0(ifelse(med >= 1, paste0(fmt_fold(med), "× over"),
                             paste0(fmt_fold(1 / med), "× under")),
                      " · n = ", n_an, " (", n_st, ")"),
         drug = factor(as.character(drug), levels = F_DRUG_ORDER),
         colour_key = case_when(source_group == "Jiao 2023" ~ "Jiao 2023",
                                source_group == "Tam lab" ~ "Tam lab",
                                source_group == "Monash group" ~ "Monash group",
                                source_group == "other" ~ "other",
                                TRUE ~ "pooled"),
         key = factor(paste(drug, label), levels = rev(paste(rep(F_DRUG_ORDER, each = length(G_ROW_LEVELS)),
                                                              rep(G_ROW_LEVELS, 2)))))

p_g <- ggplot(pm_g_dat, aes(y = key, colour = colour_key, alpha = value_basis)) +
  geom_vline(xintercept = 1, colour = "grey50", linewidth = 0.4) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.8) +
  geom_linerange(aes(xmin = sd_lo, xmax = sd_hi), linewidth = 0.3, na.rm = TRUE) +
  geom_point(aes(x = med), size = 2.2) +
  geom_text(aes(x = right, label = lab), hjust = -0.1, size = SZ_ANNOT, colour = "grey20",
            alpha = 1) +
  scale_colour_manual(values = c(pal_source, pooled = "black"), guide = "none") +
  scale_alpha_manual(values = c(reported = 1, derived = 0.5), guide = "none") +
  scale_x_continuous(trans = "log10", breaks = c(0.5, 1, 2, 4, 8, 16, 32),
                     labels = c("0.5", "1", "2", "4", "8", "16", "32"),
                     expand = expansion(mult = c(0.05, 0.75))) +
  scale_y_discrete(labels = \(k) str_remove(k, "^(Colistin \\(polymyxin E\\)|Polymyxin B) "),
                   drop = TRUE) +
  facet_grid(drug ~ species, scales = "free_y", space = "free_y", labeller = DRUG_STRIP) +
  labs(x = paste("Predicted / observed half-life (1 = exact; above 1 the model is too long).",
                 "n = animals (studies)"), y = NULL) +
  theme_panel(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))

# h: does dose move the observed half-life, or the model? -------------------------------------
# Each observed half-life at the dose it was measured at (a range as a bar), by route. The red
# line is the shipped model's median prediction with the dose supplied: the dose-missing
# prediction of the other panels (which the model treats as 1 mg/kg) scaled by the model's own
# response to dose (table54), and the band spans all variants. Only 3 of the 320 training rows
# carry a dose, so that response is an extrapolation, not a learned dose-response.
pm_h_obs <- pm_obs |> filter(!is.na(dose_lo_mg_kg)) |>
  mutate(panel = interaction(species, drug, sep = " · "))
pm_h_mrt <- pm_mrt |> filter(!is.na(dose_lo_mg_kg)) |>
  mutate(panel = interaction(species, drug, sep = " · "))
pm_h_model <- pm_pred |>
  group_by(drug) |>
  summarise(med = median(pred_half_life_h), lo = min(pred_half_life_h), hi = max(pred_half_life_h),
            .groups = "drop") |>
  inner_join(pm_dose_fold, by = "drug") |>
  cross_join(tibble(species = factor(PM_SPECIES, levels = PM_SPECIES))) |>
  mutate(panel = interaction(species, drug, sep = " · "))
H_LEVELS <- c("Mouse · Polymyxin B", "Rat · Polymyxin B",
              "Mouse · Colistin (polymyxin E)", "Rat · Colistin (polymyxin E)")
H_LABELS <- c("Polymyxin B, mouse", "Polymyxin B, rat", "Colistin, mouse", "Colistin, rat")
for (nm in c("pm_h_obs", "pm_h_mrt", "pm_h_model")) {
  assign(nm, get(nm) |> mutate(panel = factor(as.character(panel), levels = H_LEVELS,
                                             labels = H_LABELS)))
}

p_h <- ggplot() +
  geom_ribbon(data = pm_h_model, aes(dose_mg_kg, ymin = lo * fold, ymax = hi * fold),
              fill = COL_MODEL, alpha = 0.10) +
  geom_line(data = pm_h_model, aes(dose_mg_kg, med * fold), colour = COL_MODEL, linewidth = 0.6,
            linetype = "dashed") +
  # the dose response is not learned: only 3 of the 320 training records carry a dose
  geom_text(data = tibble(panel = factor(H_LABELS[2], levels = H_LABELS)),
            aes(x = 95, y = 0.12, label = "model line: extrapolated\n(3 dosed training records)"),
            hjust = 1, vjust = 0, size = SZ_ANNOT, colour = COL_MODEL, lineheight = 0.9) +
  geom_segment(data = pm_h_obs |> filter(dose_lo_mg_kg < dose_hi_mg_kg),
               aes(x = dose_lo_mg_kg, xend = dose_hi_mg_kg, y = half_life_h, yend = half_life_h,
                   colour = source_group), linewidth = 0.5, alpha = 0.6) +
  geom_point(data = pm_h_obs,
             aes(sqrt(dose_lo_mg_kg * dose_hi_mg_kg), half_life_h, colour = source_group,
                 shape = route_class, alpha = value_type), size = 1.7) +
  geom_point(data = pm_h_mrt, aes(dose_lo_mg_kg, half_life_h), shape = 8, colour = "grey25",
             size = 1.9, stroke = 0.8) +
  scale_colour_manual(values = pal_source, labels = SOURCE_LABELS, name = NULL, guide = "none") +
  scale_shape_manual(values = c(IV = 16, SC = 17, IP = 15), name = "route") +
  scale_alpha_manual(values = c(reported = 1, derived = 0.5), guide = "none") +
  scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100"),
                limits = c(0.1, 100)) +
  scale_y_log10(breaks = c(0.3, 1, 3, 10), labels = c("0.3", "1", "3", "10")) +
  facet_wrap(~panel, nrow = 1) +
  labs(x = "Dose (mg/kg, as reported)", y = "Half-life (h)") +
  # labs(title = "Does dose move the observed half-life or the model?") +
  theme_panel(legend.position = "bottom", legend.box.margin = margin(-4, 0, 0, 0))

# Page-width layout: (a, b) share one row, then (c) and (d) each span the full width.
# The degeneracy and encoding/drug fold panels of earlier builds are no longer drawn here;
# their underlying tables (table23 and table29) are still written by predict_polymyxins.py.
LAB <- list(size = 9, face = "bold")
# --- c: predicted against observed, shipped intravenous model, rat and mouse overlaid ---------
# One mark per congener in each species (8 distinct predictions; the literal and normalised
# N-terminal wordings are identical), placed at the intravenous observed value of its drug and
# species, so the two congeners of a drug share an x position. Rat and mouse are separate facets
# overlaid in one panel, so the identity line (dashed) is shared.
# Congeners use the PM_SHAPES of panel a, so one key serves a and c. Colour is the species colour of
# Figures 2, 3 and supplementary 2 (MV_RAT, MV_MOUSE), not the model colour.
sr_obs <- read_csv(file.path(DIR_TABS, "table49_polymyxin_observed_by_drug_species_iv.csv"),
                   show_col_types = FALSE) |>
  select(drug, species, obs_log10 = observed_log10_h)
sr_cells <- read_csv(file.path(DIR_TABS, "table20_polymyxin_predictions.csv"), show_col_types = FALSE) |>
  filter(role == "shipped") |>
  group_by(species, drug, name) |>
  summarise(pred_log10 = mean(pred_log10_half_life_h), .groups = "drop") |>
  inner_join(sr_obs, by = c("drug", "species")) |>
  mutate(fold_error = 10^abs(pred_log10 - obs_log10),
         congener = factor(CONGENER_SHORT[name], levels = unname(CONGENER_SHORT)),
         species = factor(str_to_title(species), levels = c("Rat", "Mouse")))
stopifnot(nrow(sr_cells) == 8)
sr_err <- sr_cells |> group_by(species) |>
  summarise(lab = sprintf("%s median error %.1f\u00d7", tolower(first(species)), median(fold_error)), .groups = "drop")

sr_lim <- range(c(sr_cells$obs_log10, sr_cells$pred_log10)) + c(-0.12, 0.12)
p_i <- ggplot(sr_cells, aes(obs_log10, pred_log10, shape = congener,
                            colour = species, fill = species)) +
  geom_abline(linetype = "dashed", colour = "grey60") +
  geom_point(size = PM_MARK_SIZE + 0.3, alpha = PM_MARK_ALPHA) +
  scale_colour_manual(values = c(Rat = MV_RAT, Mouse = MV_MOUSE), name = NULL) +
  scale_fill_manual(values = c(Rat = MV_RAT, Mouse = MV_MOUSE), name = NULL) +
  scale_shape_manual(values = PM_SHAPES, name = NULL, limits = names(PM_SHAPES),
                     labels = c("B1", "B2", "E1", "E2")) +
  guides(shape = guide_legend(nrow = 1, order = 1,
                              override.aes = list(colour = "grey35", fill = "grey35")),
         colour = guide_legend(nrow = 1, order = 2, override.aes = list(shape = 16)),
         fill = "none") +
  geom_text(data = sr_err, aes(label = lab, colour = species), x = Inf, y = -Inf,
            vjust = -0.6 - 1.4 * (seq_len(nrow(sr_err)) - 1), hjust = 1.06,
            inherit.aes = FALSE, size = 1.8, show.legend = FALSE) +
  coord_equal(xlim = sr_lim, ylim = sr_lim, expand = FALSE) +
  scale_x_loghl(breaks = c(0.5, 1, 2)) + scale_y_loghl(breaks = c(0.5, 1, 2)) +
  labs(x = "Observed half-life", y = "Predicted half-life") +
  theme_panel(legend.position = "top", legend.key.size = unit(5, "pt"),
              legend.text = element_text(size = 4.8),
              legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
              legend.spacing.y = unit(0, "pt"),
              legend.box = "vertical")

# Row 1 is a density (p_c) stacked over the superimposed strip (p_b) on the left, and the
# predicted-vs-observed scatter (p_i) on the right. Row 2 is the full-width study forest
# (p_a). The
# per-source (p_g), dose (p_h), degeneracy (p_d) and sensitivity (p_e) panels are built
# above and left out of the layout.
if (PM_IV) {
  fig4 <- ggarrange(
    ggarrange(
      # b is a one-dimensional strip and needs less height than the density above it.
      ggarrange(p_c, p_b, nrow = 2, heights = c(1.15, 0.85),
                labels = c("a", "b"), font.label = LAB),
      ggarrange(p_i, labels = "c", font.label = LAB),
      # Even split: a/b and c each take half the page width.
      ncol = 2, widths = c(1, 1)),
    ggarrange(p_a, labels = "d", font.label = LAB),
    nrow = 2, heights = c(88, 91)
  )
  save_fig(fig4, "fig4", width = 180, height = 171)
} else {
  # Supplementary: the same comparison against observations from EVERY route. Panel c of
  # Figure 4 is not repeated; it never depended on the comparator.
  # a and b sit side by side on the first row; the study forest takes the full width below.
  supp_pooled <- ggarrange(
    ggarrange(p_c, p_b, ncol = 2, widths = c(1, 1),
              labels = c("a", "b"), font.label = LAB),
    ggarrange(p_a, labels = "c", font.label = LAB),
    nrow = 2, heights = c(62, 110)
  )
  save_fig(supp_pooled, "supp_fig3", width = 180, height = 180)
}
