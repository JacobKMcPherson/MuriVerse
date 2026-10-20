# Supplementary Figure 1 (11 panels, 4 rows): source databases → 320-record
# modelling cohort, then half-life distributions by structural and study class.
#
# Row 1 (a-c): species composition of source databases; PEPlife2 by setting;
#              PEPlife2 half-life by species.
# Row 2 (d-e): merged build-funnel two-facet panel (records remaining AND
#               removed at each stage) in cols 1-2, species shift in col 3.
# Row 3 (f-h): three reconciliation gates (route, structure, assay) 467 → 320.
# Row 4 (i-k): half-life distributions across Engineering, Topology and Chirality
#              classes in the 176-record intravenous training split, showing both
#              the class composition and the within-class species split.
#
# Self-contained: sources only 00_setup.R and 02_data_preprocessing.R.
if (!exists("ROOT")) {
  source(file.path(system("git rev-parse --show-toplevel", intern = TRUE),
                   "figures/code/00_setup.R"))
}
if (!exists("hl_sources")) source(file.path(DIR_CODE, "02_data_preprocessing.R"))

# ---------- Load full modelling-set metadata (row 4) ----------
modelling_full <- read_csv(
  file.path(ROOT, "MuriVerse", "rodent", "strict", "halflife_meta_with_split.csv"),
  show_col_types = FALSE, guess_max = 1e5
) |>
  select(ID = peplife_id, species, route_class, split, hl_engineering,
         topology, chirality, half_life_h) |>
  mutate(
    species = str_to_title(species),
    Engineering = case_when(
      hl_engineering == "unmodified / none stated"     ~ "Unmodified",
      hl_engineering == "backbone modification only"   ~ "Backbone mod.",
      hl_engineering == "radio/fluorescent label only" ~ "Radio / fluorescent",
      hl_engineering == "lipidated"                    ~ "Lipidated",
      TRUE ~ hl_engineering
    ) |> factor(levels = c("Unmodified", "Backbone mod.",
                            "Radio / fluorescent", "Lipidated")),
    Topology  = str_to_title(topology) |>
      factor(levels = c("Linear", "Cyclic")),
    Chirality = case_when(
      chirality == "l"   ~ "L",
      chirality == "d"   ~ "D",
      chirality == "mix" ~ "Mixed",
      TRUE ~ chirality
    ) |> factor(levels = c("L", "D", "Mixed")),
    Route = str_to_upper(route_class) |>
      factor(levels = c("IV", "SC", "IP", "IM"))
  )

pal_species_mv <- c(Mouse = MV_MOUSE, Rat = MV_RAT)

# Row 4 describes the intravenous training split (176 of the 213 intravenous records), the rows
# the polymyxin baseline constants are read from; every record is intravenous, so there is no
# route panel.
iv_train <- modelling_full |> filter(route_class == "iv", split == "train")
stopifnot(nrow(iv_train) == 176)

# ---------- Palettes (row 1) ----------
sp_levels   <- levels(hl_sources$species_grp)
n_animal    <- sum(!sp_levels %in% non_animal)
pal_species_db <- setNames(
  c(ggpubr::get_palette("npg", n_animal),
    c("#B09C85", "#D9CBB8", "grey70", "grey45")[seq_len(length(sp_levels) - n_animal)]),
  sp_levels
)

PAL_KEPT  <- c(`kept by gate`    = COL_ANALYSIS,
               `dropped by gate` = "grey70")

PAL_BLOCK <- c(`PEPlife2 → eligible rodent set`                   = "#88CCEE",
               `eligible rodent set → structurally verified` = COL_MV_GREEN)

wrap_labs <- function(p, title_w = 24, sub_w = 32) {
  for (f in c("title", "subtitle")) {
    v <- p$labels[[f]]
    if (is.character(v) && length(v) == 1)
      p$labels[[f]] <- paste(strwrap(v, width = if (f == "title") title_w else sub_w),
                             collapse = "\n")
  }
  p + theme(plot.title.position = "plot",
            plot.title    = element_text(size = 10, face = "plain"),
            plot.subtitle = element_text(size = 8))
}

# =============================================================================
# Row 1: source databases
# =============================================================================
db_n <- count(hl_sources, database)

p_a <- ggplot(hl_sources, aes(y = fct_rev(database), fill = species_grp)) +
  geom_bar(position = position_fill(reverse = TRUE), width = 0.7) +
  geom_text(data = db_n, aes(x = 1.02, label = scales::comma(n), fill = NULL),
            hjust = 0, size = 2.6) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, 0.5),
                     expand = expansion(mult = c(0, 0.32))) +
  scale_fill_manual(values = pal_species_db) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
  labs(x = "Share of half-life records (n beside bar)", y = NULL, fill = "Species",
       title = "Species by database") +
  theme(legend.text    = element_text(size = 7),
        legend.key.size = unit(3, "mm"))

pep <- filter(hl_sources, database == "PEPlife2")

p_b <- ggplot(pep, aes(y = fct_rev(species_grp), fill = setting)) +
  geom_bar(width = 0.7) +
  scale_x_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = PAL_SETTING) +
  labs(x = "PEPlife2 records", y = NULL, fill = "Setting",
       title = "PEPlife2: species and setting")

pep_hl <- pep |>
  filter(half_life_h > 0, !species_grp %in% c("unrecorded", "unclassified"))

p_c <- ggplot(pep_hl, aes(fct_rev(species_grp), half_life_h, fill = species_grp)) +
  geom_violin(scale = "width", linewidth = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  scale_y_hl(breaks = c(1 / 60, 5 / 60, 6, 72)) +
  scale_fill_manual(values = pal_species_db, guide = "none") +
  coord_flip() +
  labs(x = NULL, y = "Half-life",
       title = "PEPlife2: half-life by species")

# =============================================================================
# Row 2: build pipeline — merged two-facet panel (d) + species shift (e)
# =============================================================================
pipeline_short <- build_pipeline |>
  mutate(stage_label = fct_inorder(as.character(stage)))

# Merge panels d and e as a two-facet panel: "Records remaining" | "Records removed"
pipeline_de <- pipeline_short |>
  select(stage_label, block, records, removed) |>
  pivot_longer(c(records, removed),
               names_to  = "metric",
               values_to = "count") |>
  mutate(metric = factor(
    ifelse(metric == "records", "Records remaining", "Records removed"),
    levels = c("Records remaining", "Records removed")
  ))

p_d <- ggplot(pipeline_de, aes(count, fct_rev(stage_label), fill = block)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = scales::comma(count)), hjust = -0.12, size = 1.9) +
  scale_x_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    expand = expansion(mult = c(0, 0.28))
  ) +
  scale_fill_manual(values = PAL_BLOCK) +
  facet_wrap(~metric, scales = "free_x", nrow = 1) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "Count", y = NULL, fill = NULL,
       title = "Build funnel: records at each pipeline stage") +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 6.5),
    legend.key.size  = unit(3, "mm"),
    axis.text.y      = element_text(size = 6),
    strip.text       = element_text(size = 8, face = "bold")
  )

p_e <- ggplot(species_shift, aes(stage, n, fill = species)) +
  geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = scales::comma(n)),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 2.6, colour = "white") +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = pal_species_mv) +
  labs(x = NULL, y = "Records", fill = "Species",
       title = "Mouse / rat across stages") +
  theme(legend.position  = "bottom",
        axis.text.x      = element_text(size = 7))

# =============================================================================
# Row 3: reconciliation gates 467 → 320
# =============================================================================
p_f <- ggplot(route_sum, aes(n, fct_rev(route_pretty), fill = kept_by_gate)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), hjust = -0.2, size = 2.6) +
  scale_x_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.22))) +
  scale_fill_manual(values = c(`TRUE` = COL_ANALYSIS, `FALSE` = "grey70"),
                    labels  = c(`TRUE` = "kept (parenteral)", `FALSE` = "dropped")) +
  labs(x = "Records (n = 467)", y = NULL, fill = NULL,
       title = "Route gate",
       subtitle = "467 → 383; keeps iv / sc / ip / im") +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(3, "mm"))

p_g <- ggplot(tier_sum, aes(n, fct_rev(structure_tier), fill = kept_by_gate)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), hjust = -0.2, size = 2.6) +
  scale_x_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.25))) +
  scale_fill_manual(values = c(`TRUE` = COL_ANALYSIS, `FALSE` = "grey70"),
                    labels  = c(`TRUE` = "kept (complete)", `FALSE` = "dropped")) +
  labs(x = "Records (n = 383)", y = NULL, fill = NULL,
       title = "Structure gate",
       subtitle = "383 → 326; keeps complete SMILES") +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(3, "mm"))

p_h <- ggplot(assay_326, aes(half_life_h, fct_rev(assay_method), fill = kept_by_gate)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.3, width = 0.6, alpha = 0.9) +
  scale_x_hl(breaks = c(1 / 60, 5 / 60, 0.25, 1, 6, 24)) +
  scale_fill_manual(values = c(`TRUE` = COL_ANALYSIS, `FALSE` = "grey70"),
                    labels  = c(`TRUE` = "kept (named)", `FALSE` = "dropped")) +
  labs(x = "Half-life", y = NULL, fill = NULL,
       title = "Assay gate",
       subtitle = "326 → 320; drops 6 unknown-method") +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(3, "mm"))

# =============================================================================
# Row 4: half-life distributions by class in the 320-record analysis set
# =============================================================================

# Shared helper: violin + overlaid boxplot for a class variable.
# fill = species so the within-class species split is visible.
hl_violin <- function(data, x_var, pal_class = NULL) {
  p <- ggplot(data, aes(.data[[x_var]], half_life_h, fill = species)) +
    geom_violin(
      aes(group = interaction(.data[[x_var]], species)),
      position = position_dodge(width = 0.8),
      scale = "width", linewidth = 0.2, alpha = 0.75
    ) +
    geom_boxplot(
      aes(group = interaction(.data[[x_var]], species)),
      position = position_dodge(width = 0.8),
      width = 0.18, outlier.shape = NA, fill = "white",
      linewidth = 0.25, colour = "grey30"
    ) +
    scale_y_hl(breaks = c(1 / 60, 5 / 60, 0.25, 1, 6, 24, 72)) +
    scale_fill_manual(values = pal_species_mv) +
    labs(x = NULL, y = "Half-life", fill = "Species") +
    theme(
      axis.text.x     = element_text(size = 7),
      legend.position = "bottom",
      legend.text     = element_text(size = 7),
      legend.key.size = unit(3, "mm")
    )
  p
}

# Count annotations: n per class (pooled across species)
make_n_df <- function(data, x_var) {
  data |>
    count(.data[[x_var]], name = "n") |>
    rename(x = 1) |>
    mutate(y = min(data$half_life_h, na.rm = TRUE) * 0.35)
}

p_i_base <- hl_violin(iv_train, "Engineering") +
  labs(title = "Engineering class (IV training split)")

engin_n <- iv_train |> count(Engineering) |>
  mutate(half_life_h = min(iv_train$half_life_h, na.rm = TRUE) * 0.3)
p_i <- p_i_base +
  geom_text(data = engin_n, aes(Engineering, half_life_h,
                                 label = paste0("n=", n), fill = NULL),
            size = 2.3, colour = "grey40", vjust = 1)

topo_n <- iv_train |> count(Topology) |>
  mutate(half_life_h = min(iv_train$half_life_h, na.rm = TRUE) * 0.3)
p_j <- hl_violin(iv_train, "Topology") +
  geom_text(data = topo_n, aes(Topology, half_life_h,
                                label = paste0("n=", n), fill = NULL),
            size = 2.3, colour = "grey40", vjust = 1) +
  labs(title = "Topology (IV training split)")

chir_n <- iv_train |> count(Chirality) |>
  mutate(half_life_h = min(iv_train$half_life_h, na.rm = TRUE) * 0.3)
p_k <- hl_violin(iv_train, "Chirality") +
  geom_text(data = chir_n, aes(Chirality, half_life_h,
                                label = paste0("n=", n), fill = NULL),
            size = 2.3, colour = "grey40", vjust = 1) +
  labs(title = "Chirality (IV training split)")


# =============================================================================
# Assembly
# =============================================================================
LAB_S1 <- list(size = 8, face = "bold")

row1 <- ggarrange(
  wrap_labs(p_a, 40, 50), wrap_labs(p_b, 40, 50), wrap_labs(p_c, 40, 50),
  ncol = 3, widths = c(1.3, 1.05, 1.0),
  labels = c("a", "b", "c"), font.label = LAB_S1
)

row2 <- ggarrange(
  wrap_labs(p_d, 50, 60), wrap_labs(p_e, 40, 50),
  ncol = 2, widths = c(2.1, 0.9),
  labels = c("d", "e"), font.label = LAB_S1
)

row3 <- ggarrange(
  wrap_labs(p_f, 40, 50), wrap_labs(p_g, 40, 50), wrap_labs(p_h, 40, 50),
  ncol = 3, widths = c(1.1, 1.1, 1.2),
  labels = c("f", "g", "h"), font.label = LAB_S1
)

row4 <- ggarrange(
  wrap_labs(p_i, 30, 40), wrap_labs(p_j, 30, 40), wrap_labs(p_k, 30, 40),
  ncol = 3,
  labels = c("i", "j", "k"), font.label = LAB_S1
)

supp_fig1 <- ggarrange(row1, row2, row3, row4,
                        nrow = 4, heights = c(0.9, 1.2, 0.9, 1.0))

save_fig(supp_fig1, "supp_fig1", width = 320, height = 440)
