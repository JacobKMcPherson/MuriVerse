# Figure 2 (8 panels): the final MuriVerse rodent half-life training set, how it was built,
# and what it looks like. Supplementary Figure 1 is self-contained (preprocessing pipeline
# from the three source databases down to the 320 modelling rows) and does not source this
# script; the panel objects defined below (`p_1a`-`p_4b`) feed Figure 2 only.
#
# Row 1: a selection funnel, b effect of filtering, c peptide sequence space (UMAP, 2x)
# Row 2: d analysis set overview, e half-life distribution, f peptide length, g route
# Row 3: h half-life by route
#
if (!exists("ROOT")) {
  source(file.path(system("git rev-parse --show-toplevel", intern = TRUE), "figures/code/00_setup.R"))
}
if (!exists("sv_meta")) source(file.path(DIR_CODE, "04b_data_muriverse_strict.R"))

# Every panel is drawn from the final analysis set (sv_meta), so one record means the same
# thing throughout; only the funnel (b) and the filtering contrast (c) look beyond it.
# Every record in the set has a SMILES that was decoded back to its sequence, so the tiers
# below describe how much of the dosed molecule that SMILES captures, not whether one exists.
# Tier colours use the shared non-species categorical palette: strongest signal (complete)
# in maroon, middle in yellow, weakest/default in grey.
COL_CAT_YELLOW <- "#E69F00"
COL_CAT_MAROON <- "#882255"
COL_CAT_PURPLE <- "#332288"
COL_CAT_GREY   <- "grey70"
TIER_LEVELS  <- c("complete", "partly resolved", "modification not drawn")
TIER_COLOURS <- c(complete = COL_CAT_MAROON,
                  `partly resolved` = COL_CAT_YELLOW,
                  `modification not drawn` = COL_CAT_GREY)
stopifnot(all(sv_meta$structure_tier %in% TIER_LEVELS),
          !anyNA(sv_meta$smiles), all(sv_meta$smiles != ""))
sv_meta <- sv_meta |> mutate(structure_tier = factor(structure_tier, levels = TIER_LEVELS))

# Color palettes.
# Mouse is royal blue (MV_MOUSE), rat is medium sea green (MV_RAT); the per-species shades
# are used in panel a's funnel (very light at the eligible end, deepening to the full tone at
# the analysis-set end) and in panel d (the eligible-vs-analysed boxplots are faceted by
# species, each in its own species shade). Every other panel uses the solid tone.
MOUSE_SOLID <- MV_MOUSE
MOUSE_LIGHT <- "#E6EEFA"
RAT_SOLID   <- MV_RAT
RAT_LIGHT   <- "#E6F4EA"
pal_species <- c(Mouse = MOUSE_SOLID, Rat = RAT_SOLID)
pal_species_light <- c(Mouse = MOUSE_LIGHT, Rat = RAT_LIGHT)

# Non-species categorical palette. Four tones is the ceiling - everything above that
# collapses the figure into a soup of hues - and the three used by three-level panels are the
# same three used by other panels (yellow = fig3 validation split, maroon = source-this-arm,
# grey = default/weakest). Purple is kept for the one four-level panel (residue composition).
pal_prov <- c(`PEPlife2 transcription (no curated full text)` = "grey55",
              `curated full text available (ambiguous match, not auto-applied)` = COL_EX_VIVO,
              `matches curated full text` = COL_MV_GREEN,
              `corrected to curated full text` = COL_ANALYSIS)

h_breaks <- 10^seq(-3, 2, 1)

# Half-life is read in clock time rather than powers of ten on panels c, e and f. The axis runs
# from 1 minute to 72 hours: anything shorter than 1 min is shown at the 1 min end and
# anything at or above 72 h at the 72 h end, so the end labels are inequalities. This is a
# display cap only - the cohort, the model and every other figure keep the unclipped values.
# (13 of the 320 records are under 1 min, 50 under 5 min and 121 under 15 min; none reach 72 h, the longest is 62 h.)
HL_TICKS_H   <- c(1 / 60, 5 / 60, 0.25, 0.5, 1, 2, 4, 6, 12, 24, 48, 72)
HL_TICK_LABS <- c("\u2264 1 m", "5 m", "15 m", "30 m", "1 h", "2 h", "4 h", "6 h", "12 h", "24 h",
                  "48 h", "\u2265 72 h")
cap_hl <- function(h) pmin(pmax(h, min(HL_TICKS_H)), max(HL_TICKS_H))

# Mouse vs rat half-life, two-sided Wilcoxon rank-sum on the uncapped values (the cap is a display
# device only). Rank-based because half-lives are heavy-tailed; needs at least two records per
# species, otherwise there is nothing to test and no label is drawn. Shown in light grey: it is
# context for the reader, not a headline result, and it is not corrected for the many levels.
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "p<0.001",
                                           ifelse(p < 0.01, "p<0.01", sprintf("p=%.2f", p))))
species_p <- function(h, sp) {
  h <- h[!is.na(sp)]; sp <- sp[!is.na(sp)]
  if (sum(sp == "Mouse") < 2 || sum(sp == "Rat") < 2) return(NA_real_)
  suppressWarnings(wilcox.test(h ~ droplevels(factor(sp)), exact = FALSE)$p.value)
}
P_GREY <- "grey68"

# Figure 2 is drawn at the 180 mm page width the rest of the pipeline targets, so a panel in
# a four-column row is about 45 mm across. Type is named once here instead of being repeated
# in seventeen per-panel theme() calls: it was that repetition that let individual panels
# drift to 5.3 pt on a canvas twice this size, which is under 3 pt once reduced to the page.
SZ_TITLE   <- 6.6
SZ_SUB     <- 5.4
SZ_AXTITLE <- 5.8
SZ_AXTEXT  <- 5.2
SZ_LEGEND  <- 5.0
SZ_STRIP   <- 5.4
SZ_NUM     <- 1.55   # geom_text count labels, in mm (~4.4 pt)

theme_panel <- function(...) {
  theme(plot.title = element_text(size = SZ_TITLE, face = "plain"),
        plot.subtitle = element_text(size = SZ_SUB),
        axis.title = element_text(size = SZ_AXTITLE),
        axis.text = element_text(size = SZ_AXTEXT),
        strip.text = element_text(size = SZ_STRIP),
        legend.text = element_text(size = SZ_LEGEND),
        legend.title = element_text(size = SZ_LEGEND),
        plot.margin = margin(7, 4, 2, 4),
        ...)
}

# Legends that belong to one panel only (residue composition, structure tier) stay with their
# panel; the mouse/rat key does not, because seven panels carried an identical copy of it.
# Panels suppress it with SPECIES_LEGEND, and the figure keys the two species once, by hand, in
# the funnel's empty bottom-right corner (panel a). That gives every other panel its legend
# strip back as plotting area.
SPECIES_LEGEND <- theme(legend.position = "none")

PANEL_LEGEND <- function(nrow = 1) list(
  guides(fill = guide_legend(nrow = nrow)),
  theme(legend.position = "bottom", legend.key.size = unit(5, "pt"),
        legend.text = element_text(size = SZ_LEGEND),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt")))

# ============================================================================
# ROW 1: Overview panels (analysis set, HL distribution, length, seq quality)
# ============================================================================

# Figure 2 panel b: Analysis set counts by unit
p_1a <- sv_meta |>
  group_by(species) |>
  summarise(records = n(), sequences = n_distinct(sequence),
            clusters = n_distinct(cluster), .groups = "drop") |>
  pivot_longer(-species, names_to = "unit", values_to = "n") |>
  mutate(unit = fct_recode(unit, `homology\nclusters` = "clusters"),
         unit = fct_relevel(unit, "records", "sequences", "homology\nclusters"),
         species = fct_reorder(species, n, .fun = sum, .desc = TRUE)) |>
  ggplot(aes(unit, n, fill = species)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(0.8), vjust = -0.3, size = 1.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  scale_fill_manual(values = pal_species) +
  labs(x = NULL, y = "Count", fill = NULL) +
  # Title and subtitle hidden from rendering; the text now lives in fig2_caption.md (panel b).
  # labs(title = "Analysis set overview",
  #      # Records partition by species. Input sequences and homology clusters can occur in
  #      # both species, so their species-level counts are not additive. State the number of
  #      # shared clusters explicitly because grouped bars otherwise read as additive.
  #      subtitle = local({
  #        shared <- sum(sapply(split(sv_meta$species, sv_meta$cluster),
  #                             \(s) n_distinct(s) > 1))
  #        sprintf("%d records, %d sequences, %d clusters (%d span both species)",
  #                nrow(sv_meta), n_distinct(sv_meta$sequence),
  #                n_distinct(sv_meta$cluster), shared)
  #      })) +
  theme_panel(legend.position = "none",
              axis.text.x = element_text(size = SZ_AXTEXT))

# Figure 2 panel c: Half-life distribution
sp_n <- count(sv_meta, species)
p_1b <- ggplot(sv_meta, aes(species, cap_hl(half_life_h), fill = species)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  stat_summary(fun = median, fun.min = \(x) quantile(x, 0.25), fun.max = \(x) quantile(x, 0.75),
               geom = "pointrange", size = 0.25, linewidth = 0.4, colour = "grey15") +
  scale_y_log10(breaks = HL_TICKS_H, labels = HL_TICK_LABS, expand = expansion(mult = 0.04)) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = 1.9, colour = P_GREY,
           label = fmt_p(species_p(sv_meta$half_life_h, as.character(sv_meta$species)))) +
  # The record counts sit in the axis labels: at the top of the panel the mouse violin reaches
  # the 72 h end and the counts ran into it.
  scale_x_discrete(labels = setNames(paste0(sp_n$species, "\nn = ", sp_n$n), sp_n$species)) +
  scale_fill_manual(values = pal_species, guide = "none") +
  labs(x = NULL, y = "Half-life") +
  # labs(title = "Half-life distribution",
  #      # The median is exactly 0.345 h, which %.2f and the caption rounded in opposite
  #      # directions (0.34 against 0.35). Three significant figures removes the tie.
  #      subtitle = sprintf("median %.3g h", median(sv_meta$half_life_h))) +
  theme_panel()

# Peptide length distribution (defined here but not currently used by Figure 2 or
# Supplementary Figure 1; available as an audit panel for ad-hoc use)
p_1c <- ggplot(sv_meta, aes(seq_length, fill = species)) +
  geom_histogram(bins = 22, position = "identity", colour = NA) +
  scale_fill_manual(values = pal_species, guide = "none") +
  facet_wrap(~species, ncol = 1, scales = "free_y") +
  labs(x = "Residues", y = "Records") +
  # labs(title = "Peptide length distribution",
  #      subtitle = sprintf("median %d residues", as.integer(median(sv_meta$seq_length)))) +
  theme_panel()

# Panel: Residue composition (natural L only, D residues, non-standard residues)
res_data <- sv_meta |>
  mutate(has_d  = chirality %in% c("d", "mix"),
         has_ns = n_nonstandard > 0,
         residues = case_when(has_d & has_ns ~ "D and non-standard",
                              has_d          ~ "contains D-residues",
                              has_ns         ~ "contains non-standard residues",
                              TRUE           ~ "natural L-residues only"),
         residues = factor(residues, levels = c("natural L-residues only", "contains D-residues",
                                                "contains non-standard residues",
                                                "D and non-standard"))) |>
  count(species, residues)
p_1d <- ggplot(res_data, aes(n, fct_rev(species), fill = residues)) +
  geom_col(width = .55, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = if_else(n >= 12, as.character(n), "")),
            position = position_stack(vjust = .5, reverse = TRUE), size = 1.6) +
  scale_x_continuous(expand = expansion(mult = c(0, .02))) +
  # Palette avoids green/teal (reserved for the rat species colour in species-split panels).
  scale_fill_manual(values = c(`natural L-residues only` = "grey80",
                               `contains D-residues` = "#DDCC77",       # yellow
                               `contains non-standard residues` = "#E69F00",   # orange
                               `D and non-standard` = "#882255"),       # maroon
                    labels = c(`natural L-residues only` = "natural L only",
                               `contains D-residues` = "D-residues",
                               `contains non-standard residues` = "non-standard",
                               `D and non-standard` = "D + non-standard")) +
  labs(x = "Records", y = NULL, fill = NULL) +
  # labs(title = "Residue composition",
  #      subtitle = sprintf("%d of %d records are all natural L-residues",
  #                         sum(res_data$n[res_data$residues == "natural L-residues only"]), nrow(sv_meta))) +
  guides(fill = guide_legend(nrow = 2, override.aes = list(linewidth = 0))) +
  theme_panel() + PANEL_LEGEND(nrow = 2)

# ============================================================================
# ROW 2: Structure encoding (features, completeness, terminals, by HL)
# ============================================================================

# Panel: the whole-molecule requirement, drawn over the KEPT records rather than the
# analysed ones. Every analysed record is now `complete` by construction - the cohort
# requires it - so counting tiers inside the analysis set would draw one bar and say
# nothing. What the tiers still describe is the step that produced that set: of the records
# that passed every other strict gate, this is how many carried a SMILES that is the whole
# dosed molecule and how many did not, and the latter are the records the requirement
# removed. Panel o breaks the removed ones down by reason.
tier_data <- sv_hl |> count(species, structure_tier)
p_2a <- ggplot(tier_data, aes(n, fct_rev(species), fill = structure_tier)) +
  geom_col(width = .55, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = if_else(n >= 8, as.character(n), "")),
            position = position_stack(vjust = .5, reverse = TRUE), size = 1.6) +
  scale_x_continuous(expand = expansion(mult = c(0, .02))) +
  scale_fill_manual(values = TIER_COLOURS, drop = FALSE,
                    labels = c(complete = "complete", `partly resolved` = "partly resolved",
                               `modification not drawn` = "mod. not drawn")) +
  labs(x = "Records", y = NULL, fill = NULL) +
  # labs(title = "Whole-molecule requirement",
  #      subtitle = sprintf("%d of %d kept records are the whole molecule; the %d that are not are not analysed",
  #                         sum(sv_hl$structure_tier == "complete"), nrow(sv_hl),
  #                         sum(sv_hl$structure_tier != "complete"))) +
  # Fills are structure tiers here, not species, so this panel keeps its own key: it is the
  # only place the tier colours are named, and panel q reuses them without a key of its own.
  theme_panel() + PANEL_LEGEND(nrow = 2)

# Panel: Why a structure is not complete (the `unresolved` note the SMILES builder recorded)
reason_of <- function(x) case_when(
  is.na(x) | x == ""                                   ~ "none",
  str_detect(x, "^substitution")                       ~ "substitution unconfirmed",
  str_detect(x, "^topology|; topology")                ~ "cyclization unspecified",
  str_detect(x, "^non-standard residue")               ~ "non-standard not drawn",
  str_detect(x, "^modification text not recognised")   ~ "mod. not recognised",
  str_detect(x, "^backbone surrogate")                 ~ "backbone surrogate",
  str_detect(x, "^side-chain conjugate|^source reports") ~ "annotation vs sequence",
  TRUE                                                 ~ "other")
reason_data <- sv_hl |>
  filter(structure_tier != "complete") |>
  mutate(reason = reason_of(unresolved)) |>
  count(reason, species) |>
  mutate(reason = fct_reorder(reason, n, .fun = sum))
reason_totals <- reason_data |> group_by(reason) |> summarise(total = sum(n), .groups = "drop")
p_2b <- ggplot(reason_data, aes(n, reason, fill = species)) +
  geom_col(width = 0.65, position = position_stack(reverse = TRUE)) +
  geom_text(data = reason_totals, aes(x = total, label = total, fill = NULL), hjust = -0.25, size = 1.5) +
  scale_fill_manual(values = pal_species, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.34))) +
  labs(x = "Records", y = NULL, fill = NULL) +
  # labs(title = "Why a structure is not whole",
  #      subtitle = sprintf("%d records removed by the requirement in n", sum(reason_data$n))) +
  theme_panel() + SPECIES_LEGEND

# Species-stacked bars (h, i, j): mouse and rat counts stacked within each bar, in the species
# colours used throughout Figures 2 and 3, with the bar total at the end.
stack_species_layers <- function(totals) list(
  geom_col(width = 0.65, position = position_stack(reverse = TRUE)),
  geom_text(data = totals, aes(x = total, label = total, fill = NULL), hjust = -0.25, size = 1.5),
  scale_fill_manual(values = pal_species, drop = FALSE),
  scale_x_continuous(expand = expansion(mult = c(0, 0.34))))

# Panel: Terminal modifications, stacked by species and split by terminus
term_mods <- sv_meta |>
  transmute(id, species = as.character(species), N = n_ter_mod, C = c_ter_mod) |>
  pivot_longer(c(N, C), names_to = "terminus", values_to = "modification") |>
  filter(!is.na(modification), !modification %in% c("free", "not reported")) |>
  count(terminus, modification, species) |>
  mutate(terminus = factor(recode(terminus, N = "N-terminus", C = "C-terminus"),
                           levels = c("N-terminus", "C-terminus")),
         species = factor(species, levels = names(pal_species)),
         modification = fct_reorder(modification, n, .fun = sum))
term_totals <- term_mods |> group_by(terminus, modification) |>
  summarise(total = sum(n), .groups = "drop")
p_2c <- ggplot(term_mods, aes(n, modification, fill = species)) +
  stack_species_layers(term_totals) +
  facet_wrap(~terminus, nrow = 1) +
  labs(x = "Annotated records", y = NULL, fill = NULL) +
  # labs(title = "Terminal modifications",
  #      # nrow(term_totals) is the number of terminus x modification combinations, not a record
  #      # count, so this read "9 with reported mods" for a panel whose bars total far more than
  #      # nine. Counted over records, which is what the bars and the rest of the figure show.
  #      subtitle = local({
  #        bad <- c("free", "not reported")
  #        n <- sum((!is.na(sv_meta$n_ter_mod) & !sv_meta$n_ter_mod %in% bad) |
  #                 (!is.na(sv_meta$c_ter_mod) & !sv_meta$c_ter_mod %in% bad))
  #        sprintf("%d of %d records carry a terminal modification", n, nrow(sv_meta))
  #      })) +
  SPECIES_LEGEND +
  theme_panel()

# ============================================================================

# Panel: where each record's dose route came from. The route is a model input and half the
# cohort's `Dose Route` column is blank, so those records were resolved against the primary
# sources instead (MuriVerse/halflife_condition_curation.py, with the ledger and the
# verbatim quotes in MuriVerse/data/curation/route_curation.csv). This panel is the audit of
# that: `primary_arm` is a source that names the route for this very arm - a single-route
# study, or a half-life matched to its own row in the source's PK table - `primary_source` is
# a source that names the routes it used, all parenteral, with the record attributed at
# source level, and `database route column` is a route PEPlife2 transcribed, which no source
# was re-read for. Drawn as a share of each quintile because the question it answers is
# whether the weaker provenance is concentrated at one end of the label range.
PROV_LEVELS  <- c("primary_arm", "primary_source", "database_route_column")
PROV_LABELS  <- c(primary_arm = "source, this arm", primary_source = "source, study level",
                  database_route_column = "database column")
# Palette avoids green/teal (reserved for the rat species colour).
PROV_COLOURS <- c(`source, this arm` = "#882255",       # maroon
                  `source, study level` = "#E69F00",    # orange
                  `database column` = "grey62")
stopifnot(all(sv_meta$route_evidence_level %in% PROV_LEVELS))
p_2d <- sv_meta |>
  mutate(quintile = ntile(half_life_h, 5),
         quintile = factor(quintile, levels = 5:1,
                           labels = c("longest", "4th", "3rd", "2nd", "shortest")),
         provenance = factor(PROV_LABELS[route_evidence_level],
                             levels = unname(PROV_LABELS))) |>
  count(quintile, provenance) |>
  ggplot(aes(n, quintile, fill = provenance)) +
  geom_col(width = .65, position = position_fill(reverse = TRUE)) +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, .02))) +
  scale_fill_manual(values = PROV_COLOURS, drop = FALSE) +
  labs(x = "Share of records", y = "Half-life quintile", fill = NULL) +
  # labs(title = "Where the route came from",
  #      subtitle = sprintf("%d of %d records' routes read back to a primary source",
  #                         sum(sv_meta$route_evidence_level != "database_route_column"),
  #                         nrow(sv_meta))) +
  theme_panel() + PANEL_LEGEND(nrow = 2)

# Panel: Topology and chirality, stacked by species
axes_species <- sv_meta |>
  select(species, topology, chirality) |>
  mutate(chirality = recode(chirality, l = "all-L", d = "all-D", mix = "mixed L/D")) |>
  pivot_longer(c(topology, chirality), names_to = "axis", values_to = "level") |>
  count(axis, level, species) |>
  mutate(axis = factor(axis, levels = c("topology", "chirality")),
         level_id = fct_reorder(paste(axis, level, sep = "\r"), n, .fun = sum))
axes_totals <- axes_species |> group_by(axis, level_id) |> summarise(total = sum(n), .groups = "drop")
p_3a <- ggplot(axes_species, aes(n, level_id, fill = species)) +
  stack_species_layers(axes_totals) +
  scale_y_discrete(labels = \(x) sub("^.*\r", "", x)) +
  facet_wrap(~axis, scales = "free_y", ncol = 1) +
  labs(x = "Records", y = NULL, fill = NULL) +
  # labs(title = "Topology and chirality",
  #      subtitle = sprintf("%d cyclic, %d with D-residues", sum(sv_meta$topology == "cyclic"),
  #                         sum(sv_meta$chirality %in% c("d", "mix")))) +
  SPECIES_LEGEND +
  theme_panel()

# Panel: cohort composition and half-life, over every categorical attribute the figure
# describes, as two facets sharing one row ordering. This was two panels - route on its own,
# and a separate eligible-vs-analysed contrast over engineering, topology and chirality - which
# meant the reader could not ask "how many records, and what half-life" of the same attribute
# in one place. Merged, every row answers both at once: the left facet is how much evidence the
# level carries and in which animal, the right facet is what that evidence says.
#
# One cohort throughout: the 320-row analysis set. The filtering contrast the engineering
# section used to carry needed the eligible and 467-row structure-verified tiers, which the
# counts facet cannot show at the same time, so the merged panel describes the set the model
# actually trains on. Four of the nine eligible engineering classes survive the strict gates,
# which is why that section is shorter here than in the funnel.
#
# Levels are ordered by record count inside each section, so the row carrying the most evidence
# sits at the top of its section. The two columns are a count and a half-life and share no scale,
# so they are two plots joined side by side rather than two facets of one: one plot can then
# carry a linear count axis and the other a true log half-life axis with its own clock-time
# ticks, which a single facet grid cannot. Everything is the 320-row analysis set; the
# pre-filtering (eligible) layer this panel once overlaid in pale tones has been removed.
ATTR_SECTIONS <- c("Engineering", "Route", "Topology", "Chirality")
ATTR_TIERS <- c("eligible", "analysed")


# Chirality is taken from the raw PEPlife2 `Chirality` column on the eligible side because the
# curated `chirality_from_sequence` exists only on the analysed side and a contrast needs both;
# the two agree on all 320 analysed rows. "N.A." and the one composite entry are dropped, as is
# topology "not reported" and any route outside the four the cohort keeps.
attr_prep <- function(d, tier) {
  d |>
    transmute(species = as.character(species), half_life_h,
              Route       = as.character(route_class),
              Engineering = as.character(hl_engineering),
              Topology    = as.character(topology),
              Chirality   = recode(str_to_lower(as.character(
                              if ("Chirality" %in% names(d)) Chirality else chirality)),
                                   d = "all-D", l = "all-L", mix = "mixed L/D",
                                   .default = NA_character_)) |>
    pivot_longer(all_of(ATTR_SECTIONS), names_to = "section", values_to = "level") |>
    filter(!is.na(level)) |>
    mutate(tier = tier)
}
attr_long <- attr_prep(sv_meta, "analysed") |>
  filter(!(section == "Route" & !level %in% c("iv", "sc", "ip", "im")),
         !(section == "Topology" & !level %in% c("linear", "cyclic")))

# One factor over every level in the figure. The levels are unique across sections, so a single
# ordering can carry all four; free_y then shows each section only its own rows.
attr_order <- attr_long |>
  count(section, level) |>
  mutate(section = factor(section, levels = ATTR_SECTIONS)) |>
  arrange(desc(section), n) |>
  pull(level)
attr_long <- attr_long |>
  mutate(level = factor(level, levels = attr_order),
         section = factor(section, levels = ATTR_SECTIONS),
         species = factor(species, levels = c("Mouse", "Rat")),
         tier = factor(tier, levels = ATTR_TIERS))

attr_counts <- attr_long |> count(section, level, species)

ATTR_DODGE <- position_dodge(width = 0.85)

# One attribute = one panel: a records column and a half-life column sharing the level labels.
# `pad` adds empty slots at the top so a panel with fewer levels keeps the same row height - and
# so the same bar and box thickness - as its neighbours, instead of stretching to fill its cell.
attr_panel <- function(section, pad = 0, show_x = TRUE) {
  n_d  <- filter(attr_counts, section == !!section)
  hl_d <- filter(attr_long,   section == !!section)
  lv   <- intersect(attr_order, as.character(unique(n_d$level)))
  lim  <- c(lv, strrep(" ", seq_len(pad)))

  # Records by species. Bars grow leftwards from the edge they share with the half-life column,
  # towards the level labels; hjust is screen-aligned, so the count is right-aligned to sit
  # beyond the bar's end. All four panels share one count range, so a bar is the same length for
  # the same count anywhere in the figure; the limit leaves room past the longest (291).
  p_n <- ggplot(n_d, aes(n, level, fill = species)) +
    geom_col(width = 0.8, position = ATTR_DODGE) +
    geom_text(aes(label = n), position = ATTR_DODGE, hjust = 1.2, size = 1.6) +
    scale_fill_manual(values = pal_species, drop = FALSE) +
    scale_x_reverse(breaks = c(0, 250), limits = c(420, 0),
                    expand = expansion(mult = c(0, 0.02))) +
    scale_y_discrete(limits = lim) +
    labs(x = if (show_x) "Records", y = NULL, fill = NULL) +
    theme_panel() +
    theme(axis.text.y = element_blank(), legend.position = "none",
          plot.margin = margin(7, 2, 6, 2))

  # Half-life by species on a true log axis, read in clock time. Rows are shared with the records
  # column, so the level labels live there and are dropped here. Labels run down and to the right
  # from their ticks at 45 degrees; all eleven ticks stay so the scale still reads as the one in
  # panel c, but only four carry a label.
  p_hl <- ggplot(hl_d, aes(cap_hl(half_life_h), level, fill = species)) +
    geom_boxplot(outlier.size = 0.2, linewidth = 0.2, position = ATTR_DODGE,
                 width = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = pal_species, drop = FALSE) +
    scale_y_discrete(limits = lim) +
    scale_x_log10(breaks = HL_TICKS_H,
                  labels = ifelse(HL_TICKS_H == 1 / 60, "\u2264 1 m", ifelse(HL_TICKS_H == 4, "4 h", ifelse(HL_TICKS_H == 72, "72 h", ""))),
                  expand = expansion(mult = 0.04)) +
    labs(x = if (show_x) "Half-life", y = NULL) +
    theme_panel() +
    # The clock-time labels run down and to the right of their ticks at 45 degrees. A negative
    # angle descends left-to-right, and hjust = 0 with vjust = 1 anchors each label's START at
    # its tick, so the text trails away below-right instead of back over the panel. The 14 pt
    # right margin above is what keeps the last label ("72h") on the page.
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),
          legend.position = "none", plot.margin = margin(7, 14, 2, 2))

  # Joined with ggarrange rather than patchwork: align = "h" lines the two plots up top and
  # bottom even though the half-life plot's turned axis labels make its bottom margin taller, so
  # every level stays on one row across both columns.
  # The level labels are their own fixed-width column rather than the records plot's y axis:
  # as an axis their width came out of the records plot, so a panel with long labels (engineering)
  # got a narrower plotting area than one with short ones (route).
  p_lab <- ggplot(tibble(level = lv), aes(1, level, label = recode(level, `unmodified / none stated` = "unmodified",
                                  `backbone modification only` = "backbone mod.",
                                  `radio/fluorescent label only` = "radio/fluor. label"))) +
    geom_text(hjust = 1, size = 1.5, lineheight = 0.9) +
    scale_y_discrete(limits = lim) +
    scale_x_continuous(limits = c(0, 1), breaks = NULL, expand = c(0, 0)) +
    labs(x = if (show_x) " ", y = NULL) +
    theme_panel() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          axis.title.x = element_text(colour = NA),
          plot.margin = margin(7, 0, 6, 2))
  # The section name sits above the plot, left-aligned, instead of running up the left margin.
  # The significance column is commented out: panels d-g now show records and half-life only.
  # It was one p-value per level (Wilcoxon, mouse vs rat), in a narrow column of its own to the
  # right of the boxes so it never collided with them or with the clock-time ticks, headed by an
  # across-level Kruskal-Wallis on half-life over the levels of this attribute, both species
  # pooled (rank-based, like the per-level tests; uncorrected), in the same light grey.
  # kw_p <- suppressWarnings(kruskal.test(half_life_h ~ droplevels(level), data = hl_d)$p.value)
  # pv <- hl_d |>
  #   group_by(level) |>
  #   summarise(p = species_p(half_life_h, as.character(species)), .groups = "drop") |>
  #   mutate(lab = fmt_p(p))
  # p_pv <- ggplot(pv, aes(0, level, label = lab)) +
  #   geom_text(hjust = 0, size = 1.4, colour = P_GREY) +
  #   scale_y_discrete(limits = lim) +
  #   scale_x_continuous(limits = c(0, 1), breaks = NULL, expand = c(0, 0)) +
  #   labs(x = if (show_x) " ", y = NULL) +
  #   # Across-groups test as a two-line header directly above the per-level values, in the margin
  #   # over the plotting area, so it cannot meet the centred title above the plot.
  #   coord_cartesian(clip = "off") +
  #   annotate("text", x = 0, y = Inf, hjust = 0, vjust = -0.35, size = 1.4, lineheight = 0.9,
  #            colour = P_GREY, label = paste0("across groups\n", fmt_p(kw_p))) +
  #   theme_panel() +
  #   theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
  #         axis.title.x = element_text(colour = NA), plot.margin = margin(7, 2, 6, 0))
  ttl <- ggpubr::text_grob(section, size = SZ_TITLE, x = 0.5, hjust = 0.5)
  ggpubr::annotate_figure(
    # ggarrange(p_lab, p_n, p_hl, p_pv, ncol = 4, widths = c(1.25, 1, 1.15, 0.8), align = "h"),
    ggarrange(p_lab, p_n, p_hl, ncol = 3, widths = c(1.25, 1, 1.15), align = "h"),
    top = ttl)
}

p_d_eng   <- attr_panel("Engineering", show_x = FALSE)
p_d_route <- attr_panel("Route", show_x = FALSE)
p_d_topo  <- attr_panel("Topology")
p_d_chir  <- attr_panel("Chirality")

# Panel: Half-life by assay method and species. The assay is a model input (study conditions),
# and it moves the median half-life by an order of magnitude, so the figure shows it.
assay_hl <- sv_meta |>
  mutate(assay_method = fct_reorder(assay_method, half_life_h, median))
assay_med <- assay_hl |> group_by(assay_method) |> summarise(m = median(half_life_h), .groups = "drop")
p_3d <- ggplot(assay_hl, aes(half_life_h, assay_method, fill = species)) +
  geom_boxplot(outlier.size = 0.25, linewidth = 0.25, position = position_dodge(0.75), width = 0.65) +
  scale_x_hl(breaks = c(1 / 60, 5 / 60, 0.25, 1, 6, 24, 72)) +
  scale_fill_manual(values = pal_species, drop = FALSE) +
  labs(x = "Half-life", y = NULL, fill = NULL) +
  # labs(title = "Half-life by assay method",
  #      subtitle = sprintf("medians span %.2f to %.2f h", min(assay_med$m), max(assay_med$m))) +
  theme_panel() + SPECIES_LEGEND

# ============================================================================
# ROW 4: Cohort assembly and data provenance
# ============================================================================

# Source concentration (data origin distribution) - audit panel, available for ad-hoc use
# Stacked bars (mouse = MV_MOUSE, rat = MV_RAT) so species composition is visible alongside
# the concentration pattern.
p_4a <- ggplot(sv_sources_sp, aes(rank, records, fill = species)) +
  geom_bar(stat = "identity", width = 1, colour = NA) +
  geom_vline(xintercept = 10.5, linetype = "dotted", colour = "grey35", linewidth = 0.3) +
  annotate("text", x = 12, y = max(sv_sources$records) * 0.72,
           label = sprintf("top 10\n= %.0f%%", 100 * sv_sources$cumulative[10]),
           hjust = 0, size = 1.8, colour = "grey20", lineheight = 0.9) +
  scale_fill_manual(values = c(Mouse = MV_MOUSE, Rat = MV_RAT)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.09))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Sources (largest first)", y = "Records") +
  theme_panel() + SPECIES_LEGEND

# Homology clusters (cross-validation fold unit) - audit panel, available for ad-hoc use
# Stacked bars (mouse = MV_MOUSE, rat = MV_RAT). Log y-axis: bars start at 0.5 so
# singletons still have visible height. Species segments are computed explicitly so the
# 0.5 floor is applied to the lower segment rather than the total.
sv_clusters_sp_wide <- sv_clusters_sp |>
  tidyr::pivot_wider(names_from = species, values_from = records, values_fill = 0) |>
  arrange(rank) |>
  mutate(
    rat_ymin   = 0.5,
    rat_ymax   = 0.5 + Rat,
    mouse_ymin = rat_ymax,
    mouse_ymax = rat_ymax + Mouse
  )
p_4b <- ggplot(sv_clusters_sp_wide) +
  geom_rect(aes(xmin = rank - 0.5, xmax = rank + 0.5, ymin = rat_ymin, ymax = rat_ymax),
            fill = MV_RAT) +
  geom_rect(aes(xmin = rank - 0.5, xmax = rank + 0.5, ymin = mouse_ymin, ymax = mouse_ymax),
            data = filter(sv_clusters_sp_wide, Mouse > 0),
            fill = MV_MOUSE) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("largest %d\n%d singletons of %d",
                           sv_clusters$records[1],
                           sum(sv_clusters$records == 1), nrow(sv_clusters)),
           hjust = 1.05, vjust = 1.3, size = 1.8, colour = "grey20", lineheight = 0.95) +
  scale_y_log10(breaks = c(1, 10, 100), limits = c(0.5, max(sv_clusters$records) * 2),
                expand = c(0, 0)) +
  annotation_logticks(sides = "l", size = 0.2, short = unit(0.8, "mm"),
                      mid = unit(1.2, "mm"), long = unit(1.8, "mm")) +
  labs(x = "Cluster (largest first)", y = "Records (log scale)") +
  theme_panel()

# Figure 2 panel a: Selection funnel (eligible records → the strict analysis set)
STAGE_PASSED <- c(
  "eligible rodent records"             = "Eligible records",
  "not curated out"                     = "Curated review pass",
  "physiological circulating half-life" = "Normal animal and matrix",
  "no engineered partner class"         = "No PEG/Fc/albumin/polymer",
  "single peptide chain"                = "Single peptide chain",
  "faithful residue sequence"           = "Residue-level sequence",
  "peptide scale (<= 50 aa)"            = "≤50 residues",
  "no other attached partner"           = "No other partner",
  "elimination half-life"               = "Elimination half-life",
  "structure verified"                  = "Verified SMILES",
  "deduplicated within publication"     = "Dedup within source",
  "deduplicated across publications"    = "Dedup across sources",
  # The three cohort requirements, which the strict builder does not apply - they come from
  # build_halflife_training_set.py and are read from the reconciliation below.
  "parenteral dose route"               = "Parenteral route",
  "whole-molecule structure"            = "Whole molecule",
  "named bioanalytical method"          = "Named assay")

relabel <- function(x, map) {
  out <- unname(map[as.character(x)])
  ifelse(is.na(out), as.character(x), out)
}

# Two files, because the selection happens in two scripts. `sv_funnel` is the strict
# builder's own funnel, which ends at the kept records; `sv_reconciliation` is the cohort
# requirements build_halflife_training_set.py applies on top of those, and together they are
# the whole path from the eligible set to the analysis set. Leaving the second half out would
# show 467 kept records turning into 320 analysed ones with no step accounting for the gap.
#
# A stage that removes nothing is bookkeeping rather than selection - the builder records a
# row for retaining observations and another for curating the dose route, neither of which is
# a filter - so only the first row and the rows that actually remove records are drawn.
#
# Per-species counts at each stage are reconstructed by subtracting per-species exclusions
# from the eligible counts, so the stacked bars split each stage into mouse and rat.
excl_wide <- sv_excl |>
  count(gate_failed, species, name = "n") |>
  pivot_wider(names_from = species, values_from = n, values_fill = 0) |>
  rename(stage = gate_failed)
if (!"Mouse" %in% names(excl_wide)) excl_wide$Mouse <- 0L
if (!"Rat"   %in% names(excl_wide)) excl_wide$Rat   <- 0L

stage_order <- sv_funnel |>
  mutate(stage = as.character(stage)) |>
  filter(!str_starts(stage, "observations retained"),
         !str_starts(stage, "dose route curated")) |>
  pull(stage)

sv_hl_lookup <- sv_hl |>
  mutate(ID = as.character(ID),
         species = as.character(species))

reconc_species <- sv_reconciliation |>
  filter(removed > 0) |>
  transmute(stage = as.character(step),
            removed_ids = as.character(removed_ids)) |>
  rowwise() |>
  mutate(ids = list(str_split(removed_ids, "\\s+")[[1]]),
         Mouse = sum(sv_hl_lookup$ID %in% ids & sv_hl_lookup$species == "Mouse"),
         Rat   = sum(sv_hl_lookup$ID %in% ids & sv_hl_lookup$species == "Rat")) |>
  ungroup() |>
  select(stage, Mouse, Rat)

mouse0 <- sum(sv_eligible$species == "Mouse")
rat0   <- sum(sv_eligible$species == "Rat")

walk_counts <- function(stage_names, exclusions, start_mouse, start_rat) {
  m <- start_mouse; r <- start_rat
  out <- tibble(stage = character(), Mouse = integer(), Rat = integer(),
                m_removed = integer(), r_removed = integer())
  for (st in stage_names) {
    hit <- exclusions |> filter(stage == st)
    dm <- if (nrow(hit) > 0) hit$Mouse[1] else 0L
    dr <- if (nrow(hit) > 0) hit$Rat[1]   else 0L
    m <- m - dm; r <- r - dr
    out <- bind_rows(out, tibble(stage = st, Mouse = m, Rat = r,
                                  m_removed = dm, r_removed = dr))
  }
  list(frame = out, mouse = m, rat = r)
}

strict_steps <- walk_counts(stage_order[-1], excl_wide, mouse0, rat0)
reconc_steps <- walk_counts(reconc_species$stage, reconc_species,
                            strict_steps$mouse, strict_steps$rat)

funnel_rows <- bind_rows(
  tibble(stage = stage_order[1], Mouse = mouse0, Rat = rat0,
         m_removed = 0L, r_removed = 0L),
  strict_steps$frame,
  reconc_steps$frame,
  tibble(stage = "Analysis set", Mouse = reconc_steps$mouse, Rat = reconc_steps$rat,
         m_removed = 0L, r_removed = 0L)) |>
  mutate(records = Mouse + Rat,
         removed = m_removed + r_removed) |>
  # Drop the strict-gate stages that remove no records (keeping the eligible row and the
  # analysis-set row).
  filter(stage == stage_order[1] | stage == "Analysis set" | removed > 0) |>
  mutate(stage = fct_inorder(relabel(stage, STAGE_PASSED)))

# Per-species gradient: light at the eligible end, deepening to the full tone at the
# analysis-set end, one shade per stage.
funnel_mouse_cols <- setNames(colorRampPalette(c(MOUSE_LIGHT, MOUSE_SOLID))(nrow(funnel_rows)),
                              levels(funnel_rows$stage))
funnel_rat_cols   <- setNames(colorRampPalette(c(RAT_LIGHT,   RAT_SOLID  ))(nrow(funnel_rows)),
                              levels(funnel_rows$stage))

funnel_long <- funnel_rows |>
  select(stage, Mouse, Rat, records, removed) |>
  pivot_longer(c(Mouse, Rat), names_to = "species", values_to = "n") |>
  mutate(stage_species = interaction(species, stage, drop = TRUE))

stage_species_cols <- c(
  setNames(funnel_mouse_cols[as.character(levels(funnel_rows$stage))],
           paste("Mouse", levels(funnel_rows$stage), sep = ".")),
  setNames(funnel_rat_cols[as.character(levels(funnel_rows$stage))],
           paste("Rat",   levels(funnel_rows$stage), sep = ".")))

bar_totals <- funnel_rows |>
  mutate(label = if_else(removed > 0,
                         paste0(records, " (−", removed, ")"),
                         as.character(records)))

# The three named sets of the dataset chain, placed to the RIGHT of the
# bar total so the set name reads as a tag on the number ("1159 Eligible set", "467 Strict
# peptide set", "320 Analysis set"). Matched on the stage label rather than on a row number,
# so a change to the funnel drops the mark instead of mislabelling a different bar.
SET_MARKS <- c("Eligible records" = "Eligible set",
               "Verified SMILES"  = "Strict peptide set",
               "Analysis set"     = "Analysis set")
# The number-and-set-name string is pre-composed so one geom_text draws both and the set
# name's x position is anchored to the right edge of the number rather than guessed.
bar_totals <- bar_totals |>
  mutate(set_name = unname(SET_MARKS[as.character(stage)]),
         number_pad = strrep(" ", nchar(label)),
         set_tag    = if_else(is.na(set_name), "",
                              paste0(number_pad, "   ", set_name)))
stopifnot(sum(!is.na(bar_totals$set_name)) == length(SET_MARKS))

p_4d <- ggplot(funnel_long,
               aes(n, fct_rev(stage), fill = stage_species)) +
  geom_col(width = 0.65, position = position_stack(reverse = TRUE)) +
  geom_text(data = bar_totals, aes(x = records, y = fct_rev(stage), label = label),
            hjust = -0.08, size = 1.5, inherit.aes = FALSE) +
  # Set name to the right of the number, same baseline, bold so it reads as a tag.
  geom_text(data = filter(bar_totals, !is.na(set_name)),
            aes(x = records, y = fct_rev(stage), label = set_tag),
            hjust = -0.08, size = 1.5, fontface = "bold", colour = "grey15",
            inherit.aes = FALSE) +
  # Default discrete padding is 0.6 either side.
  scale_y_discrete(expand = expansion(add = c(0.6, 0.6))) +
  scale_x_continuous(breaks = c(0, 1000), expand = expansion(mult = c(0, 0))) +
  # The view is pinned rather than left to the data: the key below is an annotation, and an
  # annotation lengthens the scale, so without a fixed range moving it right just widened the
  # axis and it never reached the edge.
  coord_cartesian(xlim = c(0, 1900)) +
  scale_fill_manual(values = stage_species_cols, guide = "none") +
  # The one mouse/rat key for the whole figure. It is drawn by hand in the funnel's bottom-right
  # corner, which the bars never reach (the lowest stages end well left of x = 1300),
  # so it costs the panel no space. Bars are shaded by stage, so the key uses the solid tones
  # every other panel draws the two species in.
  annotate("point", x = 1050, y = c(2, 1), shape = 15, size = 2.2,
           colour = c(MOUSE_SOLID, RAT_SOLID)) +
  annotate("text", x = 1170, y = c(2, 1), label = c("Mouse", "Rat"), hjust = 0,
           size = 1.9) +
  labs(x = "Records", y = NULL) +
  theme_panel()

# Figure 2 panel h: Peptide sequence space (the peptide analogue of chemical space)
# Each unique sequence is placed by a two-dimensional UMAP of its mean-pooled
# ESM-2 embedding, the same sequence representation the model reads. The embeddings are
# memoised in MuriVerse/esm2_embedding_cache.npz (a zip of .npy files), which is read here
# directly so the figure needs no Python.
read_npy <- function(con) {
  readBin(con, "raw", 8)
  hdr <- rawToChar(readBin(con, "raw", readBin(con, "integer", 1, size = 2, endian = "little")))
  list(hdr = hdr, n = as.integer(str_match(hdr, "\\((\\d+),")[, 2]))
}
read_esm_cache <- function(path) {
  td <- tempfile(); dir.create(td); on.exit(unlink(td, recursive = TRUE))
  unzip(path, exdir = td)
  con <- file(file.path(td, "sequences.npy"), "rb")
  h <- read_npy(con)
  w <- as.integer(str_match(h$hdr, "<U(\\d+)")[, 2])
  codes <- matrix(readBin(con, "integer", h$n * w, size = 4, endian = "little"),
                  nrow = h$n, ncol = w, byrow = TRUE)
  close(con)
  seqs <- apply(codes, 1, \(r) intToUtf8(r[r > 0]))
  con <- file(file.path(td, "vectors.npy"), "rb")
  h <- read_npy(con)
  stopifnot(str_detect(h$hdr, "<f4"), !str_detect(h$hdr, "'fortran_order': True"))
  d <- as.integer(str_match(h$hdr, ",\\s*(\\d+)\\)")[, 2])
  vec <- matrix(readBin(con, "numeric", h$n * d, size = 4, endian = "little"),
                nrow = h$n, ncol = d, byrow = TRUE)
  close(con)
  rownames(vec) <- seqs
  vec
}

seq_space <- sv_meta |>
  group_by(sequence) |>
  summarise(records = n(), log_hl = median(log10(cap_hl(half_life_h))), cluster = first(cluster),
            # species of the sequence; a sequence seen in both animals is tagged "Both" so
            # the panel still colours it, but it does not appear as a species-specific point
            species = {
              u <- unique(as.character(species))
              if (length(u) > 1) "Both" else u
            },
            .groups = "drop") |>
  mutate(species = factor(species, levels = c("Mouse", "Rat", "Both")))
esm_vec <- read_esm_cache(file.path(ROOT, "MuriVerse", "esm2_embedding_cache.npz"))
stopifnot(all(seq_space$sequence %in% rownames(esm_vec)))
# UMAP on cosine distance between embeddings. Seeded and single-threaded so the layout is
# reproducible (uwot is only deterministic with n_threads = 1).
set.seed(1986)
seq_umap <- uwot::umap(esm_vec[seq_space$sequence, ], n_neighbors = 15, min_dist = 0.3,
                       metric = "cosine", n_threads = 1)
seq_space <- bind_cols(seq_space, tibble(UMAP1 = seq_umap[, 1], UMAP2 = seq_umap[, 2]))

big_cluster <- sv_clusters$cluster[1]
seq_space <- seq_space |>
  mutate(in_big = cluster == big_cluster) |>
  arrange(in_big)   # draw the largest-cluster points on top

# The third axis is half-life: a landscape of kernel-smoothed median log10 half-life over the
# UMAP1/UMAP2 plane. Each grid node is the Gaussian-weighted mean of the sequences around it, and
# nodes with too little sequence density are left out so the surface is not extrapolated into
# empty space. The surface is drawn as projected quads in ggplot (back to front), which keeps
# the panel in the same device and theme as the rest of the figure.
LAND_BW   <- 0.07 * mean(c(diff(range(seq_space$UMAP1)), diff(range(seq_space$UMAP2))))   # Gaussian bandwidth, UMAP units
LAND_MASS <- 0.6    # minimum summed kernel weight (in sequences) for a node to be drawn
LAND_N    <- 45
AZ <- 38 * pi / 180   # view azimuth
EL <- 30 * pi / 180   # view elevation
Z_SCALE <- 1.15       # height of the z axis relative to the floor width; raised from 0.75 so the
                      # landscape fills its cell instead of leaving a band of blank above and below it

xr <- range(seq_space$UMAP1) + c(-1, 1) * 0.04 * diff(range(seq_space$UMAP1))
yr <- range(seq_space$UMAP2) + c(-1, 1) * 0.04 * diff(range(seq_space$UMAP2))
gx <- seq(xr[1], xr[2], length.out = LAND_N)
gy <- seq(yr[1], yr[2], length.out = LAND_N)
nodes <- expand_grid(i = seq_along(gx), j = seq_along(gy)) |>
  mutate(x = gx[i], y = gy[j])
w <- exp(-(outer(nodes$x, seq_space$UMAP1, "-")^2 + outer(nodes$y, seq_space$UMAP2, "-")^2) /
           (2 * LAND_BW^2))
nodes <- nodes |>
  mutate(mass = rowSums(w), z = if_else(mass >= LAND_MASS, as.vector(w %*% seq_space$log_hl) / mass, NA_real_))

# The z axis runs from the 15 min cap up to the first clock-time tick that clears the surface;
# carrying it on to 72 h would leave most of the box empty, since no smoothed median gets near it.
z_hi <- log10(HL_TICKS_H[HL_TICKS_H >= 10^max(nodes$z, na.rm = TRUE)][1])
z_lim <- c(log10(min(HL_TICKS_H)), z_hi)
z_in <- log10(HL_TICKS_H) >= z_lim[1] - 1e-9 & log10(HL_TICKS_H) <= z_lim[2] + 1e-9
# Unit box: x and y span 0-1 (aspect of the PC plane is not preserved), z spans 0-Z_SCALE
box_xy <- function(x, y) list(x = (x - xr[1]) / diff(xr), y = (y - yr[1]) / diff(yr))
box_z  <- function(z) (z - z_lim[1]) / diff(z_lim) * Z_SCALE
project <- function(x, y, z) {
  bx <- box_xy(x, y); X <- bx$x - 0.5; Y <- bx$y - 0.5
  depth <- X * sin(AZ) + Y * cos(AZ)
  tibble(u = X * cos(AZ) - Y * sin(AZ), v = depth * sin(EL) + box_z(z) * cos(EL), depth = depth)
}

# Surface quads: drawn only where all four corners exist, painted far to near
node_z <- matrix(nodes$z, nrow = LAND_N, byrow = TRUE)   # [i, j]: nodes is i-major
quads <- expand_grid(i = 1:(LAND_N - 1), j = 1:(LAND_N - 1)) |>
  filter(!is.na(node_z[cbind(i, j)]), !is.na(node_z[cbind(i + 1, j)]),
         !is.na(node_z[cbind(i, j + 1)]), !is.na(node_z[cbind(i + 1, j + 1)])) |>
  mutate(id = row_number(), zmean = (node_z[cbind(i, j)] + node_z[cbind(i + 1, j)] +
                                      node_z[cbind(i, j + 1)] + node_z[cbind(i + 1, j + 1)]) / 4)
quad_pts <- bind_rows(
  quads |> mutate(k = 1, ci = i,     cj = j),
  quads |> mutate(k = 2, ci = i + 1, cj = j),
  quads |> mutate(k = 3, ci = i + 1, cj = j + 1),
  quads |> mutate(k = 4, ci = i,     cj = j + 1)) |>
  mutate(x = gx[ci], y = gy[cj], z = node_z[cbind(ci, cj)])
quad_pts <- bind_cols(quad_pts, select(project(quad_pts$x, quad_pts$y, quad_pts$z), u, v)) |>
  group_by(id) |> mutate(depth = mean(project(x, y, z)$depth)) |>
  ungroup() |>
  arrange(desc(depth), id, k)

# Floor rectangle, the vertical z axis on the left-most corner, and axis titles on the two
# front edges
floor_xy <- tibble(x = xr[c(1, 2, 2, 1)], y = yr[c(1, 1, 2, 2)])
floor <- bind_cols(floor_xy, project(floor_xy$x, floor_xy$y, rep(z_lim[1], 4)))
corner_z <- floor |> slice_min(u, n = 1)
z_axis_x <- corner_z$x; z_axis_y <- floor_xy$y[which.min(floor$u)]
# the 3D axis and its colour bar label only every few ticks, so the clock-time labels stay legible
z_lab <- z_in & HL_TICKS_H %in% c(1 / 60, 0.25, 1, 4)
z_ticks <- log10(HL_TICKS_H[z_lab])
z_axis <- bind_cols(tibble(z = z_ticks, lab = HL_TICK_LABS[z_lab]),
                    project(rep(z_axis_x, length(z_ticks)),
                            rep(z_axis_y, length(z_ticks)), z_ticks))
z_top <- project(rep(z_axis_x, 2), rep(z_axis_y, 2), z_lim)   # axis line spans the full z range
front <- which.max(-floor$depth)
edge_label <- function(from, to, text) {
  m <- (floor[c(from, to), c("u", "v")] |> summarise(across(everything(), mean)))
  d <- floor[to, c("u", "v")] - floor[from, c("u", "v")]
  ang <- atan2(d$v, d$u) * 180 / pi
  if (ang > 90) ang <- ang - 180
  if (ang < -90) ang <- ang + 180
  centre <- floor |> summarise(u = mean(u), v = mean(v))
  out <- c(m$u - centre$u, m$v - centre$v); out <- out / sqrt(sum(out^2)) * 0.07
  tibble(u = m$u + out[1], v = m$v + out[2], angle = ang, label = text)
}
nb <- c((front %% 4) + 1, ((front - 2) %% 4) + 1)   # the two corners adjacent to the front corner
edge_names <- sapply(nb, \(k) if (floor_xy$y[k] == floor_xy$y[front]) "x" else "y")
edge_labs <- bind_rows(lapply(seq_along(nb), \(m) edge_label(front, nb[m],
  if (edge_names[m] == "x") "UMAP1" else "UMAP2")))

# Sequences sit on the surface at their position, coloured by their own median half-life
seq_pts <- seq_space |>
  mutate(zs = vapply(seq_len(n()), \(r) {
    ww <- exp(-((UMAP1[r] - UMAP1)^2 + (UMAP2[r] - UMAP2)^2) / (2 * LAND_BW^2))
    sum(ww * log_hl) / sum(ww)
  }, numeric(1))) |>
  (\(d) bind_cols(d, project(d$UMAP1, d$UMAP2, d$zs)))() |>
  arrange(desc(depth))

# The smoothed landscape is drawn as a continuous cyan ramp keyed to the median log10
# half-life axis: pale where the surface is low (short half-life), deepening to saturated
# blue-green at the top of the z range. The surface is held at 50% alpha so the species
# points, which share the panel, read through it rather than over it. Points are coloured by
# species in the royal blue / green pair used across Figure 2; sequences dosed in both
# species sit in dark grey so the one or two bi-species points are visible without pretending
# they belong to either. Two fill scales (one for the surface, one for the points) are kept
# separate with ggnewscale.
suppressPackageStartupMessages(library(ggnewscale))
pal_species_umap <- c(Mouse = MOUSE_SOLID, Rat = RAT_SOLID, Both = "grey25")
p_6a <- ggplot() +
  geom_polygon(data = floor, aes(u, v), fill = "grey96", colour = "grey60", linewidth = 0.3) +
  geom_polygon(data = quad_pts, aes(u, v, group = id, fill = zmean), alpha = 0.85,
               colour = scales::alpha("white", 0.35), linewidth = 0.08) +
  # Grey scale from very light grey to very dark grey. Species points keep their own colours on
  # top, so a neutral surface keeps the hue free for them.
  scale_fill_gradientn(colours = c("#F5F5F5", "#BDBDBD", "#7A7A7A", "#3D3D3D", "#0D0D0D"),
                       limits = z_lim,
                       oob = scales::squish,
                      name = "median\nhalf-life",
                      breaks = log10(HL_TICKS_H[z_lab]), labels = HL_TICK_LABS[z_lab],
                      guide = guide_colourbar(
                        barwidth = unit(4, "pt"), barheight = unit(24, "pt"), order = 3)) +
  new_scale_fill() +
  # Every point gets the same thin white separation ring; outline colour carries no data.
  geom_point(data = dplyr::filter(seq_pts, !in_big),
             aes(u, v, fill = species, size = records),
             shape = 21, colour = "white", stroke = 0.25, alpha = 0.95) +
  geom_point(data = dplyr::filter(seq_pts, in_big),
             aes(u, v, fill = species, size = records),
             shape = 21, colour = "white", stroke = 0.25, alpha = 0.95) +
  annotate("segment", x = z_axis$u[1], xend = z_axis$u[1], y = z_top$v[1], yend = z_top$v[2],
           colour = "grey20", linewidth = 0.3) +
  geom_segment(data = z_axis, aes(x = u - 0.012, xend = u, y = v, yend = v), colour = "grey20",
               linewidth = 0.3) +
  geom_text(data = z_axis, aes(x = u - 0.02, y = v, label = lab), hjust = 1, size = 1.7) +
  annotate("text", x = min(z_axis$u) - 0.16, y = mean(z_top$v), angle = 90, size = 1.85,
           label = "Smoothed median half-life") +
  geom_text(data = edge_labs, aes(u, v, label = label, angle = angle), size = 1.85) +
  # Mouse and rat are keyed once, in panel a. "Both" is not a species and appears nowhere else,
  # so it stays: without it the dark grey points have no meaning.
  scale_fill_manual(values = pal_species_umap, name = NULL, breaks = "Both",
                    guide = guide_legend(order = 1)) +
  # The three keys form one column to the right of the plot, in the order "Both", records, then
  # the grey colourbar. The size key's symbols are black: left to inherit the points' white
  # separation ring they drew as empty white circles on a white page.
  scale_size_area(max_size = 4.5, breaks = c(1, 10, 50), name = "records",
                  guide = guide_legend(
                    order = 2,
                    override.aes = list(shape = 21, fill = "black", colour = "black"))) +
  coord_equal(clip = "off") +
  theme_void(base_size = 10) +
  # The keys are a column to the right of the plot rather than inside it: placed inside they
  # overlapped the right-hand peak. The plot is pushed left and the column is separated from it
  # by legend.box.spacing.
  theme(plot.title = element_text(size = SZ_TITLE, face = "plain", hjust = 0.5),
        plot.subtitle = element_text(size = SZ_SUB, hjust = 0.5),
        plot.margin = margin(7, 2, 2, 14),
        legend.position = "right", legend.justification = c(0, 1), legend.box = "vertical",
        legend.box.just = "left",
        legend.title = element_text(size = SZ_LEGEND),
        legend.text = element_text(size = SZ_LEGEND),
        legend.key.size = unit(6, "pt"), legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(10, "pt"), legend.spacing.y = unit(8, "pt"))

# ============================================================================
# ASSEMBLE MERGED FIGURE
# ============================================================================

# ggplot does not wrap plot titles, so on a 45 mm panel the longer ones ran past the edge and
# were cut mid-word ("Selection", "Incomp"). Wrapping is applied here, once, to every panel
# rather than by shortening seventeen individual labs() calls - the titles stay as written and
# the figure stays retunable from one place if the column count ever changes.
# Measured against the rendered panel: at SZ_TITLE on a four-column row about 20 title
# characters and 26 subtitle characters fit before the text runs under the next panel.
wrap_labs <- function(p, title_w = 20, sub_w = 26) {
  for (f in c("title", "subtitle")) {
    v <- p$labels[[f]]
    if (is.character(v) && length(v) == 1)
      p$labels[[f]] <- paste(strwrap(v, width = if (f == "title") title_w else sub_w),
                             collapse = "\n")
  }
  p
}

# Figure 2 (8 panels): the final MuriVerse rodent half-life training set, how it was arrived
# at, and what it looks like.
#   Row 1: a selection funnel, b analysis set overview, c half-life by species - the three
#          compact descriptions of the cohort, side by side
#   Rows 2-3, left half: d engineering, e route, f topology, g chirality - each a records
#          column and a half-life column; right half: h peptide sequence space (UMAP)
# Supplementary Figure 1 covers the complementary preprocessing story (3 source databases
# -> 467 structural -> 320 modelling rows) and is self-contained; it does not source this
# file. The panel objects defined above feed Figure 2 only.
LAB <- list(size = 8, face = "bold")
fig2_merged <- ggarrange(
  # Row 1: the funnel takes the widest share because its stage names are the longest labels in
  # the figure; b and c are small multiples of two or three bars and need little.
  ggarrange(wrap_labs(p_4d), wrap_labs(p_1a), wrap_labs(p_1b),
            ncol = 3, widths = c(1.5, 1, 0.85),
            labels = c("a", "b", "c"), font.label = LAB),
  # Rows 2-3: the four attribute panels as a 2 x 2 grid in the left half, the landscape in the
  # right half across both rows. Row 1 is 60 mm; rows 2 and 3 are 30 mm each for d-g, and h
  # spans both at 60 mm; the key of the landscape sits inside it (top right).
  ggarrange(
    ggarrange(wrap_labs(p_d_eng, 30, 40), wrap_labs(p_d_route, 30, 40),
              wrap_labs(p_d_topo, 30, 40), wrap_labs(p_d_chir, 30, 40),
              ncol = 2, nrow = 2, heights = c(1, 1),
              labels = c("d", "e", "f", "g"), font.label = LAB),
    wrap_labs(p_6a, 46, 58),
    ncol = 2, widths = c(1, 1), labels = c("", "h"), font.label = LAB),
  nrow = 2, heights = c(1, 1)
)

# 180 mm is save_fig's default and the single-page width the rest of the pipeline targets.
save_fig(fig2_merged, "fig2", width = 180, height = 120)

# Comprehensive data characterization figure (17 panels, 4 per row).
