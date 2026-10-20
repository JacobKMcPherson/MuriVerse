# Data for Supplementary Figure 1: the training data from the three source databases down
# through every preprocessing gate to the 320 modelling rows. Reads all the user-supplied
# half-life databases, the two filtering-count tables and the strict reconciliation, plus
# the eligible and strict record files so a gate can be drawn before and after it fires.
if (!exists("ROOT")) {
  source(file.path(system("git rev-parse --show-toplevel", intern = TRUE), "figures/code/00_setup.R"))
}

DIR_HL_DATA   <- file.path(ROOT, "MuriVerse", "data")
DIR_RODENT    <- file.path(ROOT, "MuriVerse", "rodent")
DIR_STRICT    <- file.path(DIR_RODENT, "strict")

# ---------- Source-database species classifier (row 1 panels) ----------
# First match wins.
species_patterns <- c(
  "guinea pig" = "guinea",
  "human"      = paste0("human|volunteer|patient|healthy (male|female|adult|subject)|",
                        "\\bsubjects?\\b|\\bmen\\b|\\bwomen\\b|infants?|children|\\bchild\\b|",
                        "\\badults?\\b|\\bhsa\\b|caco-?2|hela\\b|hek ?293"),
  "mouse"      = "\\bmouse\\b|\\bmice\\b|murine|balb|c57",
  "rat"        = "\\brats?\\b|wistar|sprague",
  "monkey"     = "monkey|cynomolgus|macaque|rhesus|marmoset|baboon",
  "rabbit"     = "rabbit",
  "dog"        = "\\bdogs?\\b|beagle|canine",
  "pig"        = "\\bpigs?\\b|porcine|swine|minipig",
  "bovine"     = "bovine|calf|cattle|\\bcows?\\b",
  "sheep"      = "sheep|ovine|lamb",
  "horse"      = "horse|equine",
  "chicken"    = "chicken|avian",
  "fish"       = "zebrafish|trout|salmon",
  "insect"     = "drosophila|insect|locust",
  "microbial"  = "e\\.? ?coli|bacteri|yeast|fungal",
  "cell-free / enzyme" = paste0("trypsin|chymotrypsin|pepsin|protease|proteinase|peptidase|",
                                "elastase|papain|thermolysin|buffer|\\bpbs\\b|saline|",
                                "simulated (gastric|intestinal)|\\bsgf\\b|\\bsif\\b|enzyme|",
                                "aminopeptidase|carboxypeptidase|dpp-?(iv|4)|nep\\b|medium|",
                                "\\bmilk\\b")
)

non_animal <- c("microbial", "cell-free / enzyme", "unrecorded", "unclassified")

classify_species <- function(text, default = "unclassified") {
  t <- str_to_lower(coalesce(text, ""))
  hit <- map(species_patterns, \(p) str_detect(t, p)) |> as_tibble() |> as.matrix()
  first <- apply(hit, 1, \(r) if (any(r)) names(species_patterns)[which(r)[1]] else NA_character_)
  case_when(
    t %in% c("", "nan", "n.a.", "not mentioned") ~ "unrecorded",
    !is.na(first) ~ first,
    TRUE ~ default
  )
}

peplife <- read_csv(file.path(DIR_HL_DATA, "peplife2", "peplife2_records.csv"),
                    show_col_types = FALSE, guess_max = 1e5) |>
  transmute(
    database = "PEPlife2",
    sample   = `Test Sample`,
    species  = classify_species(sample),
    setting  = str_to_lower(str_squish(`Vivo/Vitro`)),
    setting  = case_when(
      str_detect(setting, "vivo/vitro") ~ "other",
      str_detect(setting, "^in vivo")   ~ "in vivo",
      str_detect(setting, "^in vi(t)?ro") ~ "in vitro",
      str_detect(setting, "^ex vivo")   ~ "ex vivo",
      TRUE ~ "other"
    ),
    half_life_h = suppressWarnings(as.numeric(`half-life(in Sec)`)) / 3600
  )

# ThPDB2 and PepTherDia2 are third-party databases and are not redistributed here
# (see DATA_LICENSES.md). With the raw files in place (scripts/fetch_thirdparty_data.sh) they are
# read and classified directly, and the species counts the figure needs are rewritten to
# results/database_comparison/. Without them, those committed counts stand in for the raw rows;
# they come from the code below, so Supplementary Figure 1 is unchanged either way.
FILE_THPDB       <- file.path(DIR_HL_DATA, "thpdb2", "ThPDB2.csv")
FILE_PEPTHERDIA  <- file.path(DIR_HL_DATA, "peptherdia2", "PepTherDia2.csv")
FILE_TP_COUNTS   <- file.path(ROOT, "results", "database_comparison",
                              "thpdb2_peptherdia2_species_counts.csv")

if (file.exists(FILE_THPDB) && file.exists(FILE_PEPTHERDIA)) {
  # Approved-drug databases report human clinical PK unless the text names another species
  thpdb <- read_csv(FILE_THPDB, show_col_types = FALSE, guess_max = 1e4) |>
    filter(!is.na(`Half Life`)) |>
    distinct(`Therapeutic Name`, `Half Life`, .keep_all = TRUE) |>
    transmute(database = "ThPDB2",
              sample   = `Half Life`,
              species  = classify_species(sample, default = "human"),
              setting  = "in vivo",
              half_life_h = NA_real_)

  peptherdia <- read_csv(FILE_PEPTHERDIA, show_col_types = FALSE) |>
    filter(!is.na(terminalhalflife), terminalhalflife != "-") |>
    transmute(database = "PepTherDia2",
              sample   = terminalhalflife,
              species  = classify_species(sample, default = "human"),
              setting  = "in vivo",
              half_life_h = NA_real_)

  dir.create(dirname(FILE_TP_COUNTS), showWarnings = FALSE, recursive = TRUE)
  bind_rows(thpdb, peptherdia) |>
    count(database, species, setting, name = "records") |>
    write_csv(FILE_TP_COUNTS)
} else {
  tp <- read_csv(FILE_TP_COUNTS, show_col_types = FALSE)
  expand <- function(db) {
    tp |> filter(database == db) |> uncount(records) |>
      mutate(sample = NA_character_, half_life_h = NA_real_) |>
      select(database, sample, species, setting, half_life_h)
  }
  thpdb      <- expand("ThPDB2")
  peptherdia <- expand("PepTherDia2")
}

hl_sources <- bind_rows(peplife, thpdb, peptherdia) |>
  mutate(database = fct_relevel(database, "PEPlife2", "ThPDB2", "PepTherDia2"))

# Species with < 1% of records are pooled so the stacked bars stay readable
species_levels <- hl_sources |>
  count(species, sort = TRUE) |>
  mutate(keep = species %in% non_animal | n / sum(n) >= 0.01) |>
  filter(keep, !species %in% non_animal) |>
  pull(species)

hl_sources <- hl_sources |>
  mutate(species_grp = case_when(
    species %in% c(species_levels, non_animal) ~ species,
    TRUE ~ "other animal"
  ) |>
    factor(levels = c(species_levels, "other animal", non_animal)) |>
    fct_drop())

# ---------- Build pipelines (row 2 panels) ----------
# Stage 1: PEPlife2 (4,500) narrowed to the eligible rodent in-vivo peptide set (1,159).
# Stage 2: that set narrowed to the structurally verified set (467), then reconciled to the
# 320 modelling rows.
elig_counts <- read_csv(file.path(DIR_RODENT, "eligible_filtering_counts.csv"),
                        show_col_types = FALSE) |>
  mutate(block = "PEPlife2 → eligible rodent set")

strict_counts <- read_csv(file.path(DIR_STRICT, "filtering_counts.csv"),
                          show_col_types = FALSE) |>
  mutate(block = "eligible rodent set → structurally verified")

recon <- read_csv(file.path(DIR_STRICT, "halflife_cohort_reconciliation.csv"),
                  show_col_types = FALSE)

# The two filtering tables plus the reconciliation gates, laid out as one top-to-bottom
# pipeline for the funnel panel. Each row is a stage; `records` is what remains, `removed`
# is what the gate drops. The eligible-table's final row and the strict-table's first row
# both read 1,159 — the strict table's version is dropped to keep the single pipeline.
build_pipeline <- bind_rows(
  elig_counts,
  strict_counts |> filter(stage != "eligible rodent records")
) |>
  mutate(stage_order = row_number(),
         stage = factor(stage, levels = stage))

# Reconciliation gates (467 → 320) kept as a separate block for row 3.
recon_gates <- recon |>
  filter(removed > 0) |>
  transmute(
    gate = c("route", "structure", "assay")[seq_len(n())],
    gate_label = c("Route gate", "Structure gate", "Assay gate")[seq_len(n())],
    rows_in, rows_out, removed, why,
    removed_ids = str_squish(removed_ids)
  )

# ---------- Snapshots at the 1,159, 467 and 320 marks (row 2f and row 3 panels) ----------
elig_records <- read_csv(file.path(DIR_RODENT, "eligible_in_vivo_half_life.csv"),
                         show_col_types = FALSE, guess_max = 1e5) |>
  select(ID, species, half_life_h, sequence)

kept_467 <- read_csv(file.path(DIR_STRICT, "in_vivo_half_life.csv"),
                     show_col_types = FALSE, guess_max = 1e5) |>
  select(ID, PMID, Assay, species, half_life_h, route_class, structure_tier, assay_method = Assay)

modelling_320 <- read_csv(file.path(DIR_STRICT, "halflife_meta_with_split.csv"),
                          show_col_types = FALSE, guess_max = 1e5) |>
  select(ID = peplife_id, species, assay_method)

# Species mix shift across the pipeline (panel f): three snapshots that are known exactly
# (1,159 eligible, 467 structural, 320 modelling). PEPlife2 raw rodent counts are not shown
# because its sample text does not reliably attach a route/sample pair that could be de-
# duplicated into a comparable "pre-eligibility" rodent record count.
species_shift <- bind_rows(
  elig_records |> count(species, name = "n") |> mutate(stage = "eligible\n(1,159)"),
  kept_467     |> count(species, name = "n") |> mutate(stage = "structurally\nverified (467)"),
  modelling_320 |> count(species, name = "n") |> mutate(stage = "modelling\nset (320)")
) |>
  mutate(stage = factor(stage, levels = c("eligible\n(1,159)",
                                          "structurally\nverified (467)",
                                          "modelling\nset (320)")),
         species = str_to_title(species))

# ---------- Reconciliation-gate snapshots (row 3) ----------
# Route gate (g): the 467's route_class distribution, with each class marked by whether the
# gate keeps it. Parenteral = iv/sc/ip/im; everything else drops. The pipeline writes an
# `is_parenteral` column; recompute it here from route_class so the panel does not depend
# on a column that another step could rename.
PARENTERAL <- c("iv", "sc", "ip", "im")

route_tbl <- kept_467 |>
  mutate(route_class = coalesce(route_class, "not reported"),
         route_pretty = case_when(
           route_class == "iv" ~ "iv",
           route_class == "sc" ~ "sc",
           route_class == "ip" ~ "ip",
           route_class == "im" ~ "im",
           route_class == "oral" ~ "oral",
           route_class == "intranasal" ~ "intranasal",
           route_class == "rectal" ~ "rectal",
           route_class == "not reported" ~ "not reported",
           TRUE ~ "unresolved"
         ),
         kept_by_gate = route_class %in% PARENTERAL)

route_order <- route_tbl |>
  count(route_pretty, kept_by_gate, name = "n") |>
  arrange(!kept_by_gate, desc(n)) |>
  mutate(route_pretty = factor(route_pretty, levels = route_pretty))

route_sum <- route_tbl |>
  count(route_pretty, kept_by_gate, name = "n") |>
  mutate(route_pretty = factor(route_pretty, levels = levels(route_order$route_pretty)))

# Structure gate (h): same shape, for structure_tier on the 383 kept by the route gate.
TIER_LEVELS <- c("complete", "partly resolved", "modification not drawn")

tier_tbl <- kept_467 |>
  filter(route_class %in% PARENTERAL) |>
  mutate(structure_tier = factor(structure_tier, levels = TIER_LEVELS),
         kept_by_gate = structure_tier == "complete")

tier_sum <- tier_tbl |>
  count(structure_tier, kept_by_gate, name = "n")

# Assay gate (i): the 326 that passed route+structure, by assay_method label. The six that
# the gate drops are the ones halflife_meta_with_split tagged as `unknown`. The strict
# meta table carries the resolved assay_method for the 320 kept; for the six removed, pull
# the label and half-life from the full strict in_vivo table joined on ID.
removed_ids_assay <- recon_gates |>
  filter(gate == "assay") |>
  pull(removed_ids) |>
  str_split(" ") |>
  unlist() |>
  as.integer()

# halflife_meta_with_split carries assay_method for the 320 kept; use `in_vivo_half_life`
# for the six removed records' half-life. The assay-method text itself for those six is
# the record's own `Assay` free-text column, which is what the gate flagged as unknown.
assay_326 <- kept_467 |>
  filter(route_class %in% PARENTERAL, structure_tier == "complete") |>
  left_join(modelling_320 |> select(ID, assay_method_resolved = assay_method), by = "ID") |>
  mutate(
    assay_method = coalesce(assay_method_resolved, "unknown"),
    kept_by_gate = assay_method != "unknown"
  )

# Panel i ordering: put "unknown" last so it reads as the dropped tail. Keep medians for
# a secondary read.
assay_order <- assay_326 |>
  group_by(assay_method) |>
  summarise(median_hl = median(half_life_h, na.rm = TRUE),
            n = n(),
            kept_by_gate = all(kept_by_gate),
            .groups = "drop") |>
  arrange(kept_by_gate, desc(median_hl)) |>
  mutate(assay_method = factor(assay_method, levels = assay_method))

assay_326 <- assay_326 |>
  mutate(assay_method = factor(assay_method, levels = levels(assay_order$assay_method)))

# ---------- Export compact summary tables ----------
write_csv(hl_sources |> count(database, species_grp, setting, name = "records"),
          file.path(DIR_TABS, "supp_table1_source_database_species.csv"))

write_csv(build_pipeline |> select(block, stage, records, removed, rationale),
          file.path(DIR_TABS, "supp_table2_preprocessing_pipeline.csv"))

write_csv(route_sum |> rename(route = route_pretty) |> arrange(desc(kept_by_gate), desc(n)),
          file.path(DIR_TABS, "supp_table3_route_gate.csv"))

write_csv(tier_sum, file.path(DIR_TABS, "supp_table4_structure_gate.csv"))

write_csv(assay_326 |> count(assay_method, kept_by_gate, name = "n"),
          file.path(DIR_TABS, "supp_table5_assay_gate.csv"))
