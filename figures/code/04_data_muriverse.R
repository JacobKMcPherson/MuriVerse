# Mouse half-life cohort that MuriVerse is built from, plus its curation status.
if (!exists("ROOT")) {
  source(file.path(system("git rev-parse --show-toplevel", intern = TRUE), "figures/code/00_setup.R"))
}

DIR_HL_DATA <- file.path(ROOT, "MuriVerse", "data", "peplife2")
DIR_MV      <- file.path(ROOT, "MuriVerse", "mouse")

# The cohort, its per-record assay matrix and its recovered study design all come from
# MuriVerse/build_study_dossier.py; read them rather than re-deriving them here.
DESIGN_LEVELS <- c("in vivo PK", "ex vivo (fluid stability)", "ex vivo (tissue homogenate)",
                   "in vitro (other)", "ambiguous / no flag")

peplife_raw <- read_csv(file.path(DIR_MV, "records_provenance.csv"),
                        show_col_types = FALSE, guess_max = 1e5)

mouse_rows <- peplife_raw |>
  transmute(
    pmid      = as.character(PMID),
    pmid_valid = str_detect(pmid, "^\\d+$"),
    # 10 records have no sequence at all; keep them as "" so nchar() filters drop them
    sequence  = coalesce(sequence, ""),
    length    = suppressWarnings(as.numeric(Length)),
    year      = YEAR,
    sample    = str_to_lower(coalesce(sample, "")),
    cell_line = !in_cohort,
    setting   = factor(setting, levels = c("in vivo", "ex vivo", "in vitro", "other")),
    matrix,
    design    = factor(design, levels = DESIGN_LEVELS),
    animal_dosed,
    half_life_h,
    linear_cyclic = str_to_lower(coalesce(`Linear/Cyclic`, "not reported"))
  )

# Cell-line rows are counted in the caption but excluded from the modelling cohort
mv_cohort <- filter(mouse_rows, !cell_line)


# Curation status for the mouse cohort; extractions are the curated subset
retrieval <- read_csv(file.path(DIR_MV, "fulltext_retrieval_status.csv"), show_col_types = FALSE) |>
  transmute(pmid = as.character(PMID), pmcid, fulltext_retrieved)

extracted_pmids <- list.files(file.path(DIR_MV, "extractions"), pattern = "\\.json$") |>
  str_remove("\\.json$")

mv_pmids <- mv_cohort |>
  filter(pmid_valid) |>
  count(pmid, name = "records") |>
  left_join(retrieval, by = "pmid") |>
  mutate(
    has_pmc     = !is.na(pmcid),
    has_body    = coalesce(fulltext_retrieved, FALSE),
    is_extracted = pmid %in% extracted_pmids,
    stage = case_when(is_extracted ~ "curated extraction",
                      has_body ~ "full text retrieved",
                      has_pmc ~ "PMC record only",
                      TRUE ~ "no PMC record") |>
      factor(levels = c("curated extraction", "full text retrieved", "PMC record only",
                        "no PMC record"))
  )

# Half-life values reported per agent in the curated extractions
read_extraction <- function(f) {
  j <- jsonlite::fromJSON(f, simplifyVector = FALSE)
  agents <- j$agents %||% list()
  if (!length(agents)) return(tibble())
  tibble(
    pmid     = as.character(j$pmid),
    agent    = map_chr(agents, \(a) a$name %||% NA_character_),
    role     = map_chr(agents, \(a) a$role %||% NA_character_),
    modality = map_chr(agents, \(a) a$modality %||% NA_character_),
    has_sequence  = map_chr(agents, \(a) a$sequence %||% "not reported") |>
      str_detect("^[A-Za-z]") & !str_detect(map_chr(agents, \(a) a$sequence %||% "not reported"), "^not reported"),
    has_half_life = map_chr(agents, \(a) a$half_life %||% "not reported") |>
      (\(x) !is.na(x) & !str_detect(str_to_lower(x), "not reported"))()
  )
}

mv_agents <- list.files(file.path(DIR_MV, "extractions"), pattern = "\\.json$", full.names = TRUE) |>
  map(read_extraction) |>
  list_rbind()

# Curated setting/route/matrix nuance (MuriVerse/subset_by_nuance.py); paper-level, so it
# only replaces the record-level PEPlife2 approximation above for the 48 curated PMIDs, not the
# full cohort. Pipe-joined tag columns need splitting into one row per PMID-tag before counting.
NOT_A_ROUTE  <- c("not_administered", "not_reported")
NOT_A_MATRIX <- c("not_applicable", "not_reported")
split_tags <- function(x) str_split(coalesce(x, ""), "\\|") |> map(discard, ~ .x == "")

nuance <- read_csv(file.path(DIR_MV, "nuance_subset.csv"),
                   col_types = cols(pmid = col_character(), .default = col_guess())) |>
  filter(extraction_status == "extracted")

mv_nuance_routes <- nuance |>
  filter(is_in_vivo) |>
  transmute(pmid, route = split_tags(in_vivo_routes)) |>
  unnest(route) |>
  filter(!route %in% NOT_A_ROUTE) |>
  distinct(pmid, route) |>
  count(route, name = "publications") |>
  arrange(desc(publications))

mv_nuance_matrices <- nuance |>
  filter(is_ex_vivo) |>
  transmute(pmid, matrix = split_tags(ex_vivo_matrices)) |>
  unnest(matrix) |>
  filter(!matrix %in% NOT_A_MATRIX) |>
  distinct(pmid, matrix) |>
  count(matrix, name = "publications") |>
  arrange(desc(publications))

write_csv(mv_nuance_routes, file.path(DIR_TABS, "table7_muriverse_curated_in_vivo_routes.csv"))
write_csv(mv_nuance_matrices, file.path(DIR_TABS, "table8_muriverse_curated_ex_vivo_matrices.csv"))

# One half-life per unique sequence: repeated records of the same peptide are a single label
mv_seq_hl <- mv_cohort |>
  filter(nchar(sequence) > 0, half_life_h > 0) |>
  group_by(setting, sequence) |>
  summarise(half_life_h = median(half_life_h), records = n(), .groups = "drop")

write_csv(
  mv_cohort |>
    group_by(setting, matrix) |>
    summarise(records = n(), unique_sequences = n_distinct(sequence[nchar(sequence) > 0]),
              pmids = n_distinct(pmid[pmid_valid]),
              records_with_half_life = sum(half_life_h > 0, na.rm = TRUE),
              median_half_life_h = median(half_life_h[half_life_h > 0], na.rm = TRUE),
              .groups = "drop") |>
    arrange(setting, desc(records)),
  file.path(DIR_TABS, "table3_muriverse_mouse_cohort.csv")
)
write_csv(
  mv_seq_hl |>
    group_by(setting) |>
    summarise(unique_sequences = n(), records = sum(records),
              median_half_life_h = median(half_life_h),
              q1_half_life_h = quantile(half_life_h, 0.25),
              q3_half_life_h = quantile(half_life_h, 0.75), .groups = "drop"),
  file.path(DIR_TABS, "table4_muriverse_half_life_per_sequence.csv")
)
write_csv(
  count(mv_pmids, stage, name = "pmids"),
  file.path(DIR_TABS, "table5_muriverse_fulltext_curation.csv")
)
write_csv(
  mv_agents |>
    group_by(role, modality) |>
    summarise(agents = n(), with_half_life = sum(has_half_life),
              with_sequence = sum(has_sequence), .groups = "drop") |>
    arrange(desc(agents)),
  file.path(DIR_TABS, "table6_muriverse_extracted_agents.csv")
)
write_csv(
  mv_cohort |>
    group_by(design) |>
    summarise(records = n(), unique_sequences = n_distinct(sequence[nchar(sequence) > 0]),
              pmids = n_distinct(pmid[pmid_valid]),
              records_with_half_life = sum(half_life_h > 0, na.rm = TRUE),
              median_half_life_h = median(half_life_h[half_life_h > 0], na.rm = TRUE),
              .groups = "drop"),
  file.path(DIR_TABS, "table9_muriverse_record_level_design.csv")
)

# The analysis-ready set MuriVerse/subset_in_vivo_half_life.py writes and
# MuriVerse/build_halflife_training_set.py (fig3) reads directly - one row per surviving
# measurement, already carrying the axes fig3 trains and evaluates on.
mv_hl <- read_csv(file.path(DIR_MV, "in_vivo_half_life.csv"), show_col_types = FALSE, guess_max = 1e5) |>
  mutate(PMID = as.character(PMID))

# Attrition from the raw cohort to that set, reproducing subset_in_vivo_half_life.py's own
# filter order; the final count is read from the file itself rather than re-implementing the
# absorption-phase regex and dedup key here.
mv_hl_funnel <- tibble(
  stage = c("cohort (in_cohort)", "live animal dosed", "half-life reported",
            "design = in vivo PK", "absorption dropped\n+ deduplicated"),
  records = c(
    nrow(mv_cohort),
    sum(mv_cohort$animal_dosed == "dosed"),
    sum(mv_cohort$animal_dosed == "dosed" & !is.na(mv_cohort$half_life_h)),
    sum(mv_cohort$animal_dosed == "dosed" & !is.na(mv_cohort$half_life_h) &
        mv_cohort$design == "in vivo PK"),
    nrow(mv_hl)
  )
) |>
  mutate(stage = fct_inorder(stage))

write_csv(mv_hl_funnel, file.path(DIR_TABS, "table15_muriverse_halflife_funnel.csv"))

write_csv(
  mv_hl |>
    count(pk_phase, name = "records") |>
    mutate(share = records / sum(records)) |>
    arrange(desc(records)),
  file.path(DIR_TABS, "table16_muriverse_halflife_pk_phase.csv")
)
write_csv(
  mv_hl |>
    group_by(route_class) |>
    summarise(records = n(), median_half_life_h = median(half_life_h), .groups = "drop") |>
    arrange(desc(records)),
  file.path(DIR_TABS, "table17_muriverse_halflife_route_class.csv")
)
write_csv(
  mv_hl |>
    group_by(hl_engineering) |>
    summarise(records = n(), median_half_life_h = median(half_life_h), .groups = "drop") |>
    arrange(desc(records)),
  file.path(DIR_TABS, "table18_muriverse_halflife_engineering.csv")
)
write_csv(
  mv_hl |>
    count(half_life_provenance, name = "records") |>
    mutate(share = records / sum(records)) |>
    arrange(desc(records)),
  file.path(DIR_TABS, "table19_muriverse_halflife_provenance.csv")
)
