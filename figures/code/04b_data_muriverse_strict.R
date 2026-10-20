# The STRICT mouse + rat half-life evidence base, for Figure 2.
#
# Figure 2 describes the final modelling dataset as delivered. This script loads the
# modelling metadata plus the funnel and exclusion objects that Figure 2's selection
# panels (2b, 2c, 2n, 2o, 2q) draw.
#
# 04_data_muriverse.R describes the pooled MOUSE cohort as PEPlife2 ships it, and it stays
# that way: those tables are a faithful record of the source and the manuscript cites them.
# This script describes something narrower and newer. MuriVerse/build_strict_peptide_set.py
# applies the strict peptide-scale and structure-verification gates to the eligible rodent set.
# Every kept peptide SMILES round-trips to the sequence, while unresolved modifications are
# retained and reported in structure_tier.
#
# Figure 3 no longer trains on all kept observations. build_halflife_training_set.py applies
# three cohort requirements on top of the strict gates - a parenteral dose route, a
# whole-molecule SMILES and a bioanalytical method a source names - which take the 467 kept
# records to 320 modelling rows. That is why this script loads BOTH files and why they are
# not interchangeable: `sv_hl` is the evidence base, `sv_meta` is what the model sees, and
# the structure panels in Figure 2 draw from the first because the second is constant by
# construction. The concentration-decay phase filter still removes none.
#
# Figure 2 had been drawn from the mouse pool, and three of its panels cross-referenced
# Figure 3 by number: a funnel ending "final (Fig. 3) = 513", a PK-phase panel titled
# "Fig. 3 input", and an engineering panel titled "Fig. 3's strongest axis" whose top four
# classes - fc fusion, albumin-based, polymer conjugate, depot/particle - are precisely the
# records the strict gates remove. None of those three claims survived the move to the
# strict tier, which is why this script exists.
#
# Two consequences of the strict rule shape what the figure can show:
#
#   * `setting` and `design` are CONSTANT here. Every kept record is in vivo PK, because the
#     strict set is built from rodent/eligible_in_vivo_half_life.csv, already filtered to it.
#     The old panels a and c split on those two columns and would now draw a single bar, so
#     the figure splits on species instead - which is the axis that actually varies, and the
#     one Figure 3 reports its metrics over.
#   * the exclusions are the interesting half. 692 records were removed by a named gate and
#     every one is recorded with the gate and the verbatim text that triggered it, so "what
#     the strict rule costs" is drawable rather than assertable. The 147 records removed by
#     the three cohort requirements are recorded separately, in
#     halflife_cohort_reconciliation.csv, which is what the funnel's last three bars read.
if (!exists("ROOT")) {
  source(file.path(system("git rev-parse --show-toplevel", intern = TRUE), "figures/code/00_setup.R"))
}

DIR_STRICT <- file.path(ROOT, "MuriVerse", "rodent", "strict")
DIR_RODENT <- file.path(ROOT, "MuriVerse", "rodent")

SPECIES_LEVELS <- c("Mouse", "Rat")

# The kept records. One row per surviving measurement, carrying the axes Figure 3 trains
# on plus the provenance columns the eligible set already had.
sv_hl <- read_csv(file.path(DIR_STRICT, "in_vivo_half_life.csv"),
                  show_col_types = FALSE, guess_max = 1e5) |>
  mutate(PMID    = as.character(PMID),
         species = factor(str_to_title(species), levels = SPECIES_LEVELS))

# Every gate-excluded record, with the gate that removed it. A de-duplicated record is NOT
# here - nothing was wrong with it, it was a second copy of a structure already present at
# the same species, route and label - so kept + gate-excluded + de-duplicated is the published
# record count, and the funnel below is where the de-duplicated ones appear.
sv_excl <- read_csv(file.path(DIR_STRICT, "exclusions.csv"),
                    show_col_types = FALSE, guess_max = 1e5) |>
  mutate(PMID    = as.character(PMID),
         species = factor(str_to_title(species), levels = SPECIES_LEVELS))

# The funnel as build_strict_peptide_set.py itself recorded it, read from the file rather
# than re-implemented here - re-deriving these thirteen stages in R would be a second
# implementation of the gates that could silently disagree with the one that built the set.
sv_funnel <- read_csv(file.path(DIR_STRICT, "filtering_counts.csv"), show_col_types = FALSE) |>
  mutate(stage = fct_inorder(stage))

# The modelling rows: the strict set after the concentration-decay phase filter and the
# feature-level dedup, with the HMMER cluster and train/val split Figure 3 uses.
sv_meta <- read_csv(file.path(DIR_STRICT, "halflife_meta_with_split.csv"),
                    show_col_types = FALSE, guess_max = 1e5) |>
  mutate(PMID    = as.character(PMID),
         species = factor(str_to_title(species), levels = SPECIES_LEVELS),
         split   = factor(split, levels = c("train", "val")))

# Step-by-step reconciliation of the kept records -> modelling rows, written by
# build_halflife_training_set.py. Same argument as the funnel: read, not re-derived.
sv_reconciliation <- read_csv(file.path(DIR_STRICT, "halflife_cohort_reconciliation.csv"),
                              show_col_types = FALSE)

# The eligible rodent set is the actual input to build_strict_peptide_set.py. It includes
# records missed by the older organism-based pool and is the before population for the
# selection contrast in Figure 2's funnel and before/after panels.
sv_eligible <- read_csv(file.path(DIR_RODENT, "eligible_in_vivo_half_life.csv"),
                    show_col_types = FALSE, guess_max = 1e5) |>
  mutate(PMID    = as.character(PMID),
         species = factor(str_to_title(species), levels = SPECIES_LEVELS))

# --- source-level retrieval and curation status -----------------------------------------
# Restricted to the sources the STRICT set draws on, and pooled over both species'
# retrieval ledgers. The stage definition is the one 04_data_muriverse.R uses for the mouse
# cohort, so the two panels remain comparable: a source counts as curated only when an
# extraction JSON exists for it.
#
# 12 of the 186 sources are PATENTS, not papers - `PMID` holds a patent identifier for 109
# records, with EP2277890A1 alone accounting for 62. They get their own stage rather than
# being filtered out: dropping them would make this panel's total disagree with panel a's
# source count, and "the full-text pipeline cannot reach these at all" is a fact about the
# evidence base, not a blemish to hide.
retrieval <- c("mouse", "rat") |>
  map(\(sp) read_csv(file.path(ROOT, "MuriVerse", sp, "fulltext_retrieval_status.csv"),
                     show_col_types = FALSE) |>
        transmute(PMID = as.character(PMID), pmcid = as.character(pmcid),
                  fulltext_retrieved)) |>
  list_rbind() |>
  distinct(PMID, .keep_all = TRUE)

extracted_pmids <- c("mouse", "rat") |>
  map(\(sp) list.files(file.path(ROOT, "MuriVerse", sp, "extractions"), pattern = "\\.json$")) |>
  unlist() |>
  str_remove("\\.json$") |>
  unique()

CURATION_LEVELS <- c("curated extraction", "full text retrieved",
                     "PMC record only", "no PMC record", "patent (no PMID)")

sv_pmids <- sv_hl |>
  count(PMID, name = "records") |>
  mutate(is_pmid = str_detect(PMID, "^\\d+$")) |>
  left_join(retrieval, by = "PMID") |>
  mutate(
    stage = case_when(!is_pmid                            ~ CURATION_LEVELS[5],
                      PMID %in% extracted_pmids           ~ CURATION_LEVELS[1],
                      coalesce(fulltext_retrieved, FALSE) ~ CURATION_LEVELS[2],
                      !is.na(pmcid)                       ~ CURATION_LEVELS[3],
                      TRUE                                ~ CURATION_LEVELS[4]) |>
      factor(levels = CURATION_LEVELS)
  )

stopifnot(nrow(sv_pmids) == n_distinct(sv_hl$PMID))

# --- one half-life per unique sequence ---------------------------------------------------
# Repeated records of the same peptide are a single label. Keyed on the sequence, which is
# the unit the rest of the MuriVerse figures count in; the SMILES is what the strict gates
# verify, not what the cohort is reported in.
#
# Grouped WITHIN species, so a peptide dosed in both animals contributes one point to each
# violin. The two groups therefore sum to more than the distinct-sequence count, and the
# figure says so rather than letting the two numbers look like a discrepancy.
sv_seq_hl <- sv_hl |>
  filter(!is.na(sequence), sequence != "", half_life_h > 0) |>
  group_by(species, sequence) |>
  summarise(half_life_h = median(half_life_h), records = n(), .groups = "drop")

# --- tables -----------------------------------------------------------------------------
write_csv(sv_funnel, file.path(DIR_TABS, "table32_strict_funnel.csv"))

write_csv(
  sv_excl |>
    count(gate_failed, species, name = "records") |>
    pivot_wider(names_from = species, values_from = records, values_fill = 0) |>
    mutate(records = Mouse + Rat) |>
    arrange(desc(records)),
  file.path(DIR_TABS, "table33_strict_exclusion_gates.csv")
)

write_csv(
  sv_hl |>
    group_by(species) |>
    summarise(records = n(),
              source_sequence_strings = n_distinct(sequence), sources = n_distinct(PMID),
              median_half_life_h = median(half_life_h),
              min_half_life_h = min(half_life_h), max_half_life_h = max(half_life_h),
              .groups = "drop") |>
    left_join(
      sv_meta |>
        group_by(species) |>
        summarise(sequences = n_distinct(sequence), .groups = "drop"),
      by = "species"
    ) |>
    select(species, records, sequences, source_sequence_strings, sources,
           median_half_life_h, min_half_life_h, max_half_life_h),
  file.path(DIR_TABS, "table34_strict_cohort.csv")
)

write_csv(
  sv_hl |>
    group_by(hl_engineering) |>
    summarise(records = n(), median_half_life_h = median(half_life_h), .groups = "drop") |>
    arrange(desc(records)),
  file.path(DIR_TABS, "table35_strict_engineering.csv")
)

write_csv(
  sv_hl |>
    group_by(route_class) |>
    summarise(records = n(), median_half_life_h = median(half_life_h), .groups = "drop") |>
    arrange(desc(records)),
  file.path(DIR_TABS, "table36_strict_route_class.csv")
)

write_csv(
  sv_hl |>
    count(half_life_provenance, name = "records") |>
    mutate(share = records / sum(records)) |>
    arrange(desc(records)),
  file.path(DIR_TABS, "table37_strict_provenance.csv")
)

# --- the final dataset, as Figure 2 describes it ----------------------------------------
# Figure 2 describes `sv_meta` - the rows the model is actually fitted on - rather than
# the kept records. Any kept -> modelled step is part of the arrival story and belongs to the
# supplementary figure with the rest of it.

# The four structural axes the annotation block one-hots, in one long frame so the figure can
# draw them as a single faceted panel rather than four.
AXIS_LABELS <- c(n_ter_mod = "N terminus", c_ter_mod = "C terminus",
                 topology = "topology", chirality = "chirality",
                 structure_tier = "structure")

sv_axes <- sv_meta |>
  select(split, all_of(names(AXIS_LABELS))) |>
  mutate(chirality = recode(chirality, l = "all-L", d = "all-D", mix = "mixed L/D")) |>
  pivot_longer(all_of(names(AXIS_LABELS)), names_to = "axis", values_to = "level") |>
  mutate(axis = factor(AXIS_LABELS[axis], levels = unname(AXIS_LABELS)))

# The strict partner gate should leave no attached half-life-extension partner. Keep this
# count as an explicit audit invariant in the generated table.
sv_partners <- sv_meta |>
  count(partner_class, name = "records") |>
  arrange(desc(records))

write_csv(sv_partners, file.path(DIR_TABS, "table44_strict_partner_classes.csv"))

# How concentrated the dataset is on its sources. 186 sources sounds independent; the top one
# is a patent contributing an eighth of the records, so the figure draws the distribution
# rather than the count.
sv_sources <- sv_meta |>
  count(PMID, name = "records") |>
  arrange(desc(records)) |>
  mutate(rank = row_number(), cumulative = cumsum(records) / sum(records),
         is_patent = !str_detect(PMID, "^\\d+$"))

# Homology-cluster sizes, which is what the cross-validation folds are built from.
sv_clusters <- sv_meta |>
  count(cluster, name = "records") |>
  arrange(desc(records)) |>
  mutate(rank = row_number())

# Species-resolved breakdowns for the source-concentration and homology-cluster audit
# panels in fig2.R (kept available for ad-hoc use).
sv_sources_sp <- sv_meta |>
  count(PMID, species, name = "records") |>
  left_join(select(sv_sources, PMID, rank), by = "PMID")

sv_clusters_sp <- sv_meta |>
  count(cluster, species, name = "records") |>
  left_join(select(sv_clusters, cluster, rank), by = "cluster")

write_csv(
  sv_meta |>
    group_by(species) |>
    summarise(records = n(), sequences = n_distinct(sequence),
              clusters = n_distinct(cluster), sources = n_distinct(PMID),
              median_half_life_h = median(half_life_h),
              min_half_life_h = min(half_life_h), max_half_life_h = max(half_life_h),
              median_length = median(seq_length), .groups = "drop"),
  file.path(DIR_TABS, "table42_strict_modelling_set.csv")
)

write_csv(
  sv_sources |> slice_head(n = 20) |> select(PMID, records, rank, cumulative, is_patent),
  file.path(DIR_TABS, "table43_strict_source_concentration.csv")
)
