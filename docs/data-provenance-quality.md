# MuriVerse data provenance, quality, and origin analysis

## Purpose
MuriVerse combines peptide half-life evidence from heterogeneous sources. This document defines the minimum provenance contract and the analyses required before pooled modeling.

## Source classes
- **PEPlife 2.0**: database/API-derived peptide half-life records. Preserve PEPlife ID, PMID, source row, original half-life text/unit, test sample, assay, route/dose when present, and the API/export artifact used.
- **PepTherDia 2.0**: curated therapeutic-peptide database. Preserve native record identifier, source file/version, cited publication identifier, species/matrix, and the exact source field from which a half-life value is taken.
- **ThPDB 2.0**: therapeutic peptide/protein database. Preserve native record identifier, source file/version, cited publication identifier, species/matrix, and the exact source field from which a half-life value is taken.
- **Primary literature / PubMed-linked full text**: preserve PMID/PMCID/DOI where available, extraction JSON path, quoted evidence, article section, and whether the value was directly reported or transformed.
- **Derived/model-ready datasets**: every row must remain traceable to one or more upstream records; preprocessing must never erase provenance.

## Required row-level provenance fields
Recommended canonical fields:
`record_id, source_database, source_record_id, source_file, source_file_sha256, source_version, pmid, pmcid, doi, reference, species, strain, sex, age, matrix, assay_context, study_type, peptide_name, sequence, sequence_modified, half_life_reported, half_life_unit_reported, half_life_seconds, half_life_transform, dose, route, bioanalytical_method, extraction_method, evidence_locator, evidence_quote, cohort_rule, inclusion_reason, exclusion_reason, duplicate_group_id, quality_flags`.

Do not merge records solely by peptide name or sequence. Distinct experiments, matrices, species, routes, publications, and database records are distinct observations.

## Data-quality dimensions
1. **Identity**: sequence completeness, canonical/non-canonical residues, terminal modifications, conjugates, cyclic status.
2. **Outcome integrity**: reported half-life vs derived midpoint/range parsing; unit consistency; censored values (`>`, `<`); implausible or nonnumeric values.
3. **Biological context**: mouse/rat/human, plasma/serum/whole blood, in vivo vs ex vivo, strain/sex/age.
4. **Methodological context**: assay type, bioanalytical method, incubation conditions, dose/route for in-vivo PK.
5. **Provenance completeness**: source record ID, primary citation, exact field or full-text locator, extraction evidence.
6. **Duplication/leakage**: duplicate database records, the same primary experiment represented in multiple databases, same sequence in train and validation, or multiple rows from a single publication split across folds.
7. **Representativeness**: source-specific sequence lengths, chemical modifications, publication years, matrices, and half-life ranges.

## Origin-difference analysis
At minimum, compare each source database and primary-literature cohort using:
- record count and unique peptide/sequence count;
- unique PMID count and records per PMID;
- missingness by field;
- duplicate sequence and duplicate sequence+context frequency;
- half-life distribution on raw and log scales;
- peptide length distribution;
- matrix, assay/study type, route, and modification composition;
- overlap of sequences and PMIDs between sources;
- source-specific outlier rate;
- source label predictability from molecular features as a diagnostic for dataset shift;
- train/validation source composition and group leakage.

Observed differences must be interpreted as a mixture of biology, assay design, database scope, curation policy, and publication selection—not automatically as intrinsic peptide effects.

## Modeling implications
Use source-aware evaluation. Prefer grouped splits by sequence and, where feasible, publication/experiment. Report performance stratified by source, matrix, and assay context. Consider source/database as a covariate only when it will be known at deployment; otherwise use it for auditing, stratification, weighting, sensitivity analysis, or hierarchical/random-effects approaches rather than as a predictive shortcut.

## Current state (rodent rebuild)

> **Historical record.** The verification results in this section were run against a
> 1,317-record rodent table that concatenated the two per-species
> `<species>/in_vivo_half_life.csv` files. That concatenation has since been removed from
> the project: the pipeline now runs from the 1,159-record eligible set
> (`rodent/eligible_in_vivo_half_life.csv`), and the per-species tables it was built from
> are still present and unchanged. Counts below are preserved as the record of what was
> checked at the time and should not be read as current cohort sizes.

Status of this contract against the combined mouse + rat cohort actually built by
`MuriVerse/build_halflife_training_set.py --species mouse rat`. Full detail, including the
per-record corrections and the evidence for each, is in
[`../MuriVerse/rat/README.md`](../MuriVerse/rat/README.md).

**Satisfied.**

- *Biological context — species.* `species` is a column on every record and on both saved
  datasets, corrected per record where PEPlife2 named the wrong animal.
- *Provenance completeness — source class.* `source_type` distinguishes journal article,
  patent and unidentified; `source_is_review` marks secondary re-reports. 157 of 1,317 in-vivo
  records are patent-sourced and 44 come from reviews, neither of which was previously visible.
- *Outcome integrity — units.* `UNITS_CORRECTIONS` fixes values carried against the wrong unit
  (one found: minutes recorded as hours, a 60-fold error).
- *Outcome integrity — not-a-half-life.* `CURATED_NOT_A_HALF_LIFE` removes values that are not
  half-lives of the peptide: t-test p-values read out of a statistics column, and the clearance
  of a small-molecule payload carried by a peptide-decorated vehicle.
- *Methodological context — phase and route.* `pk_phase` splits absorption, distribution,
  terminal and non-concentration readouts, with `CURATED_ABSORPTION` and `PHASE_OVERRIDES` for
  records whose label is not in the value text.
- *Duplication/leakage.* HMMER homology clustering (`phmmer` all-vs-all, single linkage; it
  replaced MMseqs2 in this update) runs over the pooled set and whole clusters move together,
  so no sequence and no cluster spans train and validation in either branch (verified: 0). The
  split is balanced within each species. `CURATED_DUPLICATES` removes tables transcribed twice
  under two spellings of `Test Sample`.
- *Traceability from published artefact to training set.* Every filter between
  `in_vivo_half_life.csv` and the training set is now recorded — rows in, rows out, the reason,
  and the PEPlife2 IDs removed — in `<cohort>/halflife_cohort_reconciliation.csv`, with an
  assertion that the chain starts at the published record count and ends at the training row
  count. Before this, 513, 467 and 459 circulated as three different totals for the mouse
  cohort with nothing reconciling them.
- *Identity — non-canonical residues.* Hyphen-delimited non-natural residue codes (Aib, Nle,
  Sar, …) are collapsed to a single `X` instead of silently expanding into three standard
  residues; `sequence_nonstandard` records the ambiguity and non-standard letters that remain
  (204 of 1,212 training rows).

**Outstanding.**

- *Provenance completeness.* `source_file_sha256`, `source_version` and `duplicate_group_id`
  from the canonical field list are still not emitted.
- *Identity — stereochemistry.* **Resolved for the strict tier.** The lowercase residues that
  mark D-amino acids were never lost by the database — they are lost by the pipeline, which
  reads `sequence` (a flattened, upper-cased derivative) instead of `Sequence` (PEPlife2's
  original). The strict tier now reads the original, so a mixed-chirality peptide is drawn
  with its D positions exactly where the source puts them; the chirality gate excludes
  nothing, where it previously excluded 32 records for a `Mix` annotation the source in fact
  states precisely. The D series is also now a true enantiomer: the threonine and isoleucine
  side-chain centres were previously left in the L configuration whatever `Chirality` said,
  so a D peptide was emitted as a diastereomer rather than a mirror image. The round trip
  verifies configuration, not only constitution. The pooled tier still reads the flattened
  column.
- *Identity — the sequence often is not the dosed molecule.* Pervasive and the most consequential
  item here. Identical strings carry order-of-magnitude different labels because the
  half-life-extension chemistry (a fatty acid, a PEG, a 600-residue PAS or PA polymer, an IgG Fc)
  is not encoded. `hl_engineering` recovers the mechanism as a covariate. **Addressed for the
  strict tier** (see below), which keeps only records whose whole molecule is drawn and
  verified; it remains outstanding for the pooled tier, which is now retained for the
  provenance funnel rather than for modelling.
- *Modeling implications — source-aware evaluation.* `species`, `source_type` and
  `hl_engineering` are all now available for stratified reporting, but nothing in the current
  two-branch layout uses them, and `species` is not a feature of either representation.
- *Curation coverage.* 155 of 804 rat records (19%) are curated against primary full text,
  limited by retrievability: only 57 of 235 rat PMIDs are in the PMC OA subset and the PEPlife2
  PDF mirror is unresponsive. The mouse cohort is curated to a comparable but higher fraction.
- *Cross-source pooling.* PepTherDia 2.0 and ThPDB 2.0 are archived under `MuriVerse/data/` but
  are not yet merged into the cohort, so the overlap and dataset-shift analyses this document
  requires between source databases have not been run.

## The strict tier (peptide sequence + full structure)

`MuriVerse/build_strict_peptide_set.py` builds a second, stricter cohort for modelling, and
`build_halflife_training_set.py --tier strict` (the default) trains on it. Its contract is one
rule: **a record is kept only when the molecule that was dosed can be written down
completely** — the sequence column is a real sequence of template-bearing residues, and every
structural feature the source records is present in the emitted SMILES (topology, chirality,
both termini, and every clause of `Chem-MOD`).

1,317 pooled records → **420 kept** over 163 studies (329 rat, 91 mouse). The funnel, the per-record exclusion reasons and the
judgement calls are documented in
[`../MuriVerse/rodent/strict/README.md`](../MuriVerse/rodent/strict/README.md); the machine-
readable forms are `rodent/strict/filtering_counts.csv` and `rodent/strict/exclusions.csv`.

What this fixes against the contract above:

- *Identity — conjugates and terminal modifications.* Acetylation, amidation, methyl esters,
  side-chain lipidation with its γ-glutamyl and OEG spacers, and exactly-structured covalent
  labels are all in the string. Records whose partner cannot be drawn truthfully — Fc,
  albumin, polymer, depot and particle vehicles, PEG in any form, glycans, drug payloads,
  proprietary dyes — are excluded rather than represented by a string that is not the
  molecule. Consequently the strict tier needs no `termini` feature
  block: there is nothing attached that is not drawn.
- *Identity — cyclic status.* Head-to-tail macrocycles and positioned disulfides are in the
  string. Two under-specified forms are now resolved rather than abandoned, because in each the
  inferred ring is the only one the molecule can close: a bare `Cyclic` on a sequence with
  exactly two cysteines is the disulfide between them, and a bare `Cyclic` on a sequence with
  no cysteine at all is a backbone macrocycle (both records this applies to confirm it in their
  own sources — `c[wLwReQeR]` is written in head-to-tail bracket notation in PMID 28691076, and
  PMID 31866565 names its peptide a cyclopeptide in the title). A cyclization that remains
  under-determined still excludes the record: a bare `Bicyclic` needs two rings and four
  cysteines, and `Bicyclic(C9-C15,C15-C21)` puts one cysteine in two disulfides.
- *Outcome integrity / identity — verification.* Every kept record is checked twice: its bare
  peptide is decoded back to a sequence by an independent SMILES→sequence reader and must
  return the recorded sequence (up to rotation for a macrocycle), and the full molecule's exact
  mass minus the bare peptide's must equal the sum of the attached fragments to 0.01 Da.
- *Duplication/leakage.* De-duplication is keyed on the **structure** — species, SMILES, route,
  label — which the pooled tier could not do, because its string did not identify the molecule.

New quality findings this work surfaced, each now recorded rather than latent:

- *Two residue codes can be run together with no delimiter.* `SarLys-RPPGFSPf` is
  sarcosine-lysine followed by a heptapeptide. The tokenizer had no rule for an undelimited
  pair of codes, so it fell through to the one-letter path and expanded `SarLys` into
  Ser-Ala-Arg-Leu-Tyr-Ser — one dipeptide became six invented residues and a 10-mer became a
  14-mer, with spurious D assignments from the lower-case letters. The round-trip check could
  not catch it: the wrong sequence decodes consistently back to itself, so the string verified
  cleanly against the sequence it was built from. It is the SEQUENCE that was wrong, not the
  drawing of it. Both affected records are now refused with that reason rather than silently
  mis-expanded; splitting undelimited codes would need the tokenizer to guess where one ends.
- *A glossary entry can be mistaken for an unapplied modification.* `Chem-MOD` often carries a
  legend explaining the sequence's own notation — "O denotes L-ornithine", "pGlu =
  Pyroglutamate" — which adds no atoms, because the residue is already in the sequence and
  already templated. The blocking rule for untemplated residues matched the legend's text and
  demoted ten records to `modification not drawn` although their structures were complete. The
  escape that existed for this matched residue CODES only, and a legend spells the residue
  out; it also omitted sarcosine and pyroglutamate, which are templated through the backbone
  builder rather than the side-chain table. Both gaps are closed, and the escape still requires
  the sequence to actually carry the residue, so a legend for a residue that is absent, or a
  clause that conjugates rather than defines, keeps blocking.
- *A conjugate recipe can lose a component silently.* The spacer pattern required a word
  boundary before the spacer name, which `2xOEG` does not provide — the multiplier is glued to
  the token. The spacer was therefore never matched, and because an unmatched token was simply
  skipped, `C18 diacid-gGlu-2xOEG` assembled as C18-diacid-gGlu and reported SUCCESS: a
  conjugate two OEG units short of the molecule, presented as complete. This is the worst
  failure mode available to that function, because it does not refuse, it under-reports. No
  published structure was affected — none of the 13 drawn conjugates uses multiplier notation,
  and the records that do were blocked for other reasons — so the defect was latent. The
  multiplier is now separated before matching, and an audit refuses any recipe that names a
  spacer it could not place, rather than assembling what it understood.
- *A half-life can be a legal claim rather than a measurement.* WO2012138867A2 contributed a
  record reading 1.5 h for an intermedin/adrenomedullin analogue, annotated `Rat plasma
  (Intravenous)` with an ELISA. The patent contains no rat pharmacokinetic experiment at all:
  its only rodent dosing is intraperitoneal for a blood-pressure readout, and the sole "1.5
  hours" in the document is claim 1 — "A biologically active intermedin or adrenomedullin
  peptide having a serum half-life of greater than 1.5 hours" — a patentable threshold the
  invention must clear. No structural gate can see this; the peptide draws and verifies
  cleanly. Removed through `CURATED_EXCLUSIONS`, where the reason is recorded.
- *Methodological context — dose route, curated from the sources.* `Dose Route` is blank on
  half the cohort, which left `route_class` reading "not reported" for 241 modelling rows.
  Those are now resolved against the primary sources — full text where reachable, patent text,
  abstract, and NLM MeSH indexing — in `MuriVerse/data/curation/route_curation.csv`, one row
  per source with the verbatim quote and its origin. `route_evidence_level` keeps arm-level
  evidence (`primary_arm`), study-level evidence (`primary_source`) and a transcribed database
  value (`database_route_column`) distinguishable in the shipped metadata, and a curated route
  without its quote is refused at build time rather than shipped. The model cohort then
  requires a parenteral route, because an enteral or mucosal terminal slope is an absorption
  rate under flip-flop kinetics and not an elimination half-life.
- *Methodological context — bioanalytical method.* `Assay` names a separation rather than a
  detector on some records ("TLC"), the readout rather than the method ("Serum decay
  experiment", "Plasma stability assay"), or nothing. The vocabulary now covers the named
  methods it was missing (EIA, LC/MS, LCHRMS, AlphaLISA, western blot) and carries `bioassay`
  as its own level, for a plasma concentration read off a cell-based reporter or an inhibition
  titre rather than measured physically — a quantity with different error behaviour, since an
  inactive metabolite an immunoassay would count does not register at all. What the column
  cannot say is taken from the sources in `MuriVerse/data/curation/assay_curation.csv`.
  Records where nothing names a method are removed: what quantity decayed is not a detail of a
  half-life, it is the half-life. The six that reached that state had a median label of 23 h
  against the cohort's 0.427 h and held three of the four longest values in the set.
- *A free-text annotation can contradict the sequence column.* `Chem-MOD` substitution clauses
  were previously assumed to be already applied to the sequence. They are now checked against
  it, and 28 records are excluded: those the sequence refutes, and those it cannot confirm
  either way because the stated position already holds the product residue. Three liraglutide
  rows are the second kind — `Lys-30 is replaced by Arg` checks out, yet the stored sequence
  is still native GLP-1(7-37), giving `C172H265N41O51` against the drug's `C172H265N43O51`.
  Catching that needs the parent sequence, which the row does not carry, so the class is
  excluded rather than kept with a flag.
- *A missing sequence could be drawn as a peptide.* When `build_smiles` is handed pandas'
  `NaN` for an empty sequence cell, its string form upper-cases to `"NAN"` — three valid
  residue codes — and the builder drew the tripeptide Asn-Ala-Asn and reported it
  `complete`. This is a live defect on any caller that passes the raw column, which is how
  `build_rodent_smiles.py` and `build_strict_peptide_set.py` both call it. It did **not**
  corrupt the published artefacts: `rodent/smiles_all_records.csv` reports all 17
  sequence-less records as `ambiguous residue` with no SMILES, because on that path the
  value arrived as an empty string, which the existing `not sequence` guard caught. No
  training row was ever affected either, since the pooled builder filters sequence-less
  records out before drawing. Fixed at the source: `build_smiles` now refuses a blank or
  `NaN` sequence explicitly, and with that fix it reproduces the published `smiles_coverage`
  column exactly (1317/1317) where before it disagreed on precisely those 17 rows.
- *A label can be a stability readout rather than a half-life, with nothing structural to
  show it.* PEPlife2 ID 5555 (`Free-ADI`) is an 83-residue fragment of a ~410-residue enzyme
  whose 1445 h "Activity Half Life" was the largest label in the cohort by 8.6×, and it drew
  and verified perfectly. It is curated out by ID in `CURATED_EXCLUSIONS`; the ceiling drops
  to 168 h.
- *A carrier can be named only in `Name`.* `Au‐AR pep‐PROTAC` is a gold-nanoparticle conjugate
  at 26.3 h — roughly 100× the strict set's median — whose every annotation column reads `Free`
  and `unmodified`, and which `hl_engineering` also misses. `Name` is now screened; `Origin` is
  deliberately not, being provenance prose about the sequence's ancestry rather than anything
  attached.
- *Radiotracer labelling is a measurement artefact, not a structure.* 76 kept records carry
  one, flagged in `tracer_label`. Isotopic substitution (³H, ¹⁴C, ¹³C, ¹⁵N) changes mass, not
  constitution. Radioiodination is treated the same way on a dosed-material argument — tracer
  specific activities run far below one iodine per peptide — rather than a drawing-difficulty
  one. A radiometal named without its chelator is excluded, the chelator being a real covalent
  adduct of unknown identity.

`smiles_coverage` semantics were also contradictory across the working tree and are now
settled, in favour of the definition the docstrings and the published column already agreed
on: coverage describes the **peptide's own axes** — constitution, topology, configuration —
and an attached partner never demotes it. Partner presence travels on
`n_partner`/`c_partner`, the terminal-encoded flags and the `termini` feature block instead.
The orphan `peptide_core` value, which the `build_halflife_training_set` docstring records as
an abandoned design because it "measured what RDKit could draw rather than the molecule", is
gone; the vocabulary is `complete | partial | none` as published. Note this makes "complete"
a statement about the peptide, *not* about the whole dosed molecule — which is precisely the
reading the strict tier exists to provide, and why its completeness is a gated, verified
property rather than a flag.

`--tier pooled` runs again. `peptide_smiles.describe_termini` and its supporting
`PARTNER_MASS_KDA` table and PEG-size parser were recovered from commit `982b854` rather than
reconstructed, and the recovered code reproduces all five published columns
(`n_ter_class`, `c_ter_class`, `partner_class`, `peg_kda`, `partner_kda`) exactly, on all
1,317 rows.

Still outstanding for this tier:

- 198 records are excluded for an unaccounted `Chem-MOD` clause. A meaningful share are
  legitimately undrawable (non-standard residues the template set has no entry for,
  N-methylated backbones, `structure given in paper`), but ~13 are only blocked because a dye's
  attachment site is unstated (`FITC labeled`, `TAMRA-labeled`, `Biotinylation`), and the
  proprietary cyanines (Cy5, Cy5.5, Cy7, Alexa, DiR, ZW800) are blocked for want of a published
  structure rather than for a chemical reason. Both are recoverable by extending
  `peptide_mods.LABEL`.
- PEG is now excluded in every form, including discrete oligomers whose structure is exact.
  That is a consistency choice rather than a chemical limit, and it costs one record.
- 203 records are excluded by the peptide-scale cap, which is now **50 residues**. That boundary
  is a choice, not a measurement, and it was 100 until it was tightened: 100 was the point where
  ESM-2's 1022-token window starts truncating, which is a property of the encoder rather than of
  the biology. 50 is the working peptide/small-protein boundary, past which half-life is
  increasingly set by hydrodynamic size rather than by the pharmacophore. Tightening it costs
  the kept set 38 records (8%) and moves the gate from fourth-largest to second.

## The two sequence columns

The single most consequential finding of the cleanup, and it is a pipeline bug rather than a
database defect. PEPlife2 ships two sequence columns:

- **`Sequence`** — the original. Marks a D residue by lower-casing its letter
  (`rPKPfQwFwLL`, `SCSlPQT`) and delimits a non-natural residue with hyphens
  (`H-Aib-QGTFTSDYSKYL`, `pGlu-LYENKPRRPYIL`).
- **`sequence`** — a flattened derivative: upper-cased and punctuation-stripped.

**403 of 1,300 rodent rows differ**, and the difference is not cosmetic. Flattening destroys
the D/L pattern, and it turns `pGlu-LYENKPRRPYIL` into `PGLULYENKPRRPYIL` — which the prose
detector then correctly rejects as not being a sequence. So two of the "outstanding" items in
this document were artefacts of reading the wrong column: the stereochemistry gap, and a
share of the records written off as prose.

The strict tier now reads `Sequence`, with the flattened column retained as a cross-check —
where the two disagree about which residues a peptide contains, the record is refused rather
than one of them being preferred. Parsing is deliberately conservative for delimited
sequences, because case cannot be trusted to mean D there: `pIPhe` is one residue whose
letters are each individually a valid code, so expanding it would invent five residues from
one.

## Two sources assessed and then removed from the workflow

Both were built into the pipeline, measured, and taken back out. The measurements are kept
here because they are the reason not to try again, not because the code still exists.

**PEPlife2's structure archive** (`data/peplife2/peplife2_structures.zip`, 1,837 PDB files
keyed by PEPlife2 ID) was wired in as a second opinion on which residues a record contains.
Three properties, each measured rather than assumed, limited it to that and then ruled it out:

1. *The ID-to-structure mapping is unreliable.* Over records with both a sequence and a
   single-chain model, 223 agree and 171 contradict — some wholesale, e.g. ID 6231's column
   reads `PTLHTHNRRRRR` (12 residues) while its PDB holds a 29-residue exendin-like peptide,
   and IDs 6231–6234 share one structure. The archive's own `old_rename_log.csv` records 759
   files renumbered `1xxx → 4xxx`, but the residual offsets are scattered, so this is not a
   correctable systematic renumbering.
2. *The coordinates do not encode stereochemistry.* RDKit CIP assignment from 3D, on records
   the `Chirality` column marks `D`, returns S — L — at every centre: PEPstrMOD built
   everything in the natural series regardless of annotation. Chirality cannot be read off
   these files, and the disulfide pairings and conformations are predictions too.
3. *The archive is heterogeneous.* 520 PEPstrMOD models, 330 CASP-format predictions, 588
   bare coordinate lists, and 399 downloaded experimental PDB entries — 113 of them
   multi-chain complexes carrying receptors, antibody H/L pairs, waters, metals and glycans,
   i.e. structures of the parent protein rather than the dosed peptide.

As a corroboration gate it excluded 47 records for contradicting it. Since the source that
disagrees is wrong at least 43% of the time where it can be checked, those exclusions were
not defensible, and both the gate and the module are gone.

**PepTherDia 2.0 and ThPDB 2.0** were assessed as additional record sources and **neither can
contribute rodent half-life records**: PepTherDia contains zero rodent mentions in any column
(its `terminalhalflife` is human clinical pharmacokinetics), and ThPDB has 5 of 3,380
half-life entries mentioning a rodent, duplicates of roughly 2 unique values. Pooling either
into a rodent label would be a species error. PepTherDia was briefly used instead as an
external check on the SMILES builder, and that check is also removed; what it established
while it ran is recorded in the strict set's README rather than in live code.

So the cohort-size constraint the ablation identified stands, and the sources that could
relieve it are not these. The remaining candidates are primary-literature extraction (already
partly built under `MuriVerse/*/extractions`) and the records this cleanup still excludes for
recoverable reasons — dyes blocked only for an unstated attachment site, the proprietary
cyanines, and the backbone-altering residues (Sar, N-methyl, Hyp, pGlu) that would each need
their own construction and decoding path.
