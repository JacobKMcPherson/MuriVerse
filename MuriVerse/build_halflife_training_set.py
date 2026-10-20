"""Build one pooled mouse-and-rat in-vivo half-life model input table.

Current strict-tier input comes from `rodent/strict/in_vivo_half_life.csv`, whose records
are at most 50 residues and carry an explicit `structure_tier`. ESM-2 reads the residue
sequence; PeptideCLM reads a SMILES only when one passed the independent structure check.
Missing molecular features remain NaN. Annotation and termini blocks carry study context
and named partners. Mouse and rat records are pooled before homology clustering so a
shared peptide family cannot cross the saved train/validation split.

The discussion below records the earlier pooled-tier design. For the current strict tier,
the saved metadata, reconciliation, feature layout, and training-data manifest are the
authoritative descriptions of rows and columns.

Historical design notes:

PeptiVerse keeps representations in separate models — its half-life entry in
`best_models.txt` is `Halflife, Transformer, XGB (peptideclm), Regression`, a WT branch over
ESM-2 embeddings of the bare sequence and a SMILES branch over PeptideCLM embeddings of the
whole molecule, compared but never combined, and it never concatenates covariates onto an
embedding. MuriVerse deliberately does NOT follow that layout. It emits a single dataset
whose feature vector is four blocks concatenated:

    [ ESM-2 over the sequence | PeptideCLM over the PEPTIDE's SMILES | annotation | termini ]

The division of labour between those blocks is the main design decision in this module, and
it was arrived at by measurement rather than taste:

* **The SMILES is the peptide and nothing else.** It is a pure function of the amino-acid
  sequence, the topology (head-to-tail macrocycle, positioned disulfides) and the chirality.
  No PEG, no Fc, no albumin, not even an acetyl cap is drawn into it.
* **Everything attached is a feature**, in the `termini` block: a class per terminus, the
  heaviest partner's class, and two continuous masses (`peg_kda`, `partner_kda`).

Mixing the two broke the molecular encoder in two opposite ways, both measured:

* Dropping an undrawable partner but keeping the record made a peptide, that peptide +
  20 kDa PEG and that peptide + IgG Fc collapse onto ONE identical string. 47 strings mapped
  to more than one half-life-extension class over 219 records, one shared by 31 records
  spanning 6,400-fold of half-life — identical inputs, contradictory labels.
* Drawing the partner at its true size swamped the peptide. A 20 kDa PEG is ~1,360 heavy
  atoms against a ~300-atom peptide, and PeptideCLM pools length-weighted windows, so the
  embedding became mostly generic polyether (SMILES-only Spearman 0.247 -> 0.191).

With the split, two records share a SMILES only when they are genuinely the same peptide in
the same configuration, nothing inert outvotes the pharmacophore, and partner identity and
size reach the model directly as numbers instead of having to be inferred from atom counts.

A related defect is fixed by the same change. `smiles_coverage` used to double as an
accidental "has a big partner" flag — a record read as `peptide_core` mostly because it
carried an Fc — so it measured what RDKit could draw rather than the molecule, and it decayed
as the builder improved. It is now a statement about the PEPTIDE only (sequence, topology and
chirality all encodable), and partner presence is stated explicitly in `termini`.

`evaluate_halflife_model.py` scores every combination of the four blocks under cluster- and
study-grouped cross-validation, so which ones ship is a measurement, not an assertion.
Nothing is imputed and no record is dropped: where a block is absent it is written as NaN;
`StandardScaler` ignores NaN when fitting and preserves it, and XGBoost treats NaN as missing
and learns a default split direction.

The trainer (`training_classifiers/train_ml_regression.py`) is used unmodified apart from a
seeded Optuna sampler.

Two of the annotation covariates — `route_class` and `species` — are study-design variables
rather than molecular properties, so this model answers "what half-life will this molecule
show in this animal by this route", not "what half-life does this structure imply". Both are
known whenever a study is being planned, so this is a prediction rather than a leak, but a
structure-only model would have to drop them. See MODEL_COVARIATES.

The cohort is mouse and rat pooled, not mouse alone. Rat roughly doubles the record count,
and the two cohorts are drawn from the same database under the same filters, so they are
pooled rather than modelled separately. Three consequences are handled explicitly:

* `species` is a column on every row and on the saved dataset, the train/val split is
  balanced within each species (see `assign_splits`), and it is a one-hot FEATURE in the
  study_context block, so a pooled metric is neither secretly a metric for one animal nor
  blind to the offset between them, and a prediction for a mouse is not forced to equal
  the prediction for a rat.
* HMMER clustering runs over the pooled set, because sequences measured in both a mouse and
  a rat would otherwise sit in one species' train set and the other's validation set.
* The same peptide measured in both animals is kept as two records - `species` is part of
  the de-duplication key - since those pairs are the only direct evidence of the
  between-species offset in the set.

Splitting is homology-cluster-based, but the clustering tool has CHANGED. PeptiVerse's
`training_data_cleaned/wt_data_split.py` used MMseqs2; this now uses HMMER — `phmmer`
all-vs-all plus single-linkage over the accepted hits (see `cluster_sequences`, which records
how the two compare and where they disagree). With one model there is only one split, over
the full cohort; the earlier per-branch re-split, which gave the SMILES subset its own
non-comparable validation set, is gone.

One truncation remains and is not fixed here: ESM-2 is called with max_length=1022, so the
handful of records whose sequence is longer than that — an antibody heavy chain, albumin
fusions, a 3611-residue construct — are embedded from their first 1022 residues only. They
sit at the long-half-life end, so the loss leans the same way as the SMILES coverage gap, on
~2% of the cohort. Windowing them as the SMILES side is windowed would remove it.

  python MuriVerse/build_halflife_training_set.py --species mouse rat
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from datasets import Dataset, DatasetDict
from transformers import AutoModelForMaskedLM, AutoTokenizer, EsmModel

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
# The PeptideCLM SMILES tokenizer comes from the upstream PeptiVerse repository, which is
# not vendored here: scripts/fetch_tokenizer.sh clones it at a pinned commit.
PEPTIVERSE = Path(os.environ.get("PEPTIVERSE_DIR", REPO / "external/PeptiVerse"))
if not (PEPTIVERSE / "tokenizer/my_tokenizers.py").exists():
    sys.exit(f"PeptiVerse tokenizer not found under {PEPTIVERSE}; run scripts/fetch_tokenizer.sh "
             "or set PEPTIVERSE_DIR")
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(PEPTIVERSE))

import halflife_condition_curation as curation  # noqa: E402
import peptide_physchem  # noqa: E402
import pk_observations  # noqa: E402
from peptide_smiles import (  # noqa: E402
    build_smiles,
    describe_termini,
    looks_like_prose,
)
from tokenizer.my_tokenizers import SMILES_SPE_Tokenizer  # noqa: E402
from study_conditions import study_conditions  # noqa: E402
from training_data_index import write_index  # noqa: E402

ESM_MODEL = "facebook/esm2_t33_650M_UR50D"
PEPTIDECLM_MODEL = "aaronfeller/PeptideCLM-23M-all"
PEPTIDECLM_VOCAB = PEPTIVERSE / "tokenizer/new_vocab.txt"
PEPTIDECLM_SPLITS = PEPTIVERSE / "tokenizer/new_splits.txt"
MAX_LENGTH_PEPTIDECLM = 768

SEED = 1986
TRAIN_FRACTION = 0.8

# Clustering for the train/val split is HMMER-based (see cluster_sequences). An edge is kept
# when phmmer's per-domain i-Evalue is at or below HMMER_MAX_EVALUE and the aligned region
# covers at least HMMER_MIN_COVERAGE of BOTH sequences; clusters are connected components.
# 1e-3 is the conventional homology cut-off and 0.8 reproduces the bidirectional coverage that
# the previous `mmseqs cluster -c 0.8 --cov-mode 0` required.
HMMER_MAX_EVALUE = 1e-3
HMMER_MIN_COVERAGE = 0.8
HMMER_THREADS = 4

# Half-life values that are not a concentration-decay half-life of the dosed agent.
CORE_PK_PHASES = ("unlabelled", "terminal / beta")

BLANK = re.compile(r"^(|n\.?a\.?|none|not reported|not mentioned|nan)$", re.I)

# Ordered: first match wins, so the chemistry that dominates clearance is assigned.
MOD_BUCKETS = [
    ("depot / particle", r"microsphere|nanoparticle|nanoghost|liposom|hydrogel|implant|\bplga\b|"
                         r"micelle|depot|encapsulat"),
    ("protein fusion", r"\balbumin\b|\bhsa\b|\bmsa\b|\bfc\b|\bigg\b|\bgfp\b|repebody|his-?tag|"
                       r"transferrin|\belp\b|\bzipp\b|nanobody|antibody|\babd\b|\baff\b"),
    ("pegylation", r"\bpeg\b|peg\d|peg-|pegylat|polyethylene glycol"),
    ("polymer conjugate", r"chondroitin|heparosan|polysarcosine|\bhpn\b|\bpcl\b|dextran|hyaluron|"
                          r"polyglutam|\bpas\d|xten"),
    ("lipidation", r"myristo|palmito|palmitoyl|stearo|lauro|octanoy|decanoy|fatty acid|lipid|"
                   r"cholesterol|\bc16\b|\bc18\b"),
    ("label / imaging", r"fitc|cy\d|cy5|alexa|rhodamine|fluoroph|biotin|\bdota\b|do3a|nota|"
                        r"\b\d{2,3}[a-z]{1,2}\b\s*label|125i|131i|18f|64cu|68ga|99mtc|177lu|\b3h\b|"
                        r"radiolab|zw800|dbco"),
    ("acetylation", r"acetylat|\bac\b"),
    ("amidation", r"amidat|\bnh2\b|\bconh2\b"),
    ("esterification", r"esterif|methyl ester|\bome\b|\bote\b"),
    ("free", r"^free$"),
]


COVARIATES = ["n_ter_mod", "c_ter_mod", "topology", "chirality"]

# Study design and whole-molecule class. Everything that depends on what is ATTACHED to the
# peptide now lives in TERMINI_COVARIATES instead, so the ablation can price the two apart.
#
# Every one is bucketed from PEPlife2 annotation text and none is derived from the
# half-life, so this is not circular: `hl_engineering` and `route_class` come from
# subset_in_vivo_half_life.py's ENGINEERING/ROUTE regexes over Name | Origin | Chem-MOD |
# N-ter MOD | C-ter MOD and over Dose Route respectively.
#
# Two of them are STUDY-DESIGN variables, not molecular properties: `route_class` (how the
# animal was dosed) and `species`. Including them is a deliberate change of what the model
# answers - "what half-life will this molecule show in this animal by this route" rather
# than "what half-life does this structure imply" - and it is the framing the figure and
# caption state. They are available whenever a study is being planned, so this is a
# prediction, not a leak, but a structure-only model would have to drop them.
# `route_model` is the dose route with blanks filled from the sample, assay and name text;
# `sample_matrix`, `assay_method`, `strain` and `sex` are likewise read from the design columns
# (see study_conditions.py). They are study-design covariates, not molecular properties.
#
# `species` is one of them and is listed EXPLICITLY. It was absent from this list until the
# species fix, which made the claim above false in the only place it could be tested: the
# animal reached the matrix solely through `strain`, so two rows differing only in species
# but both carrying `strain=unknown` were identical inputs, and the shipped model returned
# one number for mouse and rat. `strain` is not a substitute. It is unknown on a large
# minority of the cohort, and on a molecule the cohort does not contain - which is what the
# shipped model is for - a caller knows the animal long before it knows the vendor strain.
# A single-species cohort drops the column again automatically, because _drop_constant
# removes any level that never varies.
STUDY_CONDITIONS = ["route_model", "species", "sample_matrix", "assay_method", "strain", "sex"]
MOLECULE_COVARIATES = ["hl_engineering", "topology", "chirality"]
STUDY_CONTEXT_COVARIATES = STUDY_CONDITIONS + ["formulation_class", "pk_dose_unit"]

# What is bolted onto the peptide, as its own block. This exists because the SMILES
# deliberately no longer carries it: `peptide_smiles.build_smiles` draws the sequence, the
# topology and the chirality and nothing else, so a 40 kDa PEG or an IgG Fc has to reach the
# model some other way. Drawing partners into the string was tried and failed twice over -
# dropping an undrawable partner collapsed chemically different drugs onto one identical
# SMILES, and drawing a real-size one swamped the peptide under length-weighted pooling.
# As features they do neither.
TERMINI_COVARIATES = ["n_ter_class", "c_ter_class", "partner_class"]

# How much of the PEPTIDE the SMILES represents - whether its sequence, topology and
# chirality were all encodable. No longer a statement about attached partners, which are not
# the string's job any more.
COVERAGE_LEVELS = ["complete", "partial", "none"]


def bucket_mod(value: object) -> str:
    text = str(value).strip()
    if BLANK.match(text):
        return "not reported"
    lowered = " ".join(text.lower().split())
    for label, pattern in MOD_BUCKETS:
        if re.search(pattern, lowered):
            return label
    return "other chemical"


# The 20 residues ESM-2 was pretrained on. Everything else in this cohort is either an
# ambiguity code (B = Asx, Z = Glx, X = any) or PEPlife2's letter for a non-natural residue
# (O for ornithine or for a transcription slip, U for a non-standard residue). ESM-2's
# vocabulary accepts all of them, so nothing crashes - they are simply embedded as rare
# tokens carrying almost no pretraining signal, which is a silent quality gradient across
# 16% of the cohort unless it is recorded.
STANDARD_RESIDUES = set("ACDEFGHIKLMNPQRSTVWY")


# PEPlife2 writes non-natural residues as hyphen-delimited three-letter codes inside an
# otherwise one-letter sequence: "H-Aib-QGTFTSDYSKYL...". Upper-casing and stripping
# non-letters turns each ONE residue into THREE standard ones - Aib becomes A, I, B; Nle
# becomes N, L, E; Sar becomes S, A, R - which lengthens the sequence, invents residues the
# peptide does not contain, and corrupts both the ESM-2 embedding and the homology cluster.
# 52 records carry Aib and 46 carry Sar, so this is not a rare case. Each code is collapsed
# to a single X, which keeps the length right and marks the position as unknown rather than
# asserting three wrong residues.
NONNATURAL_CODES = ("Aib", "Nle", "Orn", "Dab", "Sar", "Cha", "Abu", "Nal", "MeLeu", "hArg")
# Only a hyphen-delimited (or terminal) token is a code; a bare "AIB" inside a one-letter
# sequence is three real residues and must be left alone.
NONNATURAL = re.compile(
    r"(?:(?<=-)|^)(?:" + "|".join(NONNATURAL_CODES) + r")(?=-|$)", re.IGNORECASE)


def clean_sequence(value: object) -> str:
    return re.sub(r"[^A-Z]", "", NONNATURAL.sub("X", str(value)).upper())


def nonstandard_residues(sequence: str) -> str:
    """The non-standard letters in a sequence, as a sorted string ("" when there are none)."""
    return "".join(sorted(set(sequence) - STANDARD_RESIDUES))


def _phmmer_edges(fasta: Path, domtbl: Path, threads: int) -> list[tuple[str, str]]:
    """All-vs-all phmmer, returning the homology edges that pass E-value and coverage.

    phmmer is run permissively (`-E 1000`) and every threshold is applied here instead, so
    changing HMMER_MAX_EVALUE or HMMER_MIN_COVERAGE does not require re-running the search.
    """
    subprocess.run(
        ["phmmer", "--domtblout", str(domtbl), "-o", os.devnull,
         "-E", "1000", "--domE", "1000", "--incE", "1000", "--incdomE", "1000",
         "--cpu", str(threads), str(fasta), str(fasta)],
        check=True, capture_output=True, text=True)

    edges = []
    for line in domtbl.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        target, t_len, query, q_len = f[0], int(f[2]), f[3], int(f[5])
        if target == query:
            continue
        i_evalue = float(f[12])
        if i_evalue > HMMER_MAX_EVALUE:
            continue
        # In phmmer the query is the profile, so the "hmm" coordinates are the query's and the
        # "ali" coordinates the target's. Requiring both keeps the bidirectional coverage test
        # that MMseqs2's --cov-mode 0 applies.
        q_cov = (int(f[16]) - int(f[15]) + 1) / q_len
        t_cov = (int(f[18]) - int(f[17]) + 1) / t_len
        if q_cov >= HMMER_MIN_COVERAGE and t_cov >= HMMER_MIN_COVERAGE:
            edges.append((query, target))
    return edges


def cluster_sequences(ids: list[str], sequences: list[str], workdir: Path) -> dict[str, str]:
    """Single-linkage clustering over HMMER homology, replacing MMseqs2.

    PeptiVerse's `wt_data_split.py` used `mmseqs cluster --min-seq-id 0.3 -c 0.8`. This uses
    HMMER instead: `phmmer` searches the unique sequences against themselves, an edge is kept
    when the domain i-Evalue is at or below HMMER_MAX_EVALUE and the aligned region covers at
    least HMMER_MIN_COVERAGE of BOTH sequences, and the clusters are the connected components
    of that graph. The representative is the lexicographically smallest member, so the
    assignment is deterministic — MMseqs2's greedy set-cover representative was not.

    On the pooled rodent cohort (565 unique sequences) the two agree closely, and where they
    differ HMMER is the more conservative choice, which is what matters for a split whose only
    job is to keep homologues apart:

        co-clustered pairs   MMseqs2 1844   HMMER 4507   agreeing 1836
        merged by MMseqs2 only     8   (HMMER splits these)
        merged by HMMER only    2671   (HMMER merges these)

    The 8 pairs HMMER splits are all fusion or tandem-repeat constructs — exenatide dimers,
    interferon variants sharing only their signal peptide — where phmmer finds a strong hit
    (i-Evalue down to 1e-23) but only 12-73% of the molecule aligns, so the coverage test
    rejects it. MMseqs2 merged them transitively through a cluster representative rather than
    on pairwise coverage. HMMER_MIN_COVERAGE is the knob: lowering it to 0.5 recovers five of
    the eight.

    Two properties of HMMER to be aware of, neither of which applies to an identity threshold:

    * E-values depend on the size of the searched database, which here is the cohort itself.
      Adding records therefore shifts every E-value slightly, so a cluster assignment is only
      reproducible for a fixed cohort. `main()` keeps the search inputs and output in
      `<cohort>/hmmer/` (`all.fasta`, `phmmer.domtbl`) so the split that was actually used can
      be re-derived or re-thresholded without re-running phmmer.
    * phmmer cannot score very short sequences. 33 of the 565 unique sequences — every one of
      6 residues or fewer — produce no self-hit at all and so become singletons by default
      rather than by evidence. That is the safe direction for a split (a singleton cannot
      leak), and MMseqs2 left almost all of them as singletons too, but it means HMMER
      contributes no information about the short-peptide end of this cohort.

    Returns id -> cluster representative id; a sequence with no accepted edge maps to itself.
    """
    if shutil.which("phmmer") is None:
        raise RuntimeError("phmmer is not on PATH; install HMMER (>= 3.1) before running this "
                           "script, e.g. `apt-get install hmmer`.")

    workdir.mkdir(parents=True, exist_ok=True)
    fasta = workdir / "all.fasta"
    with fasta.open("w") as fh:
        for sid, seq in zip(ids, sequences):
            fh.write(f">{sid}\n{seq}\n")

    edges = _phmmer_edges(fasta, workdir / "phmmer.domtbl", HMMER_THREADS)

    parent = {sid: sid for sid in ids}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            # Union toward the smaller id so the representative is the component's minimum.
            if rb < ra:
                ra, rb = rb, ra
            parent[rb] = ra

    assignment = {sid: find(sid) for sid in ids}
    members: dict[str, int] = {}
    for rep in assignment.values():
        members[rep] = members.get(rep, 0) + 1
    print(f"  HMMER: {len(edges)} homology edges (i-Evalue <= {HMMER_MAX_EVALUE:g}, coverage "
          f">= {HMMER_MIN_COVERAGE}) -> {len(members)} clusters, largest {max(members.values())}, "
          f"{sum(1 for n in members.values() if n == 1)} singletons")
    return assignment


# No cluster larger than this share of the validation target may go to validation. The old
# greedy pass put the single largest homology family (a quarter of the records) in validation
# whole, so "validation" was one peptide family; a large family stays in training instead and
# validation is a spread of many small clusters.
MAX_VAL_CLUSTER_SHARE = 0.25


def assign_splits(df: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    """Whole clusters go to validation until each species' validation target is met.

    Clusters stay whole because that is the only thing keeping homologous sequences out of
    both sides of the split, and sequences measured in BOTH a mouse and a rat are common, so
    a cluster holding both species moves as one unit; splitting it by species would put the
    same sequence in train and in val.

    The target is applied per species rather than once over the pooled set, so the pooled
    validation metric is not a measurement of one animal reported as a metric for both. A
    cluster larger than MAX_VAL_CLUSTER_SHARE of the validation target is kept in training,
    which stops one family from becoming the whole validation set.
    """
    cluster_to_rows: dict[str, list[int]] = {}
    for row, cid in enumerate(df["cluster"]):
        cluster_to_rows.setdefault(cid, []).append(row)

    cluster_ids = list(cluster_to_rows)
    rng.shuffle(cluster_ids)

    species = (df["species"] if "species" in df.columns
               else pd.Series(["all"] * len(df), index=df.index)).to_numpy()
    val_targets = {sp: (1 - TRAIN_FRACTION) * n for sp, n in Counter(species).items()}
    cluster_cap = MAX_VAL_CLUSTER_SHARE * sum(val_targets.values())
    filled: dict[str, int] = {sp: 0 for sp in val_targets}

    split = np.full(len(df), "train", dtype=object)
    for cid in cluster_ids:
        rows = cluster_to_rows[cid]
        if len(rows) > cluster_cap:
            continue
        contrib = Counter(species[r] for r in rows)
        if all(filled[sp] + k <= val_targets[sp] for sp, k in contrib.items()):
            split[rows] = "val"
            for sp, k in contrib.items():
                filled[sp] += k
    return split


def build_meta(csv_paths: dict[str, Path], all_phases: bool,
               hmmer_dir: Path | None = None, strict: bool = False,
               species_scope: set[str] | None = None,
               ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pool one `in_vivo_half_life.csv` per species into a single training table.

    Pooling happens before clustering and before splitting, which is the point: the two
    cohorts share 44 sequences, so clustering them separately and concatenating afterwards
    would let the same peptide land in the mouse train set and the rat validation set.
    """
    frames = []
    for species, path in csv_paths.items():
        one = pd.read_csv(path, low_memory=False)
        # The column is written by subset_in_vivo_half_life.py; fall back for older files.
        one["species"] = one["species"] if "species" in one.columns else species
        print(f"  {species:6s} {len(one):5d} rows  <- {path.relative_to(ROOT)}")
        frames.append(one)
    df = pd.concat(frames, ignore_index=True)

    # Curated species corrections can leave a record in its source cohort file while assigning
    # it to another animal. Keep those rows in the pooled build, but not in a single-species run.
    # Keyed on the species asked for, not on the dict keys: the strict tier arrives as one
    # already-pooled file whose key is the tier name, so using the keys here would filter
    # every row out.
    selected_species = species_scope or set(csv_paths)
    in_scope = df["species"].isin(selected_species)
    out_of_scope = df.loc[~in_scope]
    df = df.loc[in_scope].copy()

    # Every filter below is recorded: the row count after it, how many records it removed, and
    # the PEPlife2 IDs of those records. `in_vivo_half_life.csv` is a published artefact and so
    # is the training set, and until this table existed the gap between the two counts - 513
    # records in the mouse file against 459 training rows - was not reconstructible from
    # anything either script emitted. `reconciliation` is returned so main() can write it out.
    reconciliation: list[dict] = []

    # Diffed on the ID column, not the index: the de-duplication step calls reset_index, after
    # which an index comparison silently matches unrelated rows. PEPlife2 IDs are unique across
    # the pooled frame (asserted below), so they are a safe key.
    if not df["ID"].is_unique:
        raise SystemExit("PEPlife2 IDs are not unique across the pooled cohort; the "
                         "reconciliation cannot identify removed records")

    def record(step: str, before: pd.DataFrame, after: pd.DataFrame, why: str) -> None:
        dropped = set(before["ID"]) - set(after["ID"])
        reconciliation.append({
            "step": step,
            "rows_in": len(before),
            "rows_out": len(after),
            "removed": len(before) - len(after),
            "why": why,
            "removed_ids": " ".join(str(i) for i in sorted(dropped)),
        })

    reconciliation.append({
        "step": "selected species", "rows_in": len(pd.concat(frames, ignore_index=True)),
        "rows_out": len(df), "removed": len(out_of_scope),
        "why": "keep records whose curated species is in the requested cohort(s): "
               + ", ".join(sorted(selected_species)),
        "removed_ids": " ".join(str(i) for i in sorted(out_of_scope["ID"]))})

    if not all_phases:
        before = df
        df = df[df["pk_phase"].isin(CORE_PK_PHASES)]
        record("concentration-decay pk_phase", before, df,
               f"keeps pk_phase in {list(CORE_PK_PHASES)}; drops distribution/alpha and "
               f"non-concentration readouts, which answer a different question")

    before = df
    if strict:
        # The strict tier ships two sequence columns and they are both needed.
        # `sequence_tokens` is the true residue sequence, with a non-natural residue written
        # as its own hyphen-delimited code and a D residue in lower case; it is what the
        # verified SMILES was built from. `sequence_esm` is that same sequence as ESM-2 can
        # actually consume: upper-cased, and with an X wherever the residue is one ESM-2 was
        # never pretrained on. ESM-2 and HMMER both get the latter, because neither has a
        # token for Aib; `nonstandard_positions` records exactly where the X's are, so the
        # resulting quality gradient is measurable rather than hidden.
        missing = df["sequence_esm"].isna() | (df["sequence_esm"] == "")
        if missing.any():
            raise SystemExit(f"{int(missing.sum())} strict rows have no sequence_esm; "
                             f"rebuild with build_strict_peptide_set.py")
        df = df.copy()
        df["sequence"] = df["sequence_esm"]
        record("has a sequence", before, df,
               "the strict tier guarantees one, so nothing is dropped here")
    else:
        df = df[df["sequence"].notna()].copy()
        df["sequence"] = df["sequence"].map(clean_sequence)
        df = df[df["sequence"].str.len() > 0]
    if not strict:
        record("has a sequence", before, df,
               "PEPlife2 leaves Sequence empty on some records; a sequence model has no "
               "input for them")

    before = df
    df["half_life_h"] = pd.to_numeric(df["half_life_h"], errors="coerce")
    df = df[df["half_life_h"] > 0]
    record("positive half_life_h", before, df,
           "log10 of the label requires a strictly positive half-life")

    if strict:
        # --- parenteral dose route ------------------------------------------------------
        # build_strict_peptide_set.py has already curated `route_class` against the primary
        # sources wherever PEPlife2's `Dose Route` column was blank, and recorded how it
        # knows in `route_evidence_level`. What is left is the cohort decision.
        #
        # The cohort keeps intravenous, subcutaneous, intraperitoneal and intramuscular
        # records only. An oral, intranasal or rectal half-life is an absorption-limited
        # composite - flip-flop kinetics make the terminal slope the ABSORPTION rate, not
        # elimination - so it answers a different question from the one the model is fitted
        # on, and mixing the two teaches the model that a route column it can read is worth
        # more than the molecule. A record whose route the sources do not establish is
        # dropped for the same reason: nothing says which of the two it is.
        before = df
        df = df[df["route_class"].isin(curation.PARENTERAL)]
        record("parenteral dose route", before, df,
               f"keeps route_class in {list(curation.PARENTERAL)} after source curation; "
               f"drops enteral and mucosal routes, whose terminal slope is absorption, "
               f"and records whose route no primary source establishes")

        # --- whole-molecule structure -------------------------------------------------
        # Every kept record has a VERIFIED SMILES, but not every SMILES is the whole dosed
        # molecule: `partly resolved` leaves a reported ring or configuration unspecified
        # and `modification not drawn` emits the peptide WITHOUT a clause that alters it.
        # Both feed the molecular encoder a string that is not the molecule the label
        # belongs to, so the cohort now requires the complete tier. build_strict_peptide_set
        # resolves what can be resolved first; what remains is listed with its reason in
        # the kept file's `unresolved` column.
        before = df
        df = df[df["structure_tier"].eq("complete")]
        record("whole-molecule structure", before, df,
               "keeps structure_tier == complete; drops records whose SMILES omits a "
               "reported modification or leaves a reported ring or configuration "
               "unspecified, so the string is not the dosed molecule")

        # --- a named bioanalytical method ----------------------------------------------
        # The study-condition covariates are derived here rather than after the chain, so
        # the assay is available to filter on and is computed once.
        df = df.reset_index(drop=True)
        conditions = study_conditions(df, df["route_class"])
        for column in conditions:
            df[column] = conditions[column].to_numpy()
        # `Assay` names a separation rather than a detector on some records ("TLC"), or the
        # readout rather than the method ("Serum decay experiment"), or nothing at all.
        # Where the vocabulary comes out "unknown", the method the source states is filled
        # in from the assay ledger; `unknown` survives only where no source names one.
        assay = curation.resolve_assay(df, conditions["assay_method"])
        for column in assay:
            df[column] = assay[column].to_numpy()

        # What quantity decayed is not a detail of a half-life, it IS the half-life, and a
        # record whose method no source names does not say which quantity that was. The
        # cohort's own numbers make the point: the six records left reading `unknown` have
        # a median half-life of 23 h against 0.43 h for the set, and three of them are the
        # three longest labels in it (41, 58 and 76 h). Keeping them puts the largest
        # values in the label range behind the one covariate that cannot be checked, so
        # they are dropped rather than carried as a level.
        before = df
        df = df[df["assay_method"].ne("unknown")]
        record("named bioanalytical method", before, df,
               "drops records whose assay no source names; an unidentified method leaves "
               "the measured quantity unknown, and those records sit at the long end of "
               "the label range")

    df["n_ter_mod"] = df["N-ter MOD"].map(bucket_mod)
    df["c_ter_mod"] = df["C-ter MOD"].map(bucket_mod)
    df["chirality_reported"] = df["Chirality"].map(
        lambda v: "not reported" if BLANK.match(str(v).strip()) else str(v).strip().lower()
    )
    df["chirality"] = df["chirality_reported"]
    if strict:
        # The strict builder round-trip verifies stereochemistry against PEPlife2's original
        # sequence case. Use that verified configuration as the model feature: PEPlife2's
        # summary Chirality field contradicts the sequence for three Apraglutide rows.
        df["chirality"] = df["chirality_from_sequence"].str.lower()
        if df["chirality"].isna().any() or not df["chirality"].isin({"l", "d", "mix"}).all():
            raise ValueError("Strict rows need sequence-verified chirality (L, D or mixed)")

    # Preserve source observations. Even rows with the same molecule and label may be
    # different dose, route, formulation, strain or assay arms. The strict-set builder marks
    # possible repeats in metadata for review rather than inferring that they are duplicates.
    df = df.reset_index(drop=True)
    record("observations retained", df, df,
           "retain each positive-label source record; repeated-observation candidates are "
           "flagged for review and are not dropped")

    rec = pd.DataFrame(reconciliation)
    print(rec[["step", "rows_in", "rows_out", "removed"]].to_string(index=False))

    df["label"] = df["half_life_h"].astype(np.float32)
    df["log_label"] = np.log10(df["half_life_h"]).astype(np.float32)
    df["id"] = [f"hl_{i:05d}" for i in range(len(df))]

    if strict:
        # The strict tier's string arrives already built, round-trip decoded and
        # mass-checked by build_strict_peptide_set.py, and it carries the side-chain
        # conjugates that module resolves from Chem-MOD. Rebuilding it here would throw
        # those away, because this module reads only the two terminal columns.
        # Rows in the `sequence only` tier have no structure by definition. They are kept
        # because the sequence branch can still read them, and their SMILES block is zeroed
        # below rather than filled with a string that is not the molecule.
        missing = df["smiles"].isna() | (df["smiles"].astype(str).str.strip() == "")
        tier = df.get("structure_tier", pd.Series("complete", index=df.index))
        unexpected = missing & tier.ne("sequence only")
        if unexpected.any():
            raise SystemExit(f"{int(unexpected.sum())} rows outside the 'sequence only' "
                             f"tier have no SMILES; rebuild the strict set first")
        df["smiles"] = df["smiles"].fillna("")
        df["smiles_status"] = np.where(missing, "no structure generated",
                                       "verified by build_strict_peptide_set")
        # `smiles_coverage` is how much of the PEPTIDE the string describes, which is what
        # the model reads as a quality gradient over these rows.
        df["smiles_coverage"] = np.select(
            [missing, tier.isin(["partly resolved", "modification not drawn"])],
            ["none", "partial"], default="complete")
        df["smiles_stereo_unspecified"] = tier.eq("partly resolved").to_numpy()
        df["smiles_topology_encoded"] = ~missing.to_numpy()
    else:
        # NB all five arguments: the terminal columns decide the caps, and omitting them
        # raised a TypeError rather than quietly drawing a bare peptide.
        built = [build_smiles(row["sequence"], row["N-ter MOD"], row["C-ter MOD"],
                              row["Linear/Cyclic"], row["Chirality"])
                 for _, row in df.iterrows()]
        df["smiles"] = [b.smiles for b in built]
        df["smiles_status"] = [b.status for b in built]
        # Coverage of the PEPTIDE, not of the intact drug: 'complete' = sequence, topology
        # and chirality are all encoded; 'partial' = drawn but the topology could not be
        # positioned or the L/D pattern is unresolved; 'none' = no string. Attached
        # partners are not part of this judgement - they are the termini block's job.
        df["smiles_coverage"] = [b.coverage for b in built]
        df["smiles_stereo_unspecified"] = [b.stereo_unspecified for b in built]
        df["smiles_topology_encoded"] = [b.topology_encoded for b in built]

    if strict:
        # The strict set is now TWO tiers. In the `complete` tier the SMILES is the whole
        # dosed molecule and there is no partner to describe. In the `partner carried` tier
        # the SMILES is the whole PEPTIDE and something is attached to it, and that
        # something has to reach the model somehow - which is what this block is for.
        # Describing the termini for every row keeps one schema, and leaves `partner_class`
        # at "none" for the complete tier because for those rows it is true.
        termini = [describe_termini(row["N-ter MOD"], row["C-ter MOD"])
                   for _, row in df.iterrows()]
        df["n_ter_class"] = [t.n_class for t in termini]
        df["c_ter_class"] = [t.c_class for t in termini]
        carried = df.get("structure_tier", pd.Series("complete", index=df.index)) \
                    .eq("partner carried")
        df["partner_class"] = [t.partner_class if c else "none"
                               for t, c in zip(termini, carried)]
        df["peg_kda"] = [t.peg_kda if c else 0.0 for t, c in zip(termini, carried)]
        df["partner_kda"] = [t.partner_kda if c else 0.0 for t, c in zip(termini, carried)]
        # `describe_termini` reports the terminal CLASS, which for a carried row whose
        # terminus is a plain cap comes back as "amidation" or "acetylation". A cap is not a
        # partner; those rows carry theirs somewhere else - in `hl_engineering`, in `Name`,
        # or in a Chem-MOD clause - and are labelled for what is known about them rather
        # than for a terminal modification that has nothing to do with it.
        NOT_A_PARTNER = {"none", "amidation", "acetylation", "free"}
        if "structure_tier" in df:
            unnamed = carried & df["partner_class"].isin(NOT_A_PARTNER)
            df.loc[unnamed, "partner_class"] = "other partner"
        df["sequence_is_prose"] = False
        df["sequence_nonstandard"] = [nonstandard_residues(s) for s in df["sequence"]]
        # The length that matters is the PARSED residue count, which build_strict_peptide_set
        # recorded. The upstream selection could only count letters, and a delimited
        # non-standard residue is several letters, so that figure over-counts.
        if "n_residues" in df:
            df["seq_length"] = df["n_residues"].astype(int)
        df = df.reset_index(drop=True)
        # The study-condition covariates, including the curated assay method, were derived
        # in the filter chain above - the assay is one of the cohort requirements, so it has
        # to exist before the row is kept. Deriving them again here would be a second
        # evaluation of the same deterministic function on the same rows.
        missing = [c for c in ("sample_matrix", "assay_method", "strain", "sex",
                               "route_model") if c not in df.columns]
        if missing:
            raise SystemExit(f"study conditions missing from the strict frame: {missing}")
        return _cluster_and_split(df, rec, hmmer_dir)

    termini = [describe_termini(row["N-ter MOD"], row["C-ter MOD"])
               for _, row in df.iterrows()]
    df["n_ter_class"] = [t.n_class for t in termini]
    df["c_ter_class"] = [t.c_class for t in termini]
    df["partner_class"] = [t.partner_class for t in termini]
    df["peg_kda"] = [t.peg_kda for t in termini]
    df["partner_kda"] = [t.partner_kda for t in termini]
    # The sequence column itself is prose in a handful of rows, which invalidates that
    # row's length, its ESM-2 embedding and its homology cluster, not merely its SMILES.
    df["sequence_is_prose"] = [looks_like_prose(s) for s in df["sequence"]]
    df["sequence_nonstandard"] = [nonstandard_residues(s) for s in df["sequence"]]

    return _cluster_and_split(df, rec, hmmer_dir)


def _cluster_and_split(df: pd.DataFrame, rec: pd.DataFrame,
                       hmmer_dir: Path | None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Homology-cluster the pooled sequences and assign the train/val split.

    Shared by both tiers, so the two are split the same way: HMMER over the cohort being
    built, then single-linkage, then a split balanced within each species.
    """
    unique = sorted(df["sequence"].unique())
    seq_ids = [f"seq_{i:05d}" for i in range(len(unique))]
    # Kept rather than thrown away in a temp dir: phmmer E-values depend on the size of the
    # searched database, which here is the cohort itself, so a cluster assignment is only
    # reproducible for a fixed cohort. all.fasta plus phmmer.domtbl pin the exact search that
    # produced this split.
    if hmmer_dir is None:
        with tempfile.TemporaryDirectory() as tmp:
            assignment = cluster_sequences(seq_ids, unique, Path(tmp) / "hmmer")
    else:
        assignment = cluster_sequences(seq_ids, unique, hmmer_dir)
    id_to_seq = dict(zip(seq_ids, unique))
    seq_to_cluster = {seq: id_to_seq[assignment[sid]] for sid, seq in zip(seq_ids, unique)}

    df["cluster"] = df["sequence"].map(seq_to_cluster)
    df["split"] = assign_splits(df, np.random.default_rng(SEED))
    return df, rec


@torch.no_grad()
def _pool(model, inputs, special_ids: torch.Tensor | None) -> np.ndarray:
    hidden = model(**inputs).last_hidden_state
    valid = inputs["attention_mask"].bool()
    if special_ids is not None:
        valid = valid & ~torch.isin(inputs["input_ids"], special_ids)
    weights = valid.unsqueeze(-1).float()
    pooled = (hidden * weights).sum(dim=1) / weights.sum(dim=1).clamp(min=1e-9)
    return pooled.cpu().numpy().astype(np.float32)


@torch.no_grad()
def embed_esm(sequences: list[str], device: torch.device, batch_size: int,
              cache_path: Path | None = None) -> np.ndarray:
    """Mean-pooled ESM-2 last hidden state, computed once per unique sequence.

    ESM-2 650M over this cohort is several CPU-minutes and its output depends only on the
    sequence, so results are memoised to `cache_path` across runs. The key is the sequence
    itself, which makes a stale entry impossible as long as ESM_MODEL is unchanged — the
    file records the model id and is ignored if it disagrees.
    """
    cached: dict[str, np.ndarray] = {}
    if cache_path is not None and cache_path.exists():
        blob = np.load(cache_path, allow_pickle=False)
        if str(blob["model"]) == ESM_MODEL:
            cached = dict(zip(blob["sequences"].tolist(), blob["vectors"]))
            print(f"  ESM cache: {len(cached)} sequences from {cache_path.name}")

    unique = sorted(set(sequences) - set(cached), key=len)
    if unique:
        tokenizer = AutoTokenizer.from_pretrained(ESM_MODEL)
        model = EsmModel.from_pretrained(ESM_MODEL).to(device).eval()
        for start in range(0, len(unique), batch_size):
            batch = unique[start:start + batch_size]
            inputs = tokenizer(batch, padding=True, truncation=True, max_length=1022,
                              return_tensors="pt")
            inputs = {k: v.to(device) for k, v in inputs.items()}
            for seq, vec in zip(batch, _pool(model, inputs, None)):
                cached[seq] = vec
            print(f"  ESM {min(start + batch_size, len(unique))}/{len(unique)}")
        if cache_path is not None:
            keys = sorted(cached)
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            np.savez(cache_path, model=np.array(ESM_MODEL),
                     sequences=np.array(keys, dtype=object).astype("U"),
                     vectors=np.stack([cached[k] for k in keys]))

    return np.stack([cached[s] for s in sequences])


def peptideclm_tokenizer() -> SMILES_SPE_Tokenizer:
    return SMILES_SPE_Tokenizer(str(PEPTIDECLM_VOCAB), str(PEPTIDECLM_SPLITS))


@torch.no_grad()
def embed_peptideclm(smiles: list[str], device: torch.device, batch_size: int,
                     cache_path: Path | None = None) -> np.ndarray:
    """Mean-pooled PeptideCLM, over sliding windows when the molecule exceeds the context.

    PeptideCLM's context is 768 tokens and a third of the drawable molecules here are
    longer — up to ~25k tokens for the largest fusion backbones. Truncating keeps only the
    N-terminal fragment, and dropping those rows removes most of the long-half-life end of
    the label distribution. Pooling over consecutive windows and averaging by window length
    keeps every atom in the representation; it is a weaker encoding of long molecules than
    a native long-context pass, but it is not a biased sample of them.

    Memoised per unique SMILES, as `embed_esm` is per sequence: the windowed pass over the
    long fusion backbones dominates the build, and the result depends only on the string.
    The cache is keyed by SHA-1 of the SMILES rather than by the string itself - numpy pads
    every entry of a unicode array to the longest one, and the longest molecule here is
    ~25k characters, which turned a 3 MB cache into a 130 MB one.
    """
    cached: dict[str, np.ndarray] = {}
    if cache_path is not None and cache_path.exists():
        blob = np.load(cache_path, allow_pickle=False)
        if str(blob["model"]) == PEPTIDECLM_MODEL:
            cached = dict(zip(blob["keys"].tolist(), blob["vectors"]))
            print(f"  PeptideCLM cache: {len(cached)} molecules from {cache_path.name}")

    digest = {s: hashlib.sha1(s.encode()).hexdigest() for s in set(smiles)}
    todo = sorted({s for s in smiles if digest[s] not in cached}, key=len)
    if todo:
        fresh = _embed_peptideclm_uncached(todo, device, batch_size)
        cached.update({digest[s]: v for s, v in fresh.items()})
        if cache_path is not None:
            keys = sorted(cached)
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            np.savez(cache_path, model=np.array(PEPTIDECLM_MODEL),
                     keys=np.array(keys, dtype="U40"),
                     vectors=np.stack([cached[k] for k in keys]))

    return np.stack([cached[digest[s]] for s in smiles])


@torch.no_grad()
def _embed_peptideclm_uncached(smiles: list[str], device: torch.device,
                               batch_size: int) -> dict[str, np.ndarray]:
    tokenizer = peptideclm_tokenizer()
    model = AutoModelForMaskedLM.from_pretrained(PEPTIDECLM_MODEL).roformer.to(device).eval()
    ids_of = {a: getattr(tokenizer, a, None) for a in ("cls_token_id", "sep_token_id",
                                                       "pad_token_id")}
    special_ids = torch.tensor(
        sorted({getattr(tokenizer, a, None) for a in
                ("pad_token_id", "cls_token_id", "sep_token_id",
                 "bos_token_id", "eos_token_id", "mask_token_id")} - {None}),
        device=device, dtype=torch.long,
    )
    cls_id, sep_id = ids_of["cls_token_id"], ids_of["sep_token_id"]
    pad_id = ids_of["pad_token_id"] or 0
    wrap = (cls_id is not None) + (sep_id is not None)
    body = MAX_LENGTH_PEPTIDECLM - wrap

    def windows(text: str) -> list[list[int]]:
        ids = tokenizer(text)["input_ids"]
        if cls_id is not None and ids and ids[0] == cls_id:
            ids = ids[1:]
        if sep_id is not None and ids and ids[-1] == sep_id:
            ids = ids[:-1]
        chunks = [ids[i:i + body] for i in range(0, max(len(ids), 1), body)] or [[]]
        return [([cls_id] if cls_id is not None else []) + c
                + ([sep_id] if sep_id is not None else []) for c in chunks]

    # One flat list of windows, so batching stays even across very unequal molecules.
    owner, flat = [], []
    for i, text in enumerate(smiles):
        for w in windows(text):
            owner.append(i)
            flat.append(w)

    pooled = np.zeros((len(flat), model.config.hidden_size), dtype=np.float32)
    for start in range(0, len(flat), batch_size):
        batch = flat[start:start + batch_size]
        width = max(len(w) for w in batch)
        input_ids = torch.full((len(batch), width), pad_id, dtype=torch.long)
        mask = torch.zeros((len(batch), width), dtype=torch.long)
        for r, w in enumerate(batch):
            input_ids[r, :len(w)] = torch.tensor(w, dtype=torch.long)
            mask[r, :len(w)] = 1
        inputs = {"input_ids": input_ids.to(device), "attention_mask": mask.to(device)}
        pooled[start:start + len(batch)] = _pool(model, inputs, special_ids)
        print(f"  PeptideCLM window {min(start + batch_size, len(flat))}/{len(flat)}")

    # Average each molecule's windows, weighted by how many real tokens each contributed.
    weights = np.array([len(w) - wrap for w in flat], dtype=np.float32).clip(min=1.0)
    out = np.zeros((len(smiles), pooled.shape[1]), dtype=np.float32)
    total = np.zeros(len(smiles), dtype=np.float32)
    for row, mol in enumerate(owner):
        out[mol] += pooled[row] * weights[row]
        total[mol] += weights[row]
    out = out / total[:, None].clip(min=1e-9)
    return dict(zip(smiles, out))


def _one_hot(meta: pd.DataFrame, columns: list[str]) -> tuple[list[np.ndarray], list[str]]:
    """Levels taken in sorted order from the data, so the column layout is reproducible."""
    blocks, names = [], []
    for col in columns:
        values = meta[col].astype(str).fillna("not reported")
        for level in sorted(values.unique()):
            names.append(f"{col}={level}")
            blocks.append((values == level).to_numpy(np.float32))
    return blocks, names


def _drop_constant(blocks: list[np.ndarray], columns: list[str]) -> tuple[np.ndarray, list[str]]:
    """Stack the columns, leaving out any that never vary. A one-level column carries no
    information, and in a tree model only adds a split that can never be taken."""
    keep = []
    for i, block in enumerate(blocks):
        finite = block[np.isfinite(block)]
        if finite.size and np.ptp(finite) > 0:
            keep.append(i)
    return (np.stack([blocks[i] for i in keep], axis=1).astype(np.float32),
            [columns[i] for i in keep])


def molecule_annotation_block(meta: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    """Molecular classes and coverage, kept apart from experimental context."""
    blocks, columns = _one_hot(meta, MOLECULE_COVARIATES)
    for level in COVERAGE_LEVELS:
        columns.append(f"smiles_coverage={level}")
        blocks.append((meta["smiles_coverage"] == level).to_numpy(np.float32))
    return _drop_constant(blocks, columns)


def study_context_block(meta: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    """Animal, administration, assay, formulation and reported dose/age context."""
    blocks, columns = _one_hot(meta, STUDY_CONTEXT_COVARIATES)
    for name, source in (("log10_dose_mg_kg", "pk_dose_log10_mg_kg"),
                         ("age_years", "pk_age_years"),
                         ("body_weight_g", "pk_body_weight_g")):
        values = pd.to_numeric(meta.get(source, pd.Series(np.nan, index=meta.index)),
                               errors="coerce").to_numpy(np.float32)
        columns.append(name)
        blocks.append(values)
    return _drop_constant(blocks, columns)


def termini_block(meta: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    """What is attached to the peptide: class per terminus, and the two partner masses.

    This is the block that replaces drawing partners into the SMILES. The masses are kept
    continuous and log-scaled rather than bucketed because clearance falls steeply and
    smoothly with hydrodynamic size - the difference between 1 and 5 kDa of PEG matters more
    than the difference between 35 and 40.
    """
    blocks, columns = _one_hot(meta, TERMINI_COVARIATES)
    for name, source in (("log1p_peg_kda", "peg_kda"), ("log1p_partner_kda", "partner_kda")):
        columns.append(name)
        blocks.append(np.log1p(meta[source].fillna(0.0).to_numpy(np.float32)))
    for name, source in (("has_n_mod", "n_ter_class"), ("has_c_mod", "c_ter_class")):
        columns.append(name)
        blocks.append((meta[source] != "free").to_numpy(np.float32))
    return _drop_constant(blocks, columns)


def save_dataset(meta: pd.DataFrame, features: np.ndarray, out_dir: Path,
                 extra: dict[str, list] | None = None) -> DatasetDict:
    full = Dataset.from_dict({
        "id": meta["id"].tolist(),
        "sequence": meta["sequence"].tolist(),
        "embedding": features,
        "label": meta["log_label"].tolist(),
        "raw_label": meta["label"].tolist(),
        "split": meta["split"].tolist(),
        **(extra or {}),
    })
    ds = DatasetDict({
        "train": full.filter(lambda x: x["split"] == "train").remove_columns("split"),
        "val": full.filter(lambda x: x["split"] == "val").remove_columns("split"),
    })
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    ds.save_to_disk(str(out_dir))
    return ds


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", nargs="+", default=["mouse", "rat"],
                    choices=["mouse", "rat"],
                    help="Cohorts to pool. More than one builds the rodent training set.")
    ap.add_argument("--cohort", default="",
                    help="Name for the pooled output; defaults to the species when one is "
                         "given and 'rodent' when several are.")
    ap.add_argument("--tier", default="strict", choices=["strict", "pooled"],
                    help="'strict' (the default) trains on rodent/strict/"
                         "in_vivo_half_life.csv - peptide-scale records whose whole "
                         "molecule is drawn and verified, built by "
                         "build_strict_peptide_set.py. 'pooled' trains on the per-species "
                         "in_vivo_half_life.csv files, which include Fc fusions, albumin "
                         "constructs, polymer conjugates and kDa-scale PEG, whose "
                         "clearance is the partner's rather than the peptide's.")
    ap.add_argument("--all_phases", action="store_true",
                    help="Keep distribution/alpha and non-concentration half-lives too.")
    ap.add_argument("--batch_size", type=int, default=4)
    ap.add_argument("--blocks", nargs="+", default=None,
                    choices=["sequence", "smiles", "molecule_annotation", "study_context",
                             "termini", "physchem"],
                    help="Feature blocks for the shipped dataset. The ablation dataset "
                         "always carries all six.")
    args = ap.parse_args()
    if args.blocks is None:
        args.blocks = ["sequence", "smiles", "molecule_annotation", "study_context",
                       "termini", "physchem"]

    species = sorted(dict.fromkeys(args.species))
    cohort = args.cohort or (species[0] if len(species) == 1 else "rodent")
    strict = args.tier == "strict"
    # The strict tier is one already-pooled file, so it gets its own cohort directory
    # rather than overwriting the pooled tier's split, clusters and reconciliation.
    cohort_dir = ROOT / cohort / "strict" if strict else ROOT / cohort
    cohort_dir.mkdir(parents=True, exist_ok=True)

    print(f"[Cohort] {cohort} ({args.tier}) <- {', '.join(species)}")
    if strict:
        strict_path = ROOT / "rodent" / "strict" / "in_vivo_half_life.csv"
        if not strict_path.exists():
            raise SystemExit(f"{strict_path} is missing; run "
                             f"build_strict_peptide_set.py first")
        sources = {"strict": strict_path}
    else:
        sources = {sp: ROOT / sp / "in_vivo_half_life.csv" for sp in species}
    meta, reconciliation = build_meta(sources, args.all_phases,
                                      hmmer_dir=cohort_dir / "hmmer", strict=strict,
                                      species_scope=set(species))

    # The published per-species record counts must add up to the reconciliation's first row,
    # so the chain from in_vivo_half_life.csv to the training set cannot silently drift if a
    # filter is added upstream without being accounted for here.
    published = sum(len(pd.read_csv(path, low_memory=False)) for path in sources.values())
    entry = int(reconciliation.iloc[0]["rows_in"])
    if entry != published:
        raise SystemExit(f"reconciliation starts at {entry} rows but the source files hold "
                         f"{published}; the chain is not accounted for")
    if int(reconciliation.iloc[-1]["rows_out"]) != len(meta):
        raise SystemExit("reconciliation does not end at the training-set row count")

    pk_table_path = ROOT / "rodent" / "pk_observations.csv"
    pk_table = pk_observations.write_table(ROOT, pk_table_path)
    meta = pk_observations.enrich_training_rows(meta, pk_table)
    print(f"[Full-text evidence] {int(meta['pk_observation_matched'].sum())} exact-ID "
          f"treatment-arm observations joined; "
          f"{int(meta['label_evidence_tier'].eq('full_text_conflict').sum())} label conflicts")

    recon_path = cohort_dir / "halflife_cohort_reconciliation.csv"
    reconciliation.to_csv(recon_path, index=False)
    print(f"[Reconciliation] {published} published records -> {len(meta)} training rows"
          f"  ({recon_path.relative_to(ROOT)})")

    tokenizer = peptideclm_tokenizer()
    meta["smiles_tokens"] = [len(tokenizer(s)["input_ids"]) if isinstance(s, str) else 0
                             for s in meta["smiles"]]
    # Over-context molecules are no longer reclassified and dropped: embed_peptideclm pools
    # them over sliding windows instead. The token count stays as a column so the figures
    # can still show how much of the cohort needs more than one window.
    meta["smiles_windows"] = np.ceil(
        meta["smiles_tokens"].clip(lower=1) / MAX_LENGTH_PEPTIDECLM).astype(int)

    meta_path = cohort_dir / "halflife_meta_with_split.csv"
    strict_columns = ["sequence_tokens", "sequence_esm", "nonstandard_positions",
                      "n_nonstandard", "chirality_from_sequence", "stereo_source",
                      "conjugates", "tracer_label", "structure_tier", "partner_text",
                      "partner_source", "n_residues", "unresolved",
                      "half_life_extreme_qc", "same_study_structure_condition_n",
                      "same_study_structure_condition_label_n",
                      "same_structure_label_across_sources_n", "exact_duplicate_candidate",
                      "within_context_repeated_measurement",
                      "cross_source_repeated_measurement",
                      "repeat_observation_candidate", "repeat_candidate_scope",
                      "label_evidence_tier", "pk_observation_matched",
                      "pk_dose_text", "pk_dose_log10_mg_kg", "pk_dose_unit",
                      "pk_age_years", "pk_body_weight_g", "pk_group_n",
                      "pk_condition_conflicts", "formulation_class",
                      "evidence_agent_name", "evidence_route", "evidence_formulation",
                      "evidence_matrix", "evidence_sampling_design",
                      "evidence_bioanalytical_method", "evidence_half_life",
                      "evidence_half_life_phase", "evidence_cmax", "evidence_tmax",
                      "evidence_auc", "evidence_clearance",
                      "evidence_volume_of_distribution", "evidence_mrt",
                      "evidence_bioavailability", "pk_evidence_json"]
    extra_columns = [c for c in strict_columns if c in meta.columns] if strict else []
    for col in ("half_life_qualifier", "pk_label", "route_class",
                "route_evidence_level", "route_evidence_quote",
                "assay_evidence_source", "assay_evidence_quote"):
        if col not in meta.columns:
            meta[col] = ""
    meta["peplife_id"] = meta["ID"]
    meta[["id", "peplife_id", "species", "PMID", "Name", "sequence", "seq_length", "n_ter_mod", "c_ter_mod", "topology",
          "chirality", "chirality_reported", "half_life_h", "half_life_qualifier", "half_life_provenance", "label",
          "log_label", "pk_phase", "pk_label",
          "route_class", "route_model", "sample_matrix", "assay_method", "strain", "sex",
          # How the route and the assay came to be what they are. Shipped beside the values
          # so a reader can tell a source-curated route from a transcribed one without
          # going back to the ledgers.
          "route_evidence_level", "route_evidence_quote",
          "assay_evidence_source", "assay_evidence_quote",
          "database_route_model", "database_sample_matrix", "database_assay_method",
          "database_strain", "database_sex",
          "hl_engineering", "smiles", "smiles_status", "smiles_coverage",
          "smiles_stereo_unspecified", "smiles_topology_encoded",
          "n_ter_class", "c_ter_class", "partner_class", "peg_kda", "partner_kda",
          "sequence_is_prose", "sequence_nonstandard", "smiles_tokens", "smiles_windows",
          "cluster", "split"] + extra_columns].to_csv(meta_path, index=False)

    print(f"\n[Split] {meta_path}")
    print(f"  rows {len(meta)}  sequences {meta['sequence'].nunique()}"
          f"  HMMER clusters {meta['cluster'].nunique()}")
    print(f"  {meta['split'].value_counts().to_dict()}")
    print("  split x species:")
    print(pd.crosstab(meta["species"], meta["split"]).to_string().replace("\n", "\n    "))
    for col in COVARIATES:
        print(f"  {col}: {meta[col].value_counts().to_dict()}")
    flagged = meta["sequence_nonstandard"].ne("")
    letters = Counter(c for s in meta.loc[flagged, "sequence_nonstandard"] for c in s)
    print(f"  sequences with non-standard residues: {int(flagged.sum())} rows "
          f"({dict(sorted(letters.items()))} by letter)")
    print("\n[SMILES] reconstruction outcome")
    print(meta["smiles_status"].value_counts().to_string())
    print("\n[SMILES] coverage of the peptide (partners are not the string's job)")
    print(meta["smiles_coverage"].value_counts().to_string())
    print("\n[Termini] what is attached")
    print(meta["partner_class"].value_counts().to_string())

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    base = REPO / "results/half_life_muriverse"

    print(f"\n[Encode] {ESM_MODEL} on {device}, {meta['sequence'].nunique()} sequences")
    esm_features = embed_esm(meta["sequence"].tolist(), device, args.batch_size,
                             cache_path=ROOT / "esm2_embedding_cache.npz")

    drawn = (meta["smiles_coverage"] != "none").to_numpy()
    print(f"\n[Encode] {PEPTIDECLM_MODEL} on {device}, {int(drawn.sum())} peptides "
          f"({meta.loc[drawn, 'species'].value_counts().to_dict()})")
    drawn_features = embed_peptideclm(meta.loc[drawn, "smiles"].tolist(), device,
                                      args.batch_size,
                                      cache_path=ROOT / "peptideclm_embedding_cache.npz")

    # NaN, not zero, on the records with no structure. Zero is a value the model would read
    # as a real position in PeptideCLM space; NaN is passed through StandardScaler (which
    # disregards it when fitting) to XGBoost, which treats it as missing and learns a default
    # split direction instead of being told something false about the molecule.
    clm_features = np.full((len(meta), drawn_features.shape[1]), np.nan, dtype=np.float32)
    clm_features[drawn] = drawn_features

    mol_features, mol_columns = molecule_annotation_block(meta)
    context_features, context_columns = study_context_block(meta)
    ter_features, ter_columns = termini_block(meta)
    # Interpretable physicochemistry of the dosed molecule, read off the verified SMILES.
    # This is the block that lets a configuration distinguish two peptides of the same
    # engineering class, which the bucketed annotation columns cannot; see peptide_physchem.
    phys_features, phys_columns = peptide_physchem.matrix(meta)
    print(f"\n[Encode] {mol_features.shape[1]} molecular annotation + "
          f"{context_features.shape[1]} study context + {ter_features.shape[1]} termini + "
          f"{phys_features.shape[1]} physchem columns")

    available = {
        "sequence": (esm_features, ESM_MODEL, 0, None, None),
        "smiles": (clm_features, PEPTIDECLM_MODEL, int((~drawn).sum()),
                   "NaN; StandardScaler disregards it and XGBoost treats it as a missing "
                   "value", None),
        "molecule_annotation": (mol_features, "one-hot", 0, None, mol_columns),
        "study_context": (context_features, "one-hot + numeric", 0, "NaN for unreported "
                          "dose/age/weight; XGBoost handles missing values", context_columns),
        "termini": (ter_features, "one-hot + log1p mass", 0, None, ter_columns),
        "physchem": (phys_features, "RDKit descriptors + sequence composition",
                     int(np.isnan(phys_features).any(axis=1).sum()),
                     "NaN where the SMILES does not parse; XGBoost treats it as missing",
                     phys_columns),
    }

    def write(out: Path, names: list[str], extra: dict[str, list]) -> DatasetDict:
        blocks, start, described = [], 0, []
        for name in names:
            values, encoder, missing, missing_encoding, named = available[name]
            blocks.append(values)
            entry = {"name": name, "encoder": encoder, "start": start,
                     "dim": int(values.shape[1]), "missing": missing}
            if missing_encoding:
                entry["missing_encoding"] = missing_encoding
            if named:
                entry["columns"] = named
            described.append(entry)
            start += values.shape[1]
        matrix = np.concatenate(blocks, axis=1)
        ds = save_dataset(meta, matrix, out, extra=extra)
        (out / "feature_names.json").write_text(json.dumps(
            {"label": "log10(half_life_h)", "dim": int(matrix.shape[1]),
             "blocks": described}, indent=2) + "\n")
        print(f"[Saved] train={len(ds['train']):4d} val={len(ds['val']):4d}  "
              f"{matrix.shape[1]:5d} features ({' + '.join(names)})\n  {out}")
        return ds

    shared = {"species": meta["species"].tolist(),
              "smiles_coverage": meta["smiles_coverage"].tolist(),
              "cluster": meta["cluster"].tolist(),
              "pmid": meta["PMID"].astype(str).tolist()}

    print()
    suffix = f"{cohort}_strict" if strict else cohort
    write(base / f"halflife_{suffix}_with_embeddings", args.blocks,
          {**shared, "smiles": meta["smiles"].fillna("").tolist()})

    # Every block, whatever the model uses, so evaluate_halflife_model.py can take the
    # matrix apart and show why the shipped combination was chosen rather than asserting it.
    write(base / f"halflife_{suffix}_ablation_with_embeddings",
          ["sequence", "smiles", "molecule_annotation", "study_context", "termini",
           "physchem"],
          {**shared, "smiles": meta["smiles"].fillna("").tolist()})

    if strict and cohort == "rodent" and species == ["mouse", "rat"]:
        write_index(meta_path, cohort_dir / "training_data")


if __name__ == "__main__":
    main()
