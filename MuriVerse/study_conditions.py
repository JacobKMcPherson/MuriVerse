"""Study-condition covariates read from PEPlife2's free-text design columns.

Half-life depends on how it was measured as much as on the molecule: held-out experiments
showed the sample matrix (serum against plasma shifts the median log10 half-life by 0.7), the
assay method and the animal strain separate sources far better than the peptide itself, and
`Dose Route` is blank for half of the records although the sample text often names it. These
are STUDY-DESIGN covariates, in the same sense as `route_class` and `species`: they answer
"what half-life will this molecule show under these conditions", and each is known when a
study is planned.

Every field is a small closed vocabulary with an explicit `unknown` level, so a record whose
text says nothing is not mistaken for one that says "other".
"""

from __future__ import annotations

import pandas as pd

TEXT_COLUMNS = ["Test Sample", "Assay", "Name", "Origin", "Nature", "Activity"]


def _has(series: pd.Series, pattern: str) -> pd.Series:
    import re
    rx = re.compile(pattern)
    return series.map(lambda v: bool(rx.search(v)))


def _first(series: pd.Series, rules: list[tuple[str, str]], default: str) -> pd.Series:
    """First rule whose pattern matches, in order; `default` when none does."""
    out = pd.Series(default, index=series.index, dtype=object)
    done = pd.Series(False, index=series.index)
    for label, pattern in rules:
        hit = _has(series, pattern) & ~done
        out[hit] = label
        done |= hit
    return out


def study_conditions(frame: pd.DataFrame, route_class: pd.Series) -> pd.DataFrame:
    """sample_matrix, assay_method, strain, sex and a text-recovered dose route."""
    sample = frame["Test Sample"].fillna("").astype(str).str.lower()
    assay = frame["Assay"].fillna("").astype(str).str.lower()
    text = frame[TEXT_COLUMNS].fillna("").astype(str).agg(" ".join, axis=1).str.lower()

    matrix = _first(sample, [("plasma", r"plasma"), ("serum", r"serum"),
                             ("blood", r"blood"), ("tissue", r"brain|tissue|liver|kidney|tumou?r")],
                    "unknown")

    # Immunoassays measure the labelled or antigenic species, mass spectrometry the intact
    # peptide, radiometric counting the label, and a bioassay the activity the peptide
    # still has, which is why they are kept apart.
    #
    # A bioassay is its own level rather than part of "unknown". It reports the
    # concentration of ACTIVE peptide inferred from a dose-response - a GLP-2R reporter
    # assay, an ex vivo HIV-1 inhibition titre - so an inactive metabolite that an
    # immunoassay or a mass spectrometer would still count does not register at all. The
    # half-life it yields is an activity half-life measured in plasma, and a model is
    # better off told so than handed it as a missing value.
    method = _first(assay, [("radioimmunoassay", r"radioimmuno|\bria\b"),
                            # Matched on what makes it a bioassay - a cell, a reporter, an
                            # inhibition titre - and NOT on a bare "activity assay", which
                            # also ends "Radioactivity assay": that is a gamma counter.
                            ("bioassay",
                             r"cell[- ]based|cell based|reporter assay|bioassay|"
                             r"receptor[- ]mediated|inhibition assay"),
                            # "EIA" is an enzyme immunoassay and "AlphaLISA" a bead-based
                            # one; a western blot quantifies with an antibody too, so all
                            # three belong with ELISA rather than with "unknown".
                            ("ELISA / immunoassay",
                             r"elisa|immunoassay|\bia\b|\beia\b|alphalisa|lisa\b|"
                             r"western ?blot|immunoblot"),
                            # "LC/MS" and "LCHRMS" are the same instrument as "LC-MS"; the
                            # hyphen-only pattern read them as unknown.
                            ("LC-MS",
                             r"lc[-/ ]?h?r?ms|mass spec|ms/ms|maldi|\blcms\b|"
                             r"liquid chromatography[- ]mass"),
                            ("HPLC", r"hplc|\bu?hplc\b|high.performance liquid"),
                            ("radiometric / imaging",
                             r"radio|gamma|scintill|\bpet\b|spect|fluoresc|imaging")],
                    "unknown")

    strain = _first(sample, [("Sprague-Dawley", r"sprague|\bsd\b"), ("Wistar", r"wistar"),
                             ("C57BL/6", r"c57"), ("BALB/c", r"balb"),
                             ("ICR / Swiss / Kunming / CD-1", r"\bicr\b|kunming|swiss|cd-?1\b")],
                    "unknown")

    sex = _first(sample, [("female", r"female"), ("male", r"\bmale")], "unknown")

    # `Dose Route` is blank on half the records; the sample, assay or name text often names it.
    route = route_class.astype(object).copy()
    blank = route.eq("not reported")
    recovered = _first(text, [("iv", r"intraven|\biv\b|tail vein"), ("sc", r"subcutan|\bsc\b"),
                              ("ip", r"intraperiton|\bip\b"), ("oral", r"\boral|gavage"),
                              ("intranasal", r"intranas")], "")
    use = blank & recovered.ne("")
    route[use] = recovered[use]
    return pd.DataFrame({"sample_matrix": matrix, "assay_method": method, "strain": strain,
                         "sex": sex, "route_model": route})
