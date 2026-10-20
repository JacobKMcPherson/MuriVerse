# Data sources

What each input is, where it comes from, and what this repository redistributes. The terms of
use of each source govern reuse of the corresponding files; the CC BY 4.0 license of this
repository applies only to material created for this work.

| source | used for | in this repository |
|---|---|---|
| PEPlife2 (https://webs.iiitd.edu.in/raghava/peplife2/) | peptide half-life records that seed the cohort | the record table `MuriVerse/data/peplife2/peplife2_records.csv` and the API dumps under `MuriVerse/data/peplife2/api/`; rebuilt by `fetch_peplife2_api.py` and `build_peplife2.py`. The mirrored PDFs and site exports are not included |
| PubMed / PMC / Europe PMC (NCBI E-utilities, Europe PMC REST) | study metadata and open-access full text | PubMed metadata XML under `MuriVerse/data/pubmed/`. Retrieved full text is not included; `build_rodent_cohort.py` and `render_fulltext.py` rebuild it |
| Primary papers and patents behind the curation ledgers | dose route and assay method of each record | the ledgers `MuriVerse/data/curation/*.csv`, with the verbatim quote each verdict rests on; the source documents are not included (`fetch_curation_sources.py` re-fetches them) |
| Curated extractions (this work) | per-study half-life, dose, route and matrix values | `MuriVerse/{mouse,rat}/extractions/*.json` |
| Polymyxin B / colistin studies (this work) | held-out application | `MuriVerse/data/polymyxin_*.csv` |
| ThPDB2 and PepTherDia2 | source-database comparison in Supplementary Figure 1 | not included. Per-database species counts derived from them are in `results/database_comparison/`. To regenerate, place `ThPDB2.csv` at `MuriVerse/data/thpdb2/` and `PepTherDia2.csv` at `MuriVerse/data/peptherdia2/` and run `figures/code/02_data_preprocessing.R` |
| PeptiVerse tokenizer (https://huggingface.co/ChatterjeeLab/PeptiVerse) | PeptideCLM SMILES tokenization | not included; downloaded by `scripts/fetch_tokenizer.sh` at a pinned commit |
| ESM-2 (`facebook/esm2_t33_650M_UR50D`), PeptideCLM (`aaronfeller/PeptideCLM-23M-all`) | sequence and SMILES embeddings | not included; downloaded from Hugging Face on first use |
