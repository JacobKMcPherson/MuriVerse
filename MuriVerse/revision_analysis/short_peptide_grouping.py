"""Does adding identical-sequence and identical-structure keys change the CV groups of the IV cohort?

Run: python MuriVerse/revision_analysis/short_peptide_grouping.py
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "MuriVerse"))
from evaluate_halflife_model import union_groups  # noqa: E402

meta = pd.read_csv(ROOT / "MuriVerse/rodent/strict/halflife_meta_with_split.csv", low_memory=False)
iv = meta[meta["route_class"].astype(str) == "iv"].reset_index(drop=True)
print("IV records", len(iv))

canon = iv["sequence"].astype(str).str.upper().str.replace(r"[^A-Z]", "", regex=True)
struct = iv["smiles"].astype(str)
short = iv["n_residues"].astype(float) <= 6
print("short (<=6 residues):", int(short.sum()), "ids", iv.loc[short, "id"].tolist())

base = union_groups(iv["cluster"].to_numpy(), iv["PMID"].astype(str).to_numpy())
new = union_groups(iv["cluster"].to_numpy(), iv["PMID"].astype(str).to_numpy(),
                   canon.to_numpy(), struct.to_numpy())

def partition(g):
    return {frozenset(idx) for idx in pd.Series(range(len(g))).groupby(g).apply(list)}

print("groups baseline", len(set(base)), "with identity keys", len(set(new)))
print("partition unchanged:", partition(base) == partition(new))
print("duplicate canonical sequences among short:",
      int(canon[short].duplicated(keep=False).sum()),
      "; duplicate structures among short:", int(struct[short].duplicated(keep=False).sum()))
# Short records that share a canonical sequence or structure with ANY record in a different baseline group.
b = pd.Series(base)
cross = []
for i in iv.index[short]:
    same = iv.index[((canon == canon[i]) | (struct == struct[i])) & (iv.index != i)]
    if any(b[j] != b[i] for j in same):
        cross.append(iv.loc[i, "id"])
print("short records whose identical peptide sits in another baseline group:", cross)
print(iv.loc[short, ["id", "Name", "sequence", "PMID", "species", "cluster"]].to_string())
