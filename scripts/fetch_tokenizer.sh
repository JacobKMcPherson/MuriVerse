#!/usr/bin/env bash
# Download the PeptideCLM SMILES-SPE tokenizer files that MuriVerse/build_halflife_training_set.py
# imports. They belong to the upstream PeptiVerse repository, so they are fetched at a pinned
# commit and checked against SHA-256 sums instead of being copied into this repository.
set -euo pipefail

COMMIT=8cf0b21dae356278ae96b414a088e4360357d16c
BASE="https://huggingface.co/ChatterjeeLab/PeptiVerse/resolve/${COMMIT}/tokenizer"
DEST="${PEPTIVERSE_DIR:-$(cd "$(dirname "$0")/.." && pwd)/external/PeptiVerse}/tokenizer"

mkdir -p "$DEST"
for f in my_tokenizers.py new_splits.txt new_vocab.txt; do
  curl -fsSL "$BASE/$f" -o "$DEST/$f"
done

(cd "$DEST" && sha256sum -c - <<'EOF'
753dd37948caaa3ebeceedf8cfc5313477ae3fb1717f522df7a01ae1a39582e9  my_tokenizers.py
90446eadc9a1d8722be58b6b8e31c2756e0d321cd448cb5fae7538a33e3e3af4  new_splits.txt
96af9b4e0aa1fee93a9123720c1a46f977269ecc4ff81ac1cd4e0385b8f4aa2e  new_vocab.txt
EOF
)
echo "tokenizer files are in $DEST"
