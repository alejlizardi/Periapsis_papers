#!/usr/bin/env bash
# Build the arXiv upload bundle for the hexagonal Hardin family paper.
# Produces arxiv_upload.zip containing ONLY what arXiv needs:
#   - hexagonal_hardin_family.tex   (bibliography is inline; no .bib/.bbl)
#   - assets/lockup.png, assets/glyph_square.png   (the two referenced images)
#   - anc/verify_family.py          (ancillary reproducibility script)
# Excludes the compiled PDF (arXiv rebuilds), unused mark.png, and all docs/aux.
set -euo pipefail
cd "$(dirname "$0")"

STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

mkdir -p "$STAGE/assets" "$STAGE/anc"
cp hexagonal_hardin_family.tex "$STAGE/"
cp assets/lockup.png assets/glyph_square.png "$STAGE/assets/"
cp code/verify_family.py "$STAGE/anc/"

rm -f arxiv_upload.zip
( cd "$STAGE" && zip -r -X "$OLDPWD/arxiv_upload.zip" . >/dev/null )

echo "Wrote arxiv_upload.zip:"
unzip -l arxiv_upload.zip
