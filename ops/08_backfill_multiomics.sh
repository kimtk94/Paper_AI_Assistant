#!/usr/bin/env sh
set -eu

REPO_DIR=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$REPO_DIR"

INPUT_DIR=${1:-examples/backfill_inputs}
DATA_DIR=${2:-data/multiomics}

if [ ! -d "$INPUT_DIR" ]; then
  echo "input dir not found: $INPUT_DIR"
  exit 1
fi

count=0
for f in "$INPUT_DIR"/*.json; do
  [ -e "$f" ] || continue
  echo "ingesting: $f"
  python3 src/multiomics_ingest.py --input "$f" --data-dir "$DATA_DIR"
  count=$((count + 1))
done

echo "done. ingested files: $count"
