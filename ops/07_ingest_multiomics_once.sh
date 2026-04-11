#!/usr/bin/env sh
set -eu

REPO_DIR=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$REPO_DIR"

INPUT_PATH=${1:-examples/multiomics_records.example.json}
DATA_DIR=${2:-data/multiomics}

python3 src/multiomics_ingest.py --input "$INPUT_PATH" --data-dir "$DATA_DIR"
