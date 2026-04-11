#!/usr/bin/env sh
set -eu

REPO_DIR=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
SERVICE_SRC="$REPO_DIR/ops/multiomics_ingest.service"
TIMER_SRC="$REPO_DIR/ops/multiomics_ingest.timer"
SERVICE_DST="/etc/systemd/system/multiomics_ingest.service"
TIMER_DST="/etc/systemd/system/multiomics_ingest.timer"

sudo cp "$SERVICE_SRC" "$SERVICE_DST"
sudo cp "$TIMER_SRC" "$TIMER_DST"
sudo systemctl daemon-reload
sudo systemctl enable --now multiomics_ingest.timer
sudo systemctl list-timers --all | grep multiomics_ingest || true
