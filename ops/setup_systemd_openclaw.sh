#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SERVICE_SRC="$REPO_DIR/ops/openclaw.service"
SERVICE_DST="/etc/systemd/system/openclaw.service"

sudo cp "$SERVICE_SRC" "$SERVICE_DST"
sudo systemctl daemon-reload
sudo systemctl enable --now openclaw
sudo systemctl status openclaw --no-pager -l
