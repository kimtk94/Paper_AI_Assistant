#!/usr/bin/env sh
set -eu

REPO_DIR=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$REPO_DIR"

echo "[1/3] Pull latest code"
git pull --ff-only

echo "[2/3] Ensure OpenClaw is installed"
if ! command -v openclaw >/dev/null 2>&1; then
  echo "openclaw not found. installing..."
  curl -fsSL https://openclaw.ai/install.sh | sh
fi

echo "[3/3] Restart service"
sudo systemctl restart openclaw
sudo systemctl status openclaw --no-pager -l
