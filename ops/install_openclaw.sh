#!/usr/bin/env bash
set -euo pipefail

if command -v openclaw >/dev/null 2>&1; then
  echo "openclaw already installed: $(openclaw --version || true)"
  exit 0
fi

curl -fsSL https://openclaw.ai/install.sh | bash

echo "installed openclaw: $(openclaw --version || true)"
