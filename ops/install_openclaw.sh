#!/usr/bin/env sh
set -eu

if command -v openclaw >/dev/null 2>&1; then
  echo "openclaw already installed: $(openclaw --version 2>/dev/null || true)"
  exit 0
fi

curl -fsSL https://openclaw.ai/install.sh | sh

echo "installed openclaw: $(openclaw --version 2>/dev/null || true)"
