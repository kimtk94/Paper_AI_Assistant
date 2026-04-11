#!/usr/bin/env sh
set -eu

if command -v openclaw >/dev/null 2>&1; then
  echo "openclaw already installed: $(openclaw --version 2>/dev/null || true)"
  exit 0
fi

if command -v bash >/dev/null 2>&1; then
  curl -fsSL https://openclaw.ai/install.sh | bash
else
  echo "bash is required by OpenClaw install script."
  exit 1
fi

echo "installed openclaw: $(openclaw --version 2>/dev/null || true)"
