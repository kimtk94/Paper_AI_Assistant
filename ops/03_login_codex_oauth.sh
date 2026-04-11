#!/usr/bin/env sh
set -eu

if ! command -v codex >/dev/null 2>&1; then
  echo "codex command not found. trying to install Codex CLI..."
  sh ops/02_install_codex_cli.sh
fi

codex --login
