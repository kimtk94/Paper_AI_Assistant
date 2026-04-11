#!/usr/bin/env sh
set -eu

if ! command -v codex >/dev/null 2>&1; then
  echo "codex command not found. Please install Codex CLI first."
  exit 1
fi

codex --login
