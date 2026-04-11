#!/usr/bin/env sh
set -eu

if ! command -v codex >/dev/null 2>&1; then
  echo "codex command not found. trying to install Codex CLI..."
  sh ops/02_install_codex_cli.sh
fi

# Prefer current CLI syntax: `codex login`
if codex login --help >/dev/null 2>&1; then
  codex login
elif codex --help >/dev/null 2>&1; then
  echo "This codex version does not support 'codex login'."
  echo "Please run: codex --help"
  exit 1
else
  echo "Unable to detect codex login command."
  exit 1
fi
