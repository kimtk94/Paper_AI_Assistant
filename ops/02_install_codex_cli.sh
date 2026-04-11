#!/usr/bin/env sh
set -eu

if command -v codex >/dev/null 2>&1; then
  echo "codex already installed: $(codex --version 2>/dev/null || true)"
  exit 0
fi

if command -v npm >/dev/null 2>&1; then
  npm install -g @openai/codex
elif command -v brew >/dev/null 2>&1; then
  brew install --cask codex
else
  echo "Codex CLI is not installed and no supported package manager was found."
  echo "Install Node.js/npm (recommended) or Homebrew, then rerun this script."
  exit 1
fi

echo "installed codex: $(codex --version 2>/dev/null || true)"
