#!/usr/bin/env bash
set -euo pipefail

if [ "${1:-}" = "" ]; then
  echo "Usage: $0 <OPENAI_API_KEY>"
  exit 1
fi

KEY="$1"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_FILE="$REPO_DIR/.env.openclaw"

if [ ! -f "$ENV_FILE" ]; then
  cp "$REPO_DIR/.env.openclaw.example" "$ENV_FILE"
fi

if grep -q '^OPENAI_API_KEY=' "$ENV_FILE"; then
  sed -i "s|^OPENAI_API_KEY=.*|OPENAI_API_KEY=$KEY|" "$ENV_FILE"
else
  echo "OPENAI_API_KEY=$KEY" >> "$ENV_FILE"
fi

echo "OPENAI_API_KEY has been applied to $ENV_FILE"
