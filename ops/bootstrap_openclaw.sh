#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

git pull --ff-only
bash ops/install_openclaw.sh

echo "Choose auth mode:"
echo "  1) ChatGPT OAuth (no API key)"
echo "  2) OpenAI API key"
read -r -p "Select [1/2]: " mode

if [ "$mode" = "1" ]; then
  bash ops/login_codex_oauth.sh
elif [ "$mode" = "2" ]; then
  read -r -p "Enter OPENAI_API_KEY: " key
  bash ops/set_openai_key.sh "$key"
else
  echo "Invalid selection"
  exit 1
fi

bash ops/setup_systemd_openclaw.sh
