#!/usr/bin/env sh
set -eu

REPO_DIR=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$REPO_DIR"

git pull --ff-only
sh ops/01_install_openclaw.sh

echo "Choose auth mode:"
echo "  1) ChatGPT OAuth (no API key)"
echo "  2) OpenAI API key"
printf "Select [1/2]: "
read mode

if [ "$mode" = "1" ]; then
  sh ops/02_install_codex_cli.sh
  sh ops/03_login_codex_oauth.sh
elif [ "$mode" = "2" ]; then
  printf "Enter OPENAI_API_KEY: "
  read key
  sh ops/04_set_openai_key.sh "$key"
else
  echo "Invalid selection"
  exit 1
fi

sh ops/05_setup_systemd_openclaw.sh
