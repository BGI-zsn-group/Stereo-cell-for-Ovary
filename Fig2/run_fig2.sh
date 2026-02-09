#!/usr/bin/env bash
set -euo pipefail

# Run Fig.2 Harmony integration using YAML config
# From repo root, run:
#   bash Fig2/run_fig2.sh

CONFIG="Fig2/configs/fig2_harmony.yaml"
SCRIPT="Fig2/fig2_harmony_integration.R"

Rscript "$SCRIPT" --config "$CONFIG"
