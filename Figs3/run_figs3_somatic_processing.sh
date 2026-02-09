#!/usr/bin/env bash
set -euo pipefail

# Figs3: standalone processing workflow (somatic integration)
# From repo root:
#   bash Figs3/run_figs3_somatic_processing.sh
#
CONFIG="Figs3/configs/figs3_somatic_processing.yaml"
Rscript "Figs3/figs3_somatic_processing.R" --config "$CONFIG"
