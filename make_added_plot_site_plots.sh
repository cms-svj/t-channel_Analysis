#!/usr/bin/env bash
set -euo pipefail

analysis_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
site_linear_dir="${analysis_dir}/SVJ-tchannel-Run2_site/linear"
log_dir="${analysis_dir}/logs"
timestamp="$(date +%Y%m%d_%H%M%S)"

cd "${analysis_dir}"
set +u
source condor/initCondor.sh
set -u

mkdir -p \
  "${log_dir}" \
  "${site_linear_dir}/EventLevelDNN/ScoreDistribution" \
  "${site_linear_dir}/Preselection/Nminus1"

python3 -u EventDNNScore_makerusingscores.py \
  --plot-only \
  --draw-raw \
  --raw-only \
  --output EventDNNScore_plots/dnn_scores.root \
  --plot-dir "${site_linear_dir}/EventLevelDNN/ScoreDistribution" \
  |& tee "${log_dir}/site_EventDNNScore_${timestamp}.log"

python3 -u Figure2_makerusingskims.py \
  --plot-only \
  --draw-raw \
  --raw-only \
  --site-layout \
  --site-variable-class non-wnae \
  --years 2016,2017,2018 \
  --output Figure2_plots/noWNAE_selected_MC_vars_Run2.root \
  --plot-dir "${site_linear_dir}" \
  |& tee "${log_dir}/site_Figure2_noWNAE_${timestamp}.log"

python3 -u Figure2_makerusingskims.py \
  --plot-only \
  --draw-raw \
  --raw-only \
  --site-layout \
  --site-variable-class softdrop-only \
  --years 2016,2017,2018 \
  --output Figure2_plots/WNAE_dataMC_vars_Run2.root \
  --plot-dir "${site_linear_dir}" \
  |& tee "${log_dir}/site_Figure2_WNAE_softdrop_${timestamp}.log"

python3 -u Figure2_makerusingskims.py \
  --plot-only \
  --draw-raw \
  --raw-only \
  --site-layout \
  --site-variable-class wnae-only \
  --years 2016,2017,2018 \
  --output Figure2_plots/WNAE_dataMC_vars_Run2.root \
  --plot-dir "${site_linear_dir}" \
  |& tee "${log_dir}/site_Figure2_WNAE_scores_${timestamp}.log"

python3 -u NMinusOne_maker_RA2.py \
  --plot-only \
  --raw-yields \
  --output supplementry_plots_postTrimandgapveto/Nminus1plots/nminus1_RA2_allbkgs_Run2_withSignals.root \
  --plot-dir "${site_linear_dir}/Preselection/Nminus1" \
  |& tee "${log_dir}/site_NMinusOne_RA2_${timestamp}.log"
