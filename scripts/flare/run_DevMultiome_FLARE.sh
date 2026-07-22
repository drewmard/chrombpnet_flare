#!/usr/bin/env bash

# conda activate r

set -euo pipefail

# Set paths:
repo_dir="/data1/offitk/mardera1/github/FLARE"
script_dir="${repo_dir}/scripts"
variants="1kg_rare"
variants="Vuckovic_BloodCellGWAS"

# ChromBPNet matrix input:
raw_input="/data1/offitk/mardera1/chrombpnet_flare/output/summarize/${variants}.DevBloodMultiome.txt"

# Desired output files:
processed_input="${repo_dir}/data/processed/${variants}.DevBloodMultiome.FLARE.txt"
model_dir="${repo_dir}/new_models/DevBloodMultiome"
pred_file="${repo_dir}/predictions/${variants}.DevBloodMultiome.FLARE.pred.txt"
perf_file="${repo_dir}/predictions/${variants}.DevBloodMultiome.FLARE.performance.txt"

##########################################################################################

# Below: preprocessing, training, and evaluation of a FLARE model

# Default params
model="all"
truth_file="${repo_dir}/data/FLARE_training_snps.txt"

mkdir -p "$(dirname "$processed_input")" "$model_dir" "$(dirname "$pred_file")"

cd "$script_dir"

./FLARE_Preprocess.R \
  -i "$raw_input" \
  -o "$processed_input" \
  -m "$model" --no-training

# ./FLARE_Training.R \
#   -i "$processed_input" \
#   -o "$model_dir"

./FLARE_Predict.R \
  -i "$processed_input" \
  -o "$pred_file" \
  -m "$model_dir"

# ./FLARE_Eval/FLARE_lasso_weights.R \
# 	-m "$model_dir" \
# 	-o "$model_dir/lasso_weights.txt" \
# 	-p "$model_dir/lasso_weights_plot.pdf"
# 
./FLARE_Eval/FLARE_performance.R \
  -p "$pred_file" \
  -t "$truth_file" \
  -o "$perf_file"
# 
