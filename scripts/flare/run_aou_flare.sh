# Set paths:
repo_dir="/data1/offitk/mardera1/github/FLARE"
script_dir="${repo_dir}/scripts"

SET=AoU_eQTL_variants_pip_0.9_pip_0.01

# ChromBPNet matrix input:
raw_input="/data1/offitk/mardera1/chrombpnet_flare/output/summarize/$SET.Trisomy_Controls.txt"

# Desired output files:
processed_input="${repo_dir}/data/processed/$SET.Trisomy_Controls.FLARE.txt"
model_dir="${repo_dir}/new_models/Trisomy_Controls"
pred_file="${repo_dir}/predictions/$SET.Trisomy_Controls.FLARE.pred.txt"

./FLARE_Preprocess.R \
  -i "$raw_input" \
  -o "$processed_input" \
  -m "$model" \
  --no-training

./FLARE_Predict.R \
  -i "$processed_input" \
  -o "$pred_file" \
  -m "$model_dir"
