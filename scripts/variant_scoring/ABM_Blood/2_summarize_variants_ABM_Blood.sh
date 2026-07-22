#!/bin/bash

set -euo pipefail
shopt -s nullglob

# Keep these values synchronized with 1_score_variants_ABM_Blood.sh.
# variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/neuro.rare.variants.1kg.lt_0.001.tsv
# SET=1kg_rare
variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/Vuckovic_and_Chen_BloodCellGWAS.tsv
SET=Vuckovic_and_Chen_BloodCellGWAS


DATASET=DevBloodMultiome
num_folds=2

root=/data1/offitk/mardera1/chrombpnet_flare
score_dir=$root/output/variant_scores/$SET/$DATASET
merged_dir=$root/output/merged_variant_scores/$SET/$DATASET
summary_dir=$root/output/variant_summary/$SET/$DATASET
log_dir=$root/output/logs/variant_summary/$SET/$DATASET
jobscript=$root/scripts/helper_scripts/summarize_variants.jobscript.sh

partition=cpu
time=60
cpus=1
mem=30G

if [[ ! -d "$score_dir" ]]; then
    echo "ERROR: score directory does not exist: $score_dir" >&2
    exit 1
fi

expected_lines=$(wc -l < "$variants")
cluster_dirs=("$score_dir"/*)

if (( ${#cluster_dirs[@]} == 0 )); then
    echo "ERROR: no cluster directories found in $score_dir" >&2
    exit 1
fi

for cluster_dir in "${cluster_dirs[@]}"; do
    [[ -d "$cluster_dir" ]] || continue
    cluster=$(basename "$cluster_dir")
    complete=true

    echo "$cluster"
    for ((fld = 0; fld < num_folds; fld++)); do
        fold=fold_$fld
        score_files=("$cluster_dir/$fold"/*.variant_scores.tsv)
        observed_lines=0
        if (( ${#score_files[@]} > 0 )); then
            observed_lines=$(awk 'FNR == 1 && $0 ~ /variant_id|allele1/ {next} {count++} END {print count+0}' "${score_files[@]}")
        fi
        echo "  $fold: expected $expected_lines; observed $observed_lines"
        if [[ "$expected_lines" -ne "$observed_lines" ]]; then
            complete=false
        fi
    done

    if [[ "$complete" != true ]]; then
        echo "  SKIPPING: one or more folds are incomplete"
        continue
    fi

    mkdir -p "$merged_dir/$cluster" "$summary_dir" "$log_dir"
    summary_file=$summary_dir/$cluster.mean.variant_scores.tsv
    job_name=$cluster.$DATASET.$SET.summarize

    if [[ -s "$summary_file" ]]; then
        echo "  COMPLETE: $summary_file"
    elif squeue -h -u "$USER" -o "%j" | grep -Fxq "$job_name"; then
        echo "  ALREADY ACTIVE: $job_name"
    else
        echo "  SUBMITTING: $job_name"
        sbatch -J "$job_name" -t "$time" -c "$cpus" --mem="$mem" \
            -p "$partition" --requeue \
            -o "$log_dir/$cluster.summ.log.txt" \
            -e "$log_dir/$cluster.summ.err.txt" \
            "$jobscript" "$score_dir" "$merged_dir" "$summary_dir" \
            "$cluster" "$num_folds"
    fi
done
