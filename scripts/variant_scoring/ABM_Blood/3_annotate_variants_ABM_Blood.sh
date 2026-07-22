#!/bin/bash

set -euo pipefail
shopt -s nullglob

# Keep these values synchronized with 1_score_variants_ABM_Blood.sh.
variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/neuro.rare.variants.1kg.lt_0.001.tsv
SET=1kg_rare
variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/Vuckovic_and_Chen_BloodCellGWAS.tsv
SET=Vuckovic_and_Chen_BloodCellGWAS
DATASET=DevBloodMultiome

root=/data1/offitk/mardera1/chrombpnet_flare
ABM_ROOT=/data1/offitk/mardera1/data/cvejic_blood_multiome/abm/chrombpnet_abm
peaks_dir=$ABM_ROOT/bm_multiome/motif_compendium/results/ABM_all/peaks
score_dir=$root/output/variant_scores/$SET/$DATASET
summary_dir=$root/output/variant_summary/$SET/$DATASET
annotated_dir=$root/output/variant_annotations/$SET/$DATASET
log_dir=$root/output/logs/variant_annotations/$SET/$DATASET
jobscript=$root/scripts/helper_scripts/annotate_variants.jobscript.sh

partition=cpu
time=120
cpus=1
mem=30G

if [[ ! -d "$score_dir" ]]; then
    echo "ERROR: score directory does not exist: $score_dir" >&2
    exit 1
fi

mkdir -p "$annotated_dir" "$log_dir"
expected_lines=$(wc -l < "$variants")
cluster_dirs=("$score_dir"/*)

for cluster_dir in "${cluster_dirs[@]}"; do
    [[ -d "$cluster_dir" ]] || continue
    cluster=$(basename "$cluster_dir")
    summary_file=$summary_dir/$cluster.mean.variant_scores.tsv
    peaks_file=$peaks_dir/$cluster/$cluster.no_blacklist.bed
    annotated_file=$annotated_dir/$cluster.annotations.tsv
    observed_lines=0

    if [[ -s "$summary_file" ]]; then
        observed_lines=$(awk '$0 !~ /variant_id|allele1/ {count++} END {print count+0}' "$summary_file")
    fi
    echo "$cluster: expected $expected_lines; observed $observed_lines"

    if [[ "$expected_lines" -ne "$observed_lines" ]]; then
        echo "  SKIPPING: summary is missing or incomplete"
    elif [[ ! -s "$peaks_file" ]]; then
        echo "  SKIPPING: missing peak file $peaks_file"
    elif [[ -s "$annotated_file" ]]; then
        echo "  COMPLETE: $annotated_file"
    else
        job_name=$cluster.$DATASET.$SET.annotate
        if squeue -h -u "$USER" -o "%j" | grep -Fxq "$job_name"; then
            echo "  ALREADY ACTIVE: $job_name"
        else
            echo "  SUBMITTING: $job_name"
            sbatch -J "$job_name" -t "$time" -c "$cpus" --mem="$mem" \
                -p "$partition" --requeue \
                -o "$log_dir/$cluster.log" -e "$log_dir/$cluster.err" \
                "$jobscript" "$summary_dir" "$annotated_dir" \
                "$peaks_file" "$cluster"
        fi
    fi
done
