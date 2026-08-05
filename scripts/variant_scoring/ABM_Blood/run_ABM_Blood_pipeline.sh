#!/bin/bash
#SBATCH --job-name=ABM_Blood.pipeline
#SBATCH --partition=cpu
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/data1/offitk/mardera1/chrombpnet_flare/output/logs/ABM_Blood.pipeline.%j.log
#SBATCH --requeue

# Submit with:
#   sbatch run_ABM_Blood_pipeline.sh
#
# This is a lightweight controller.  The stage scripts submit the real work;
# this job periodically reruns them so that completed prerequisites trigger the
# next stage and incomplete outputs are resubmitted after failed jobs vanish
# from squeue.

set -uo pipefail
shopt -s nullglob

ROOT=/data1/offitk/mardera1/chrombpnet_flare
SCRIPT_DIR=$ROOT/scripts/variant_scoring/ABM_Blood

# These must match the values in the four stage scripts.
DATASET=ABM_Blood
NUM_FOLDS=2
POLL_SECONDS=${POLL_SECONDS:-300}
# SET=1kg_rare
# VARIANTS=$ROOT/variant_lists/neuro.rare.variants.1kg.lt_0.001.tsv
VARIANTS=$ROOT/variant_lists/Vuckovic_and_Chen_BloodCellGWAS.tsv
SET=Vuckovic_and_Chen_BloodCellGWAS

ABM_ROOT=/data1/offitk/mardera1/data/cvejic_blood_multiome/abm/chrombpnet_abm
MODEL_DIR=$ABM_ROOT/chrombpnet-smk-abm/results/chrombpnet
PEAKS_DIR=$ABM_ROOT/bm_multiome/motif_compendium/results/ABM_all/peaks
ANNOTATION_DIR=$ROOT/output/variant_annotations/$SET/$DATASET
SUMMARY_DIR=$ROOT/output/variant_summary/$SET/$DATASET
FINAL_OUTPUT=$ROOT/output/summarize/$SET.$DATASET.txt

stage1=$SCRIPT_DIR/1_score_variants_ABM_Blood.sh
stage2=$SCRIPT_DIR/2_summarize_variants_ABM_Blood.sh
stage3=$SCRIPT_DIR/3_annotate_variants_ABM_Blood.sh
stage4=$SCRIPT_DIR/4_process_scores_ABM_Blood.sh

for script in "$stage1" "$stage2" "$stage3" "$stage4"; do
    if [[ ! -r "$script" ]]; then
        echo "ERROR: cannot read $script" >&2
        exit 1
    fi
done

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

job_active() {
    local job_name=$1 names
    if ! names=$(squeue -h -u "$USER" -o '%j'); then
        return 2
    fi
    grep -Fxq "$job_name" <<< "$names"
}

data_rows() {
    awk '$0 !~ /variant_id|allele1/ {count++} END {print count+0}' "$1"
}

# Failed jobs can leave a nonempty partial file, while the original stage
# submitters use -s as their completion check. Preserve those files for
# debugging and clear the canonical pathname so the stage can be retried.
quarantine_incomplete_outputs() {
    local expected cluster file rows job_name suffix active_status
    expected=$(wc -l < "$VARIANTS")
    suffix=${SLURM_JOB_ID:-manual}.$(date +%s)

    while IFS= read -r cluster; do
        file=$SUMMARY_DIR/$cluster.mean.variant_scores.tsv
        job_name=$cluster.$DATASET.$SET.summarize
        job_active "$job_name"
        active_status=$?
        if [[ -s "$file" && "$active_status" -eq 1 ]]; then
            rows=$(data_rows "$file")
            if [[ "$rows" -ne "$expected" ]]; then
                mv "$file" "$file.incomplete.$suffix"
                echo "[$(timestamp)] Preserved incomplete summary ($rows/$expected rows): $file.incomplete.$suffix"
            fi
        fi

        file=$ANNOTATION_DIR/$cluster.annotations.tsv
        job_name=$cluster.$DATASET.$SET.annotate
        job_active "$job_name"
        active_status=$?
        if [[ -s "$file" && "$active_status" -eq 1 ]]; then
            rows=$(data_rows "$file")
            if [[ "$rows" -ne "$expected" ]]; then
                mv "$file" "$file.incomplete.$suffix"
                echo "[$(timestamp)] Preserved incomplete annotation ($rows/$expected rows): $file.incomplete.$suffix"
            fi
        fi
    done < <(eligible_clusters)
}

# Return the clusters for which all inputs needed by stage 1 exist.  This
# mirrors stage 1's skip rules, so intentionally absent models do not prevent
# the pipeline from finishing.
eligible_clusters() {
    local cluster_dir cluster fold
    for cluster_dir in "$MODEL_DIR"/*; do
        [[ -d "$cluster_dir" ]] || continue
        cluster=$(basename "$cluster_dir")
        [[ -s "$PEAKS_DIR/$cluster/$cluster.no_blacklist.bed" ]] || continue
        for ((fold = 0; fold < NUM_FOLDS; fold++)); do
            [[ -s "$cluster_dir/fold_$fold/models/${cluster}_fold_$fold._chrombpnet_nobias.h5" ]] || break
        done
        (( fold == NUM_FOLDS )) && printf '%s\n' "$cluster"
    done
}

annotations_complete() {
    local cluster found=0
    while IFS= read -r cluster; do
        found=1
        if [[ ! -s "$ANNOTATION_DIR/$cluster.annotations.tsv" ]]; then
            return 1
        fi
    done < <(eligible_clusters)
    (( found == 1 ))
}

echo "[$(timestamp)] Starting $DATASET pipeline controller (poll every ${POLL_SECONDS}s)"

while true; do
    quarantine_incomplete_outputs

    echo "[$(timestamp)] Checking scoring jobs"
    if ! bash "$stage1"; then
        echo "[$(timestamp)] WARNING: stage 1 submitter failed; retrying next pass" >&2
    fi

    echo "[$(timestamp)] Checking summary jobs"
    if ! bash "$stage2"; then
        echo "[$(timestamp)] WARNING: stage 2 submitter failed; retrying next pass" >&2
    fi

    echo "[$(timestamp)] Checking annotation jobs"
    if ! bash "$stage3"; then
        echo "[$(timestamp)] WARNING: stage 3 submitter failed; retrying next pass" >&2
    fi

    if annotations_complete; then
        echo "[$(timestamp)] All annotations complete; running final processing"
        if bash "$stage4" && [[ -s "$FINAL_OUTPUT" ]]; then
            echo "[$(timestamp)] Pipeline complete: $FINAL_OUTPUT"
            exit 0
        fi
        echo "[$(timestamp)] WARNING: final processing failed; retrying next pass" >&2
    else
        echo "[$(timestamp)] Upstream work remains incomplete"
    fi

    sleep "$POLL_SECONDS"
done
