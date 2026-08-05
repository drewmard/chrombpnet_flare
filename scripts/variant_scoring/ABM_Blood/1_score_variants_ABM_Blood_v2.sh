#!/bin/bash

set -e
set -u
set -o pipefail

PRECOMPUTEDSCORES=/data1/offitk/mardera1/chrombpnet_flare/output/variant_scores/msk_example/domcke_2020/domcke_2020.fetal_brain.Astrocytes/fold_0/domcke_2020.fetal_brain.Astrocytes.fold_0.msk_example.variant_scores.shuffled.tsv

VARIANTS=${1:-${VARIANTS:-/data1/offitk/mardera1/chrombpnet_flare/variant_lists/neuro.rare.variants.1kg.lt_0.001.tsv}}
SET=${2:-${SET:-1kg_rare}}
variants=$VARIANTS
# variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/Vuckovic_and_Chen_BloodCellGWAS.tsv
# SET=Vuckovic_and_Chen_BloodCellGWAS

DATASET=ABM_Blood

ABM_ROOT=/data1/offitk/mardera1/data/cvejic_blood_multiome/abm/chrombpnet_abm
model_dir=$ABM_ROOT/chrombpnet-smk-abm/results/chrombpnet
peaks_dir=$ABM_ROOT/bm_multiome/motif_compendium/results/ABM_all/peaks

echo $DATASET

score_dir=/data1/offitk/mardera1/chrombpnet_flare/output/variant_scores/$SET/$DATASET
log_dir=/data1/offitk/mardera1/chrombpnet_flare/output/logs/variant_scores/$SET/$DATASET
mkdir -p $score_dir
mkdir -p $log_dir

hg38_ref_fasta=/data1/offitk/mardera1/data/neuro_variants_notsynapse/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
hg38_chrom_sizes=/data1/offitk/mardera1/data/neuro_variants_notsynapse/GRCh38_EBV.chrom.sizes.tsv

jobscript=/data1/offitk/mardera1/chrombpnet_flare/scripts/helper_scripts/score_variants.jobscript.sh
scoring_script=/data1/offitk/mardera1/github/variant-scorer/src/variant_scoring.per_chrom.py

# time=2:00:00
# cpus=2
# mem=60G
# partitions=gpushort

time=48:00:00
cpus=2
mem=128G
partitions=gpu

for clust in $model_dir/*; do
    cluster=$(basename $clust)

    echo
    echo $cluster
    echo

    mkdir -p $log_dir/$cluster

    # for fld in {0..4}; do
    for fld in {0..1}; do

        fold=fold_$fld
        echo
        echo $fold
        echo
        
        # Each cell type has fold-specific models and a cell-type peak file.
        peaksFile="$peaks_dir/$cluster/$cluster.no_blacklist.bed"
        chrombpnetFile="$model_dir/$cluster/$fold/models/${cluster}_${fold}._chrombpnet_nobias.h5"

        if [[ ! -s "$peaksFile" ]]; then
            echo "SKIPPING $cluster: missing peak file $peaksFile"
            continue
        fi

        if [[ ! -s "$chrombpnetFile" ]]; then
            echo "SKIPPING $cluster $fold: missing model file $chrombpnetFile"
            continue
        fi
        
        mkdir -p $score_dir/$cluster/$fold
        
        expected_lines=$(wc -l < $variants)
        if [[ -n $(ls -A $score_dir/$cluster/$fold) ]]; then
            observed_lines=$(cat $score_dir/$cluster/$fold/*.variant_scores.tsv | grep -v variant_id | wc -l)
        else
            observed_lines=0
        fi
        
        echo Expected Lines: $expected_lines
        echo Observed Lines: $observed_lines
        
        job_name="${cluster}.${fold}.${SET}"
        
        if [[ "$expected_lines" -eq "$observed_lines" ]]; then
            echo "COMPLETE: $job_name"
        
        elif squeue -h -u "$USER" -o "%j" | grep -Fxq "$job_name"; then
            echo "ALREADY ACTIVE: $job_name"
        
        else
            echo "SUBMITTING: $job_name"
        
            sbatch \
                -J "$job_name" \
                -t "$time" \
                -c "$cpus" \
                --mem="$mem" \
                -p "$partitions" \
                --gpus=1 \
                --requeue \
                -o "$log_dir/$cluster/$fold.log.txt" \
                -e "$log_dir/$cluster/$fold.err.txt" \
                "$jobscript" "$scoring_script" \
                    -l "$variants" \
                    -g "$hg38_ref_fasta" \
                    -s "$hg38_chrom_sizes" \
                    -m "$chrombpnetFile" \
                    -p "$peaksFile" \
                    -pg "$hg38_ref_fasta" \
                    -ps "$hg38_chrom_sizes" \
                    -o "$score_dir/$cluster/$fold/$cluster.$fold.$SET" \
                    --shuffled_scores "$PRECOMPUTEDSCORES" \
                    -sc chrombpnet \
                    --no_hdf5
        fi
    done
done
