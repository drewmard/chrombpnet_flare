#!/bin/bash

set -euo pipefail
shopt -s nullglob

PRECOMPUTEDSCORES=/data1/offitk/mardera1/chrombpnet_flare/output/variant_scores/msk_example/domcke_2020/domcke_2020.fetal_brain.Astrocytes/fold_0/domcke_2020.fetal_brain.Astrocytes.fold_0.msk_example.variant_scores.shuffled.tsv

# variants=/data1/offitk/mardera1/data/neuro_variants_notsynapse/example_scoring_file.tsv
# SET=example
variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/myb_motif_snps.tsv
SET=myb_motif_snps
# variants=/data1/offitk/mardera1/chrombpnet_flare/output/snp_lists/msk_example_list.tsv
# SET=msk_example
DATASET=trevino_2021

variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/myb_motif_snps.tsv
SET=myb_motif_snps
DATASET=trevino_2021

variants=/data1/offitk/mardera1/data/genetics/3333177.AGCTAAGCGG-ATTAATACGC.hard-filtered.gnomad-filtered.andrew_filter1.scoring_file2_chrfilter.txt
SET=ryl1
DATASET=trevino_2021

variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/padhi_predict.tsv
SET=padhi_predict

variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/AoU_eQTL_variants_pip_0.9_pip_0.01.sites.tsv
SET=AoU_eQTL_variants_pip_0.9_pip_0.01

variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/RBC_Trans_Credible_Sets.tsv
SET=RBC_Trans_Credible_Sets

variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/tp53_driver_gwas_list.tsv
SET=tp53_driver_gwas_list
DATASET=Trisomy_Controls

variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/neuro.rare.variants.1kg.lt_0.001.tsv
SET=1kg_rare

variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/Vuckovic_and_Chen_BloodCellGWAS.tsv
SET=Vuckovic_and_Chen_BloodCellGWAS
DATASET=Trisomy_Controls

# HEADDIR=/data1/offitk/mardera1/data/neuro_variants
# for DATASET in domcke_2020 trevino_2021; do
# model_dir=$HEADDIR/chrombpnet_models/$DATASET
# peaks_dir=$HEADDIR/peaks/$DATASET

HEADDIR=/data1/offitk/mardera1/data/Trisomy/Control
model_dir=$HEADDIR/chrombpnet_models/K562_bias
peaks_dir=$HEADDIR/macs_peaks
for DATASET in Trisomy_Controls; do

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

if [[ ! -s "$variants" ]]; then
    echo "ERROR: variant file is missing or empty: $variants" >&2
    exit 1
fi

cluster_dirs=("$model_dir"/*)
if (( ${#cluster_dirs[@]} == 0 )); then
    echo "ERROR: no cluster directories found in $model_dir" >&2
    exit 1
fi

for clust in "${cluster_dirs[@]}"; do
    [[ -d "$clust" ]] || continue
    cluster=$(basename "$clust")

    # skip if it's the SYNAPSE_METADATA_MANIFEST.tsv file
    if [[ "$cluster" == "SYNAPSE_METADATA_MANIFEST.tsv" ]]; then
        continue
    fi

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
        
        # For new datasets:
        peaksFile=$peaks_dir/${cluster}_peaks.narrowPeak.filter_blacklist.bed
        chrombpnetFile=$model_dir/$cluster/$fold/models/chrombpnet_nobias.h5
        
        # # for Marderstein Kundu et al
        # peaksFile=$peaks_dir/$cluster.overlap.peaks.bed.gz
        # chrombpnetFile=$model_dir/$cluster/$fold/chrombpnet_nobias.h5
        
        if [[ ! -s "$peaksFile" ]]; then
            echo "SKIPPING $cluster: missing peak file $peaksFile"
            continue
        fi
        if [[ ! -s "$chrombpnetFile" ]]; then
            echo "SKIPPING $cluster $fold: missing model file $chrombpnetFile"
            continue
        fi

        mkdir -p "$score_dir/$cluster/$fold"
        
        expected_lines=$(wc -l < "$variants")
        score_files=("$score_dir/$cluster/$fold"/*.variant_scores.tsv)
        observed_lines=0
        if (( ${#score_files[@]} > 0 )); then
            observed_lines=$(awk 'FNR == 1 && $0 ~ /variant_id|allele1/ {next} {count++} END {print count+0}' "${score_files[@]}")
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
            sbatch -J "$job_name" -t "$time" -c "$cpus" --mem="$mem" \
                -p "$partitions" --gpus=1 --requeue \
                -o "$log_dir/$cluster/$fold.log.txt" \
                -e "$log_dir/$cluster/$fold.err.txt" \
                "$jobscript" "$scoring_script" \
                    -l "$variants" -g "$hg38_ref_fasta" \
                    -s "$hg38_chrom_sizes" -m "$chrombpnetFile" \
                    -p "$peaksFile" -pg "$hg38_ref_fasta" \
                    -ps "$hg38_chrom_sizes" \
                    -o "$score_dir/$cluster/$fold/$cluster.$fold.$SET" \
                    --shuffled_scores "$PRECOMPUTEDSCORES" \
                    -sc chrombpnet --no_hdf5
        fi
    done
done

done

#                 -t 1000000 \
