#!/bin/bash

set -e
set -u
set -o pipefail

processed_dir=/data1/offitk/mardera1/data/Trisomy/Control
tagalign_dir=$tagalign_dir

model_input_dir=/data1/offitk/mardera1/data/Trisomy/Control
filtered_peak_dir=$model_input_dir/macs_peaks
negatives_dir=$model_input_dir/negative_peaks

bias_name=K562
bias_model=/data1/offitk/mardera1/data/ENCSR868FGK_bias_fold_0.h5
model_dir=$model_input_dir/${bias_name}_bias
log_dir=/data1/offitk/mardera1/chrombpnet_flare/output/logs/model_training/${bias_name}_bias

hg38_ref_fasta=/data1/offitk/mardera1/data/hg38_refs/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
hg38_chrom_sizes=/data1/offitk/mardera1/data/hg38_refs/GRCh38_EBV.chrom.sizes.tsv
hg38_splits_dir=/data1/offitk/mardera1/data/hg38_refs/splits

jobscript=/data1/offitk/mardera1/chrombpnet_flare/scripts/helper_scripts/train_chrombpnet_models.from_tagaligns.jobscript.sh

time=24:00:00
cpus=2
mem=60G
partitions=gpu

for sample in $negatives_dir/*; do
    sample=$(basename $sample)
    echo $sample

    mkdir -p $log_dir/$sample

    for fld in {0..4}; do
        fold=fold_$fld
        echo $fold

        mkdir -p $model_dir/$sample/$fold

        report=$model_dir/$sample/$fold/evaluation/overall_report.html

        [[ -f $report ]] || \
        sbatch -J $sample.$fold -t $time -c $cpus --mem=$mem \
            -p $partitions --gpus 1 --requeue \
            -o $log_dir/$sample/$fold.log.txt \
            -e $log_dir/$sample/$fold.err.txt \
            $jobscript \
                $hg38_ref_fasta \
                $hg38_chrom_sizes \
                $hg38_splits_dir/$fold.json \
                $filtered_peak_dir/$sample.filtered.peaks.bed.gz \
                $negatives_dir/$sample/${fold}_negatives.bed \
                $tagalign_dir/$sample.tagAlign.gz \
                $bias_model \
                $model_dir/$sample/$fold
    done
done

