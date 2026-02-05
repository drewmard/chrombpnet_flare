#!/bin/bash

set -e
set -u
set -o pipefail

variants=/data1/offitk/mardera1/chrombpnet_flare/output/snp_lists/msk_example_list.tsv
SET=msk_example
variants=/data1/offitk/mardera1/data/genetics/3333177.AGCTAAGCGG-ATTAATACGC.hard-filtered.gnomad-filtered.andrew_filter1.scoring_file2_chrfilter.txt
SET=ryl1
DATASET=trevino_2021

for DATASET in domcke_2020 trevino_2021; do
echo $DATASET

score_dir=/data1/offitk/mardera1/chrombpnet_flare/output/variant_scores/$SET/$DATASET
log_dir=/data1/offitk/mardera1/chrombpnet_flare/output/logs/$SET/$DATASET

merged_dir=/data1/offitk/mardera1/chrombpnet_flare/output/merged_variant_scores/$SET/$DATASET
summary_dir=/data1/offitk/mardera1/chrombpnet_flare/output/variant_summary/$SET/$DATASET

processed_dir=/data1/offitk/mardera1/data/neuro_variants/peaks/$DATASET
overlap_peak_dir=$processed_dir

annotated_dir=/data1/offitk/mardera1/chrombpnet_flare/output/variant_annotations/$SET/$DATASET
log_dir=/data1/offitk/mardera1/chrombpnet_flare/output/logs/variant_annotations/$SET/$DATASET

jobscript=/data1/offitk/mardera1/chrombpnet_flare/scripts/helper_scripts/annotate_variants.jobscript.sh

time=120
cpus=1
mem=30G
partitions=cpu # gpushort

for clust in $score_dir/*; do
    cluster=$(basename $clust)
    echo
    echo $cluster
    echo

    mkdir -p $annotated_dir
    mkdir -p $log_dir

    summary_file=$summary_dir/$cluster.mean.variant_scores.tsv

    expected_lines=$(wc -l < $variants)
    if [[ -f $summary_file ]]; then
        observed_lines=$(cat $summary_file | grep -v variant_id | wc -l)
    else
        observed_lines=0
    fi

    echo Expected Lines: $expected_lines
    echo Observed Lines: $observed_lines

    annotated_file=$annotated_dir/$cluster.annotations.tsv

    if [[ $expected_lines -eq $observed_lines ]]; then
        [[ -f $annotated_file ]] || \
        sbatch -J $cluster.$DATASET.$SET -t $time -c $cpus --mem=$mem \
            -p $partitions --requeue \
            -o $log_dir/$cluster.log \
            -e $log_dir/$cluster.err \
            $jobscript \
                $summary_dir \
                $annotated_dir \
                $overlap_peak_dir/$cluster.overlap.peaks.bed.gz \
                $cluster
    fi
done

done
