#!/bin/bash

set -e
set -u
set -o pipefail

# micromamba activate chrombpnet2

variants=/data1/offitk/mardera1/chrombpnet_flare/output/snp_lists/msk_example_list.tsv
SET=msk_example
variants=/data1/offitk/mardera1/data/genetics/3333177.AGCTAAGCGG-ATTAATACGC.hard-filtered.gnomad-filtered.andrew_filter1.scoring_file2_chrfilter.txt
SET=ryl1

DATASET=trevino_2021

variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/padhi_predict.tsv
SET=padhi_predict
DATASET=Trisomy_Controls

# for DATASET in domcke_2020 trevino_2021; do
echo $DATASET

score_dir=/data1/offitk/mardera1/chrombpnet_flare/output/variant_scores/$SET/$DATASET
log_dir=/data1/offitk/mardera1/chrombpnet_flare/output/logs/$SET/$DATASET

merged_dir=/data1/offitk/mardera1/chrombpnet_flare/output/merged_variant_scores/$SET/$DATASET
summary_dir=/data1/offitk/mardera1/chrombpnet_flare/output/variant_summary/$SET/$DATASET

jobscript=/data1/offitk/mardera1/chrombpnet_flare/scripts/helper_scripts/summarize_variants.jobscript.sh

partitions=cpu #gpushort
time=60
cpus=1
mem=30G
# num_folds=5
num_folds=1
latest_fold=fold_$((num_folds - 1))

for clust in $score_dir/*; do
    cluster=$(basename $clust)
    echo
    echo $cluster
    echo

    expected_lines=$(wc -l < $variants)
    if [[ -n $(ls -A $score_dir/$cluster/$latest_fold) ]]; then
        observed_lines=$(cat $score_dir/$cluster/$latest_fold/*.variant_scores.tsv | grep -v variant_id | wc -l)
    else
        observed_lines=0
    fi

    echo Expected Lines: $expected_lines
    echo Observed Lines: $observed_lines

    summary_file=$summary_dir/$cluster.mean.variant_scores.tsv

    if [[ $expected_lines -eq $observed_lines ]]; then

        mkdir -p $merged_dir/$cluster
        mkdir -p $summary_dir
        mkdir -p $log_dir

        [[ -f $summary_file ]] || \
        sbatch -J $cluster.$SET -t $time -c $cpus --mem=$mem \
            -p $partitions --requeue \
            -o $log_dir/$cluster/$fold.summ.log.txt \
            -e $log_dir/$cluster/$fold.summ.err.txt \
            $jobscript \
                $score_dir \
                $merged_dir \
                $summary_dir \
                $cluster \
                $num_folds
    fi
done

done
