#!/bin/bash

set -e
set -u
set -o pipefail

fld=0
variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/myb_motif_snps.tsv
SET=myb_motif_snps
variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/padhi_predict.tsv
SET=padhi_predict
DATASET=Trisomy_Controls
variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/RBC_Trans_Credible_Sets_v2.tsv
SET=RBC_Trans_Credible_Sets_v2
variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/ankrd26_5utr.tsv
SET=ankrd26_5utr
variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/AoU_abbrev.tsv
SET=AoU_abbrev
variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/tp53_driver_gwas_list.tsv
SET=tp53_driver_gwas_list
DATASET=Trisomy_Controls

DATASET=Trisomy_Controls
for DATASET in Trisomy_Controls; do
HEADDIR=/data1/offitk/mardera1/data/Trisomy/Control
model_dir=$HEADDIR/chrombpnet_models/K562_bias
peaks_dir=$HEADDIR/macs_peaks

# variants=/data1/offitk/mardera1/data/genetics/snp_of_interest.txt
# SET=ryl1
# DATASET=trevino_2021
# HEADDIR=/data1/offitk/mardera1/data/neuro_variants
# clust=trevino_2021.c7
# for DATASET in domcke_2020 trevino_2021; do
# model_dir=$HEADDIR/chrombpnet_models/$DATASET
# peaks_dir=$HEADDIR/peaks/$DATASET

score_dir=/data1/offitk/mardera1/chrombpnet_flare/output/variant_scores/$SET/$DATASET
log_dir=/data1/offitk/mardera1/chrombpnet_flare/output/logs/variant_shap/$SET/$DATASET
shap_dir=/data1/offitk/mardera1/chrombpnet_flare/output/variant_shap

hg38_ref_fasta=/data1/offitk/mardera1/data/neuro_variants_notsynapse/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
hg38_chrom_sizes=/data1/offitk/mardera1/data/neuro_variants_notsynapse/GRCh38_EBV.chrom.sizes.tsv

jobscript=/data1/offitk/mardera1/chrombpnet_flare/scripts/helper_scripts/score_variants.jobscript.sh
scoring_script=/data1/offitk/mardera1/github/variant-scorer/src/variant_scoring.py
shap_script=/data1/offitk/mardera1/github/variant-scorer/src/variant_shap.py

time=60
cpus=1
mem=20G
partitions=gpushort

clust=/data1/offitk/mardera1/data/Trisomy/Control/chrombpnet_models/K562_bias/DCs
for clust in $model_dir/*; do
    cluster=$(basename $clust)
    echo
    echo $cluster
    echo
    
    mkdir -p $log_dir/$cluster
    
    # for fld in {1..2}; do
    for fld in 0; do
        fold=fold_$fld
        echo
        echo $fold
        echo
        
        # mkdir -p $shap_dir/$cluster/$SET/$DATASET/$fold
        mkdir -p $shap_dir/$SET/$DATASET/$cluster/$fold/
        # for new:
        chrombpnetFile=$model_dir/$cluster/$fold/models/chrombpnet_nobias.h5
        # # for Marderstein Kundu et al
        # chrombpnetFile=$model_dir/$cluster/$fold/chrombpnet_nobias.h5
        
        #         score_output_file=$shap_dir/$cluster/$fold/$cluster.$fold.common.shap.variant_scores.tsv
        shap_output_file=$shap_dir/$cluster/$fold/$cluster.$fold.shap.variant_shap.counts.h5
        # 
        #         expected_lines=$(wc -l < $variants)
        #         if [[ -n $(ls -A "$shap_dir/$cluster/$fold" 2>/dev/null) && $(ls "$shap_dir/$cluster/$fold"/*.variant_scores.tsv 2>/dev/null) ]]; then
        #             observed_lines=$(cat "$shap_dir/$cluster/$fold"/*.variant_scores.tsv | grep -v variant_id | wc -l)
        #         else
        #             observed_lines=0
        #         fi
        # 
        #         echo Expected Lines: $expected_lines
        #         echo Observed Lines: $observed_lines
        # 
        #         [[ $expected_lines -eq $observed_lines ]] || \
        sbatch -J $cluster.$fold.score -t $time -c $cpus --mem=$mem \
            -p $partitions --gpus 1 --requeue \
            -o $log_dir/$cluster/$fold.variant_scores.log.txt \
            -e $log_dir/$cluster/$fold.variant_scores.err.txt \
            $jobscript $scoring_script \
                -l $variants \
                -g $hg38_ref_fasta \
                -s $hg38_chrom_sizes \
                -m $chrombpnetFile \
                -o $shap_dir/$SET/$DATASET/$cluster/$fold/$cluster.$fold.scores \
                -t 0 \
                -sc chrombpnet
                
        [[ -f $shap_output_file ]] || \
        sbatch -J $cluster.$fold.shap -t $time -c $cpus --mem=$mem \
            -p $partitions --gpus 1 --requeue \
            -o $log_dir/$cluster/$fold.variant_shap.log.txt \
            -e $log_dir/$cluster/$fold.variant_shap.err.txt \
            $jobscript $shap_script \
                -l $variants \
                -g $hg38_ref_fasta \
                -s $hg38_chrom_sizes \
                -m $chrombpnetFile \
                -o $shap_dir/$SET/$DATASET/$cluster/$fold/$cluster.$fold \
                -sc chrombpnet
                # -m $model_dir/$cluster/$fold/models/chrombpnet_nobias.h5 \
                
    done
done
done