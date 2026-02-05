#!/bin/bash

source ~/.bashrc

micromamba activate chrombpnet2

annotate_script=/data1/offitk/mardera1/github/variant-scorer/src/variant_annotation.py
tss_gencode_coding_file=/data1/offitk/mardera1/data/neuro_variants_notsynapse/hg38.gencode.protein_coding.tss.bed

summary_dir=$1
annotated_dir=$2
peaks=$3
sample=$4

mkdir -p $annotated_dir

out_prefix=$annotated_dir/$sample

python $annotate_script \
    -l $summary_dir/$sample.mean.variant_scores.tsv \
    -o $out_prefix \
    -p $peaks \
    -g $tss_gencode_coding_file \
    -sc chrombpnet

echo "Done Annotating"
echo

# bgzip $out_prefix.annotations.tsv

