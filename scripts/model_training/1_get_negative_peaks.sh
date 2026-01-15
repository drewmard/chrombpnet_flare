#!/bin/bash

# set -e
# set -u
# set -o pipefail

# srun --nodes=2 --pty /bin/bash
# srun --nodes=1 --cpus-per-task=6 --mem=32G --pty /bin/bash

model_input_dir=/data1/offitk/mardera1/data/Trisomy/Control
filtered_peak_dir=$model_input_dir/macs_peaks
negatives_dir=$model_input_dir/negative_peaks

hg38_exclude=/data1/offitk/mardera1/data/hg38_refs/hg38_exclusion_regions.slop_1057.bed.gz
hg38_ref_fasta=/data1/offitk/mardera1/data/hg38_refs/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
hg38_chrom_sizes=/data1/offitk/mardera1/data/hg38_refs/GRCh38_EBV.chrom.sizes.tsv
hg38_splits_dir=/data1/offitk/mardera1/data/hg38_refs/splits

conda activate chrombpnet2

for celltype in DCs Erythroid Granulocytic_lineage HSCs Immature_B Large_preB Mast_cells_2 Mast_cells Megakaryocytes MEMPs Monocytic_lineage NK_cells pDCs Pre-pro_B_cells ProB_1 ProB_2 T_cells; do
    rm -r $negatives_dir/$celltype
    mkdir -p $negatives_dir/$celltype
    infile=$filtered_peak_dir/${celltype}_peaks.narrowPeak.filter_blacklist.bed

    for fld in {0..4}
    do
        fold=fold_$fld
        chrombpnet prep nonpeaks -g $hg38_ref_fasta \
                                 -p $infile \
                                 -c $hg38_chrom_sizes \
                                 -fl $hg38_splits_dir/$fold.json \
                                 -br $hg38_exclude \
                                 -o $negatives_dir/$celltype/$fold &
    done
    wait

done

