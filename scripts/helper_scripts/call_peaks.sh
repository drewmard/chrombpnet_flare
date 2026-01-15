# Call peaks

# srun --cpus-per-task=4 --mem=32G --time=06:00:00 --pty bash

conda activate macs
celltype=Mast_cells_2

# cDCs = empty
for celltype in DCs Erythroid Granulocytic_lineage HSCs Immature_B Kupffer_cells Large_preB Low_quality_cells Mast_cells_2 Mast_cells Megakaryocytes MEMPs Monocytic_lineage NK_cells pDCs Pre-pro_B_cells ProB_1 ProB_2 Small_preB T_cells; do

inFile=/data1/offitk/mardera1/data/Trisomy/Control/$celltype.tagAlign.gz
outFile=${celltype}
outDir=/data1/offitk/mardera1/data/Trisomy/Control/macs_peaks

macs3 callpeak -t $inFile -g hs -f BED --nomodel --extsize 200 --shift -100 -q 0.05 -B -n $outFile --outdir $outDir

infile=/data1/offitk/mardera1/data/Trisomy/Control/macs_peaks/${celltype}_peaks.narrowPeak
outfile=/data1/offitk/mardera1/data/Trisomy/Control/macs_peaks/${celltype}_peaks.narrowPeak.filter_blacklist.bed
hg38_exclude=/data1/offitk/mardera1/data/etc/hg38_exclusion_regions.slop_1057.bed.gz
conda activate bedtools
bedtools intersect -v -a $infile -b $hg38_exclude | bedtools sort -i stdin > $outfile
conda deactivate

done

