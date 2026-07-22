#!/bin/bash

set -euo pipefail
shopt -s nullglob

# Keep these values synchronized with 1_score_variants_ABM_Blood.sh.
variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/neuro.rare.variants.1kg.lt_0.001.tsv
SET=1kg_rare
# variants=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/Vuckovic_and_Chen_BloodCellGWAS.tsv
# SET=Vuckovic_and_Chen_BloodCellGWAS
DATASET=ABM_Blood

root=/data1/offitk/mardera1/chrombpnet_flare
annotation_dir=$root/output/variant_annotations/$SET/$DATASET
outdir=$root/output/cbp_sub/$SET/$DATASET
summary_outdir=$root/output/summarize

annotation_files=("$annotation_dir"/*.annotations.tsv)
if (( ${#annotation_files[@]} == 0 )); then
    echo "ERROR: no annotation files found in $annotation_dir" >&2
    exit 1
fi

mkdir -p "$outdir" "$summary_outdir"

# Variant metadata are identical across cell types, so take them from the first
# annotation file. These columns match the existing stage-4 processing scripts.
metainfo_file=$summary_outdir/$SET.metainfo.tsv
cut -f1-5,32-37 "${annotation_files[0]}" > "$metainfo_file"

sub_files=()
for file in "${annotation_files[@]}"; do
    filename=$(basename "$file")
    cell=${filename%.annotations.tsv}
    out=$outdir/$filename.sub

    awk -v cell="$cell" -v dataset="$DATASET" 'BEGIN {FS=OFS="\t"}
        NR == 1 {print $8 "." cell "." dataset, $9 "." cell "." dataset, $38 "." cell "." dataset; next}
        {print $8, $9, $38}' "$file" > "$out"
    sub_files+=("$out")
    echo "Wrote $out"
done

output_file=$summary_outdir/$SET.$DATASET.txt
paste "$metainfo_file" "${sub_files[@]}" > "$output_file"

echo "Output written to $output_file"
