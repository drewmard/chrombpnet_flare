#!/bin/bash

# This script processes variant annotation files for different variant sets and datasets.
# It extracts specific columns from each file and saves the output to a designated directory.

# Loop over different variant sets:
variantSet=msk_example
dataset=trevino_2021
for variantSet in msk_example; do

# Extract chr, pos, and gene columns once and store in a file
input_file=/data1/offitk/mardera1/chrombpnet_flare/output/variant_annotations/$variantSet/trevino_2021/trevino_2021.c0.annotations.tsv
metainfofile=/data1/offitk/mardera1/chrombpnet_flare/output/summarize/$variantSet.metainfo.tsv
mkdir -p /data1/offitk/mardera1/chrombpnet_flare/output/summarize
cut -f1-5,32-37 "$input_file" > "$metainfofile"

# Loop over each dataset
# for dataset in ameen_2022 domcke_2020 encode_2024 trevino_2021 corces_2020; do
dataset=domcke_2020
for dataset in trevino_2021 domcke_2020; do
echo $dataset

# Set the 'bias' directory

# Define the input directory path where the variant annotations are stored
dir=/data1/offitk/mardera1/chrombpnet_flare/output/variant_annotations/$variantSet/$dataset

# Define the output directory path for the processed data
outdir=/data1/offitk/mardera1/chrombpnet_flare/output/cbp_sub/$variantSet/$dataset

# Create the output directory if it doesn't already exist
rm -r $outdir
mkdir -p $outdir

# Loop over each file in the input directory
for file in $dir/*; do

# Print the full path of the current file being processed
echo $file

# Extract the filename from the full path (without the directory prefix)
filename=$(basename $file)

# Define the output file path, appending '.sub' to the filename
out=$outdir/${filename}.sub

# Extract cell type identifier
prefix="${dataset}."
suffix=".annotations.tsv"
cell=$(basename "$file" | sed -n "s/${prefix}\(.*\)${suffix}/\1/p")

# Use awk to extract specific columns from the input file 
# Rename them with the cell type suffix
# And with the dataset id
# Then write to the output file:
# BEGIN {OFS="\t"} sets the output field separator to a tab character
# {print $8, $9, $34} extracts columns 8, 9, and 34 from each line of the input file
awk 'BEGIN {OFS="\t"} {print $8, $9, $38}' $file | awk -v cell="$cell" -v dataset="$dataset" 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=3;i++)$i=$i"."cell"."dataset} NR>1{for(i=1;i<=3;i++)$i=$i} {print $1, $2, $3}' > $out

# Print the full path of the output file that is now saved
echo $out

done

# Combine all the columns into a single file
output_file=/data1/offitk/mardera1/chrombpnet_flare/output/summarize/$variantSet.$dataset.txt
paste "$metainfofile" "$outdir"/*.sub > "$output_file"

echo "Output written to $output_file"

done

# Combine all datasets into a single "all_dataset" file
dir=/data1/offitk/mardera1/chrombpnet_flare/output/summarize
# paste $dir/$variantSet.ameen_2022.txt \
#     <(cut -f12- $dir/$variantSet.corces_2020.txt) \
#     <(cut -f12- $dir/$variantSet.domcke_2020.txt) \
#     <(cut -f12- $dir/$variantSet.encode_2024.txt) \
#     <(cut -f12- $dir/$variantSet.trevino_2021.txt) > $dir/$variantSet.all_dataset.txt
paste $dir/$variantSet.domcke_2020.txt \
    <(cut -f12- $dir/$variantSet.trevino_2021.txt) > $dir/$variantSet.all_dataset.txt

echo "Final output merged and written to $dir/$variantSet.all_dataset.txt"

done




