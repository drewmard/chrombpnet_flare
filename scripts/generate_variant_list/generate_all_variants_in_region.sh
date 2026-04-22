fasta=/data1/offitk/mardera1/data/hg38_refs/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
scriptfile=/data1/offitk/mardera1/data/hg38_refs/make_variant_scoring_from_fasta.py
outfile=/data1/offitk/mardera1/chrombpnet_flare/variant_lists/ankrd26_5utr.tsv
region=chr10:27100325-27100501

samtools faidx $fasta $region \
  | python $scriptfile \
  | awk '{for(i=1;i<=10;i++) print}' \
  > $outfile
