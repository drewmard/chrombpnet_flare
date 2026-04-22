
phylop=/data1/offitk/mardera1/data/hg38_refs/hg38.phyloP100way.bw
out=/data1/offitk/mardera1/tmp/phylop_region.bw
bigWigToBedGraph \
  -chrom 6 \
  -start 135322454 \
  -end 135324454 \
  $phylop \
  $out

phylop=/data1/offitk/mardera1/data/hg38_refs/hg38.phyloP100way.bw
out=/data1/offitk/mardera1/chrombpnet_flare/output/variant_shap/ankrd26_5utr/Trisomy_Controls/Megakaryocytes/fold_0/Megakaryocytes.fold_0.phylop_region.bw
bigWigToBedGraph -chrom=chr10 \
  -start=27099325 \
  -end=27101325 \
  "$phylop" \
  "$out"


bigWigToBedGraph -chrom=chr6 -start=135322454 -end=135324454 $phylop $out


chr6:135323327-135323580