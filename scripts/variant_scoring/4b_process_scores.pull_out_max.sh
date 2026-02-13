in="padhi_predict.Trisomy_Controls.txt"
out="padhi_predict.Trisomy_Controls.max_abs_logfc.tsv"

awk -F'\t' 'BEGIN{OFS="\t"}
NR==1{
  # find abs_logfc.mean.* columns + store their cell type names
  for(i=1;i<=NF;i++){
    if($i ~ /^abs_logfc\.mean\./ && $i !~ /^abs_logfc\.mean\.pval\./){
      idx[++k]=i
      ct[k]=$i
      sub(/^abs_logfc\.mean\./,"",ct[k])
      sub(/\.Trisomy_Controls$/,"",ct[k])
    }
  }
  print "row_id","max_abs_logfc.mean","max_celltype"
  next
}
{
  max=-1; maxct="NA"
  for(j=1;j<=k;j++){
    i=idx[j]
    v=$i+0
    if(v>max){ max=v; maxct=ct[j] }
  }
  print NR-1, max, maxct
}' "$in" > "$out"

cat padhi_predict.Trisomy_Controls.max_abs_logfc.tsv