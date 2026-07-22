# conda activate flare

mkdir -p /data1/offitk/mardera1/chrombpnet_flare/output/flare/inputs/
mkdir -p /data1/offitk/mardera1/chrombpnet_flare/output/flare/predictions/

variants=msk_example
variants=ryl1
variants=neanderthal_SNPs
cd /data1/offitk/mardera1/github/FLARE/scripts
input="/data1/offitk/mardera1/chrombpnet_flare/output/summarize/$variants.all_dataset.txt"
output="/data1/offitk/mardera1/chrombpnet_flare/output/flare/inputs/$variants.FLARE-fb.txt"
model="fetal_brain"
./FLARE_Preprocess.R -i $input -o $output -m $model

input="/data1/offitk/mardera1/chrombpnet_flare/output/flare/inputs/$variants.FLARE-fb.txt"
output="/data1/offitk/mardera1/chrombpnet_flare/output/flare/predictions/$variants.FLARE-fb.pred.txt"
modelpath="/data1/offitk/mardera1/github/FLARE/models/ASD.FLARE-fb"
./FLARE_Predict.R -i $input -o $output -m $modelpath

