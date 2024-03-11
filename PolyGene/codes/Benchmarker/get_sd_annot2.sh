cellname=/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/Benchmarker/POPs_MAGMA_SAIGE_eqwtd_0kb_qmatched_MAGMA
index=1
module load gcc/6.2.0
module load R

IFS="
"

cmd="Rscript get_sd_annot.R  $cellname $index"
sbatch --time=300:00 --mem=20000 --output=getsd.out --error=getsd.err -p short -c 1 --wrap="$cmd"


