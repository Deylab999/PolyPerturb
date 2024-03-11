module load gcc/6.2.0
module load R

IFS="
"                                                                                                                       

for chrom in {1..22}
do
cmd="Rscript polygene_predict_coexpressdb.R $chrom"
sbatch --time=60:00 --mem=20000 --output=polygene.out --error=polygene.err -p short -c 1 --wrap="$cmd"
done


