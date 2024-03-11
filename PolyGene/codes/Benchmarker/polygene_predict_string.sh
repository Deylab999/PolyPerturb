module load gcc/6.2.0
module load R

IFS="
"                                                                                                                       

for chrom in {1..22}
do
cmd="Rscript polygene_predict_string.R $chrom"
sbatch --time=400:00 --mem=40000 --output=polygene2.out --error=polygene2.err -p short -c 1 --wrap="$cmd"
done


