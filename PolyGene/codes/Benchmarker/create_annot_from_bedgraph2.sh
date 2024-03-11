module load gcc/6.2.0
module load conda2
source activate ldsc
bedfile_path=/n/groups/price/kushal/singlecellLDSC/data/BEDFILES/Benchmarker/PolyGene_Predict_PPI_PoPS_justMAGMA_Mean
bimfile_path=/n/groups/price/kushal/DATA/BIMS
annot_path=/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/Benchmarker/PolyGene_Predict_PPI_PoPS_justMAGMA_Mean
#ldsc_path=/n/groups/price/kushal/LDSC/ldsc

IFS="
"

#TASKFILE=/n/groups/price/kushal/LDSC-Average/data/Enames.txt
TASKFILE=/n/groups/price/kushal/singlecellLDSC/data/PolyGene_Predict_PPI_PoPS_justMAGMA_Mean.txt
for line in `cat $TASKFILE | awk '{print $1}' | sort | uniq`;
do
    name=`echo $line | awk '{print $1}'`
    if [ ! -d $annot_path/$name ]
    then
	mkdir $annot_path/$name
    fi
    for bedline in `ls $bedfile_path/$name/ | cat | sort | uniq | cut -f 1 -d '.'`;
    do
	bedname=`echo $bedline | awk '{print $1}'`
	if [ ! -d $annot_path/$name/$bedname ]
	then
	    mkdir $annot_path/$name/$bedname
	fi
	if [ ! -f $annot_path/$name/$bedname/$bedname.22.annot.gz ]
	then
	    cmd="~/.conda/envs/ldsc/bin/python  make_annot_combine_from_bedgraph.py --bedname $bedname --bedfile_path $bedfile_path/$name --bimfile_path $bimfile_path --annot_path $annot_path/$name/$bedname"
	    sbatch --time=150:00 --mem=20000 --output=annot.out --error=annot.err -p short -c 1 --wrap="$cmd"
	fi
    done
done
