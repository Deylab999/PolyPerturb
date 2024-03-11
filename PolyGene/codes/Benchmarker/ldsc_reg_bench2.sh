annot_cell=/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/Benchmarker/PolyGene_Predict_PPI_PoPS_justMAGMA_SAIGE_Mean
baseline_cell=/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/Baselines
baseline_version=baselineLD_v2.1
ldsc_path=/n/groups/price/kushal/LDSC/ldsc
weights_path=/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights
freq_path=/n/groups/price/kushal/LDSC/1000G_Phase3_frq
#sumstats_cell=/n/groups/price/ldsc/sumstats_formatted
sumstats_cell=/n/groups/price/kushal/singlecellLDSC/traits/all_sumstats
output_cell_pre=/n/groups/price/kushal/singlecellLDSC/data/LDSC_RESULTS/Benchmarker/PolyGene_Predict_PPI_PoPS_justMAGMA_SAIGE_Mean

IFS="
"

#sumstats_taskfile=/n/groups/price/kushal/singlecellLDSC/data/traits_bio.txt
#sumstats_taskfile=/n/groups/price/kushal/singlecellLDSC/data/traits_blood.txt  
annot_taskfile=/n/groups/price/kushal/singlecellLDSC/data/PolyGene_Predict_PPI_PoPS_justMAGMA_SAIGE_Mean.txt

module load conda2
source activate ldsc

if [ ! -d $output_cell_pre ]
then
    mkdir $output_cell_pre
fi

output_cell=$output_cell_pre/$baseline_version

if [ ! -d $output_cell ]
then
    mkdir $output_cell
fi

echo $output_cell
for line in `cat $annot_taskfile | awk '{print $1}' | sort | uniq`;
do
    annot_module=`echo $line | awk '{print $1}'`
    echo $annot_cell $annot_module
    if [ ! -d $annot_cell/$annot_module ]
    then
        echo "Error: annotation module directory not found" > ldsc_logfile.log
        exit 100
    fi
    if [ ! -d $output_cell/$annot_module ]
    then
        mkdir $output_cell/$annot_module
    fi
    for ll in `ls $annot_cell/$annot_module | awk '{print $1}' | sort | uniq`;
    do
	annot_dir=`echo $ll | awk '{print $1}'`
	echo $annot_dir
	if [ ! -d $annot_cell/$annot_module/$annot_dir ]
	then
            echo "Error: annotation module directory not found" > ldsc_logfile.log
            exit 101
	fi
	if [ ! -d $output_cell/$annot_module/$annot_dir ]
	then
            mkdir $output_cell/$annot_module/$annot_dir
	fi
        sumstats_file=$annot_module.sumstats
        echo $sumstats_cell $sumstats_file
        if [ ! -f $sumstats_cell/$sumstats_file ]
        then
	    echo "Error: sumstats file not found" > ldsc_logfile.log
	    exit 102
        fi
        if [ ! -f $output_cell/$annot_module/$annot_dir/$sumstats_file.results ]
        then
        cmd="/home/kkd14/.conda/envs/ldsc/bin/python $ldsc_path/ldsc.py  --h2 $sumstats_cell/$sumstats_file --ref-ld-chr $annot_cell/$annot_module/$annot_dir/$annot_dir.,$baseline_cell/$baseline_version/baselineLD.  --frqfile-chr $freq_path/1000G.EUR.QC. --w-ld-chr $weights_path/weights.hm3_noMHC. --overlap-annot --print-coefficients --print-delete-vals --out $output_cell/$annot_module/$annot_dir/$sumstats_file"
         sbatch --time=120:00 --mem=20000 --output=reg_max.out --error=reg_max.err -p short -c 1 --wrap="$cmd"
         fi
	done
done



