

projPath=/ix/jluke/Projects/sunny/TIIL-017-RoseProj

## exome
project=PDO_sm61 ## PDO_sm31B2 PDO_sm30
metadata=$projPath/sampleinfo/${project}.exome.metadata.20250102.txt
tumornormallist=$projPath/sampleinfo/${project}.exome.tumor_normal.20250102.list
genome=grch38

## only run tumors (as input for facets downstream analysis)
for sample in `awk -F"\t" 'NR>1 {print $1}' $tumornormallist | sort | uniq `; do 
	# echo $sample 

	script=scripts/batches/submit_facets.$sample.sh

	if [ $sample == 'TIIL_0165' ]; then continue; fi 

	echo $script
	sbatch $script
	sleep 2	

done 

chmod u+x scripts/batches/*.sh


