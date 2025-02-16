
projPath=/ix/jluke/Projects/sunny/TIIL-017-RoseProj

project=PDO_sm61
aligner=bwamem
caller=gatkscnv
# sampleinfo=$projPath/sampleinfo/$project.genome.tumor_normal.20241223.gatkscnv.tsv

# for datatype in exome genome; do 
for datatype in exome; do 
	echo $datatype

	if [ $datatype == 'exome' ]; then 
		if [ $project == 'PDO_sm30' ]; then 
	        sampleinfo=$projPath/sampleinfo/$project.exome.tumor_normal.20241223.gatkscnv.tsv
	    fi
	    if [ $project == 'PDO_sm61' ]; then 
	        sampleinfo=$projPath/sampleinfo/$project.exome.tumor_normal.20250102.gatkscnv.reorder.tsv
	        ## pooling both batches: sm61 (30+31)
	    fi
	fi 
	if [ $datatype == 'genome' ]; then 
		sampleinfo=$projPath/sampleinfo/$project.genome.tumor_normal.20200526.delly.tsv
	fi 
	
	# ## echo -e "GATK-SCNV detection part I: per sample "
	# ## need to run every sample regardless of tumor or control
	# for sample in `awk -F"\t" '$2=="tumor" || $2=="control" {print $1}' $sampleinfo | sort | uniq`; do

	## echo -e "GATK-SCNV detection part I and II: per tumor vs PoN"
	## only run tumor samples
	for sample in `awk -F"\t" '$2=="tumor" {print $1}' $sampleinfo | sort | uniq`; do 
		
		# echo $sample 

		script=scripts/batches/submit_gatkscnv.$datatype.$sample.sh

		if [ $sample == 'TIIL_0165' ]; then continue; fi 

		echo $script
		sbatch $script
		sleep 2	

	done 

done 

chmod u+x scripts/batches/*.sh


