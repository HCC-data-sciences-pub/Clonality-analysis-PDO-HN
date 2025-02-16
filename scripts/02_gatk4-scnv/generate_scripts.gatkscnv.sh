
projPath=/ix/jluke/Projects/sunny/TIIL-017-RoseProj

project=PDO_sm61
aligner=bwamem
caller=gatkscnv

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
		javamem=16g
	fi 
	if [ $datatype == 'genome' ]; then 
		sampleinfo=$projPath/sampleinfo/AMMelanoma.genome.tumor_normal.20200526.delly.tsv
		javamem=32g
	fi 

	for sample in `awk -F"\t" '$2=="tumor" || $2=="control" {print $1}' $sampleinfo | sort | uniq`; do
	# for sample in `awk -F"\t" '$2=="tumor" {print $1}' $sampleinfo | sort | uniq`; do 
		echo $sample 
		
		script=scripts/submit_gatkscnv.sh
		scriptNew=scripts/batches/submit_gatkscnv.$datatype.$sample.sh

		cp  $script $scriptNew
		sed -i "s|^project=.*|project=$project|;s|^sample=.*|sample=$sample|;s|^aligner=.*|aligner=$aligner|;s|^caller=.*|caller=$caller|;s|^datatype=.*|datatype=$datatype|;s|^#SBATCH --mem .*|#SBATCH --mem $javamem|" $scriptNew

	done 

done 

chmod u+x scripts/batches/*.sh


