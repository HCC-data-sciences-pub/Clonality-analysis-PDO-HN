
projPath=/ix/jluke/Projects/sunny/TIIL-017-RoseProj

## exome
project=PDO_sm61 ## PDO_sm31B2 PDO_sm30
metadata=$projPath/sampleinfo/${project}.exome.metadata.20250102.txt
tumornormallist=$projPath/sampleinfo/${project}.exome.tumor_normal.20250102.list
genome=grch38

for sample in `awk -F"\t" 'NR>1 {print $1}' $tumornormallist | sort | uniq `; do 
	echo $sample 
	tumor=''
	normal=''

	## assuming sample is tumor ...
	tumor=`awk -F"\t" -v s=$sample '$1==s {print $1}' $tumornormallist | sort | uniq | perl -wpe 's/\s+/_/g' | sed 's/_$//' ` ## each row has unique tumor id
	normal=`awk -F"\t" -v s=$sample '$1==s {print $2}' $tumornormallist | sort | uniq | perl -wpe 's/\s+/_/g' | sed 's/_$//' ` ## but same normal might be used in more than one tumor/normal pair
	if [ -z $tumor ]; then 
		## sample is normal!
		tumor=`awk -F"\t" -v s=$sample '$2==s {print $1}' $tumornormallist | sort | uniq | perl -wpe 's/\s+/_/g' | sed 's/_$//'`
		normal=`awk -F"\t" -v s=$sample '$2==s {print $2}' $tumornormallist | sort | uniq | perl -wpe 's/\s+/_/g' | sed 's/_$//'`
	fi 

	if [ $tumor == "" ] || [ -z $tumor ]; then 
		echo -e "cannot match $sample as tumor or normal. program exit"
		exit 
	fi 

	echo -e "tumor = $tumor, normal = $normal"

	script=scripts/submit_facets.sh ## only part II
	scriptsNew=scripts/batches/submit_facets.$sample.sh
	cp $script $scriptsNew
	sed -i "s|^sample=.*|sample=$sample|;s|^genome=.*|genome=$genome|;s|^tumor=.*|tumor=$tumor|;s|^normal=.*|normal=$normal|" $scriptsNew

done 


chmod u+x scripts/batches/*.sh
