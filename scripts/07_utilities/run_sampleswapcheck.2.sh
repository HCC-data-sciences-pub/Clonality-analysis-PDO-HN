
echo "hostname = $HOSTNAME"
hostname=htc

now=$(date +"%m-%d-%Y_%H:%M:%S")

## -----------------------------------------------

project=PDO_sm61 ## ATM_PAAD MSLNpilot
# sample=$1 ## MTS-1
# genome=$2 ## grch38 grch37 hg19
# part=$3 ## I or II
# aligner=$4 ## bwamem
# caller=$5 ## mutect2
# tumor=$6
# normal=$7

# echo "Sample = $sample; part = $part"
# if [ -z $sample ] || [ -z $genome ] || [ -z $aligner ]; then echo -e "Input missing. \n$0 <sample> <genome>\n"; exit; fi 

caller=somalier ## somalier, crosscheckfingerprints

## -----------------------------------------------

threads=12 ## 12 for star, change mem to 36Gb; 4 for other parts, change mem to 10Gb 

## -----------------------------------------------

projPath=/ix/jluke/Projects/sunny/TIIL-017-RoseProj

## -----------------------------------------------

refGenome=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/GDC/GRCh38.d1.vd1/GRCh38.d1.vd1.fa
refGenomeDict=$refGenome.dict

## -----------------------------------------------

bamfilelist=$projPath/sampleinfo/${project}.exome.sampleswapcheck.20250102.txt
samples=`cut -f 1 $bamfilelist | grep -v sampleid | perl -wpe 's/\s+/ /g'`

## -----------------------------------------------

cd $projPath
pwd=$PWD

echo -e "project = $project\ncaller = $caller\nbamfilelist = $bamfilelist\nPWD = $pwd\n"

## -----------------------------------------------

if [ $caller == 'somalier' ] && [ 1 -eq 1 ]; then 
	
	set -eo pipefail ## important! otherwise script will keep going despite failed

	echo -e "caller = $caller"

	export PATH=$PATH:/ix/rbao/Software/somalier_2/0.2.15

	## -----------------------------------------------

	somalierSites=/ix/rbao/Software/somalier_2/0.2.15/sites_files/sites.hg38.vcf.gz
	if [ $project == 'AMMelanoma' ]; then somalierSites=/ix/rbao/Software/somalier_2/0.2.15/sites_files/sites.hg38.rna.vcf.gz; fi 

	echo -e "somalierSites = $somalierSites\n"

	## -----------------------------------------------

	# sample1=SM21_0565_S2 ## exome
	# sample2=SM21_0566_S1 ## rnaseq

	## -----------------------------------------------

	inString=''
	# for sample in $sample1 $sample2; do 
	for sample in $samples; do

		if [ $sample == 'NA' ] || [ $sample == '' ]; then continue;fi 

		echo $sample 

		## -----------------------------------------------

		outputPath=''
		in=''
		outputPath=`awk -F"\t" -v s=$sample '{ if($1==s) print $3} ' $bamfilelist`
		in=`awk -F"\t" -v s=$sample '{ if($1==s) print $2} ' $bamfilelist`

		echo -e "outputPath = $outputPath\nin = $in"

		## -----------------------------------------------

		tmpDir=$outputPath/tmp
		if [ ! -d $tmpDir ]; then mkdir -p $tmpDir; fi 

		## https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00761-2
		dir='alignment_stats/somalier'
		if [ ! -d $outputPath/$dir ]; then mkdir -p $outputPath/$dir; fi

		## -----------------------------------------------

		## somalier: check sample swaps
		## only need to run once!
		echo "[ " `date` " ] Extracting informative sites for relateness eval: somalier"
		# in=$outputPath/alignment/$sample.$aligner.merged.rmdup.bam
		out=$outputPath/alignment_stats/somalier
		inString="$inString $out/$sample.somalier"

		if [ ! -d $out ]; then mkdir -p $out; fi 

		echo -e "out = $out\n"
		# somalier extract -d $out/ --sites $somalierSites -f $refGenome $in

	done  


	## -----------------------------------------------

	## calculate relateness 
	## can run multiple times with different somalierGroups file!
	echo "[ " `date` " ] relateness calculation: somalier "
	somalierGroups=$projPath/sampleinfo/$project.exome.somalier_group.20250102.txt
	# in=`ls $outputPath/../*/alignment_stats/somalier/*.somalier`
	# out=$outputPath/../$project.somalier.relateness
	out=$outputPath/../$project.exome.somalier.relateness

	echo -e "somalierGroups = $somalierGroups\ninString = $inString\nout = $out\n"
	somalier relate --groups $somalierGroups --min-depth 7 --min-ab 0.3 --infer -o $out $inString

fi 

## -----------------------------------------------

if [ $caller == 'crosscheckfingerprints' ] && [ 0 -eq 1 ]; then 
	
	set -eo pipefail ## important! otherwise script will keep going despite failed

	echo -e "caller = $caller"

	export PATH=$PATH:/ix/rbao/Software/gatk/4.2.6.1

	# module load picard/2.18.12
	PICARD=/ihome/crc/install/picard/2.18.12/picard.jar

	## -----------------------------------------------

	## https://gatk.broadinstitute.org/hc/en-us/articles/360041696232-Detecting-sample-swaps-with-Picard-tools
	## https://docs.google.com/spreadsheets/d/1eFc4y6Km7mhZ5blMwevk86btraA1CUmJ/edit#gid=2008786599
	## https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false
	## https://twitter.com/yfarjoun/status/1324391908253200384
	# haplotypeSNP=/ix/jluke/Projects/sunny/TIIL-001-Augustin-Proj/scripts/broad_haplotypeSNPlist.txt ## this file does not work! missing header
	# haplotypeSNP=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/GDC/broad-references/hg38/v0/hg38_v0_Homo_sapiens_assembly38.haplotype_database.txt ## complained header does not match!!
	haplotypeSNP=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/GDC/broad-references/hg38/v0/hg38_v0_Homo_sapiens_assembly38.haplotype_database.DNAseq.txt ## manually fixed the header in this file to be consistent with header from BAM files

	echo -e "haplotypeSNP = $haplotypeSNP\n"

	## -----------------------------------------------

	# sample1=SM21_0565_S2 ## exome
	# sample2=SM21_0566_S1 ## rnaseq

	## -----------------------------------------------

	inString=''
	# for sample in $sample1 $sample2; do 
	for sample in $samples; do

		if [ $sample == 'NA' ] || [ $sample == '' ]; then continue;fi 

		echo $sample 

		## -----------------------------------------------

		outputPath=''
		in=''
		outputPath=`awk -F"\t" -v s=$sample '{ if($1==s) print $3} ' $bamfilelist`
		in=`awk -F"\t" -v s=$sample '{ if($1==s) print $2} ' $bamfilelist`
		# infn=`echo -e $in | sed 's/\//\t/g' | cut -f 12`
		infn=`awk -F"\t" -v s=$sample '{ if($1==s) print $4} ' $bamfilelist`

		echo -e "outputPath = $outputPath\nin = $in\ninfn = $infn\n"

		## -----------------------------------------------

		tmpDir=$outputPath/tmp
		if [ ! -d $tmpDir ]; then mkdir -p $tmpDir; fi 

		## https://gatk.broadinstitute.org/hc/en-us/articles/360041696232-Detecting-sample-swaps-with-Picard-tools
		dir='alignment_stats/picard'
		if [ ! -d $outputPath/$dir ]; then mkdir -p $outputPath/$dir; fi

		## -----------------------------------------------

		## crosscheckfingerprints: check sample swaps
		## https://gatk.broadinstitute.org/hc/en-us/community/posts/360077548551-Where-do-I-find-ExtractFingerprint
		## only need to run once!
		echo "[ " `date` " ] Extracting informative sites for relateness eval: crosscheckfingerprints"
		# in=$outputPath/alignment/$sample.$aligner.merged.rmdup.bam
		out=$outputPath/alignment_stats/picard/$infn.fingerprint.vcf.gz
		inString="$inString I=$out"

		echo -e "out = $out\n"
		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" ExtractFingerprint -H $haplotypeSNP -I $in -O $out -R $refGenome


	done  

	## -----------------------------------------------

	## calculate relateness 
	## can run multiple times with different lists of input files!
	echo "[ " `date` " ] crosscheck fingerprint: picard "
	out=$projPath/results/$project.exome.picard.crosscheck_metrics

	echo -e "inString = $inString\nout = $out\n"
	java -XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir -jar ${PICARD} CrosscheckFingerprints HAPLOTYPE_MAP=$haplotypeSNP $inString CROSSCHECK_BY=SAMPLE EXPECT_ALL_GROUPS_TO_MATCH=true O=$out MATRIX_OUTPUT=$out.matrix R=$refGenome TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000

fi 

## -----------------------------------------------

echo "[ " `date` " ] Program finished!"