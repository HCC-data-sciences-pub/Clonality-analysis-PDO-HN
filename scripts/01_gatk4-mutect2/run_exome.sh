
echo "hostname = $HOSTNAME"
hostname=htc

now=$(date +"%m-%d-%Y_%H:%M:%S")
alignRawReads=1 ## whether to directly align raw reads or trimmed reads
skipMergeReads=1 ## whether to skip merging 3' overlapping mates
skipsoftclip=0 ## by default use softclipped reads in mutect2 and gatkhc
usepon=1 ## whether to use PON for mutect2 somatic calling
emitgerm=0 ## mutect2 4.0.3.0 or higher : Call all apparent germline site even though they will ultimately be filtered. Required as PureCN input
usegivenvcf=0 ## Only the alleles passed by the user should be considered by mutect2
rmBAM=1 ## remove bam files to save space ...
skipmarkdup=0 ## whether to skip mark dup step, for amplicon data!
mutect2BAM=1 ## whether to generate mutect2 bam file (this slows down mutect2 by a ton)
HaplotypeCallerBAM=1 ## whether to generate gatkhc bam file 

softclipInstring=''
softclipOutstring=''
if [ $skipsoftclip -eq 1 ]; then softclipInstring=' --dont-use-soft-clipped-bases'; softclipOutstring='.nosc';fi 

emitgermInstring=
emitgermOutstring=
if [ $emitgerm -eq 1 ]; then emitgermInstring=' --genotype-germline-sites'; emitgermOutstring='.wgerm';fi 

givenvcfInstring=
givenvcfOutstring=
if [ $usegivenvcf -eq 1 ]; then 
	emitgerm=0

	# userVCF=/gpfs/data/bioinformatics/Projects/CRI-BIO-630-RO-MSpiotto-rbao/results/exome/MSLNfull_samples_grch38/MSLNfull.bwamem.mutect2.pass.decomp.vcf.gz.AF0.1.srt.uniq.vcf.gz

	## I regenerated the united sommut vcf files containing  filtered variants from all 83 tumors ...
	userVCF=/gpfs/data/bioinformatics/Projects/CRI-BIO-630-RO-MSpiotto-rbao/results/exome/MSLNfull_samples_grch38/MSLNfull.tumor83.exome.variants.raw.csv.tumorDPAlt5.flt5.srt.vcf.gz.decomp.norm.uniq.vcf.gz
	echo -e "userVCF = $userVCF"

	## note gatk4.1.3.0 mutect2 does not have the option '--genotyping-mode GENOTYPE_GIVEN_ALLELES ' any more....
	givenvcfInstring=" --alleles $userVCF"
	givenvcfOutstring='.given'
fi 

markdupString='.dedup'
if [ $skipmarkdup -eq 1 ]; then markdupString=''; fi 

## --------------------------------------
## command line options 
## --------------------------------------

project=PDO_sm31B2
sample=$1 ## MTS-1
genome=$2 ## grch38 grch37 hg19
part=$3 ## I or II
aligner=$4 ## bwamem
caller=$5 ## mutect2
tumor=$6
normal=$7
annotator=$8
chr=$9

echo "Sample = $sample; part = $part"
if [ -z $sample ] || [ -z $genome ] || [ -z $aligner ]; then echo -e "Input missing. \n$0 <sample> <genome>\n"; exit; fi 

# genome=$1 ## USA300_FPR3757
# sample=$2 ## mrsaA_F1
# run=$3 ## single multi
# aligner=$4 ## bwamem novo 
# caller=$5 ## mpileup freebayes gatkhc 

## --------------------------------------
## Tools
## --------------------------------------

if [ $hostname == 'htc' ]; then 
	## 9/13/2024: gatk/4.6.0.0 requires higher java version; load java/21.0.2-openjdk resolves this
	modules='gcc/8.2.0 htslib/1.9 bwa/0.7.17 fastqc/0.11.7 bedtools/2.29.0 bedops/2.4.35 bzip2/1.0.6 xz/5.2.3 r/4.1.0 sambamba/0.6.8 samtools/1.9 igv/2.4.16 java/21.0.2-openjdk'

	## annovar/20180416 lancet/1.0.6 abra/2.12 flash/1.2.11 
	## 9/13/2024: update gatk to 4.6.0.0; also update funcotator bundle 
	export PATH=/ix/rbao/Software/UCSCtools:/ix/rbao/Software/star/2.7.11b/Linux_x86_64_static:/ix/rbao/Software/kallisto/0.46.1:/ix/rbao/Software/salmon/1.0.0/bin:/ix/rbao/Software/bamreadcount/0.8.0/bin:/ix/rbao/Software/mosdepth/0.2.6:/ix/rbao/Software/gatk/4.6.0.0:/ix/rbao/Software/annovar/20180416:/ix/rbao/Software/somalier_2/0.2.15:/ix/rbao/Software/flash/2.2.00:/ix/rbao/Software/pigz/2.4:$PATH
	PICARD=/ihome/crc/install/picard/2.18.12/picard.jar
	TRIMMOMATIC=/ihome/crc/install/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar
	HUMANDB=/ix/rbao/Software/annovar/humandb

	## fastqc returned error:
	## /usr/bin/perl: symbol lookup error: /ihome/crc/install/gcc-8.2.0/perl/5.28.0/lib/5.28.0/x86_64-linux-thread-multi/auto/Cwd/Cwd.so: undefined symbol: Perl_xs_handshake
	## solution: https://github.com/bioconda/bioconda-recipes/issues/4390
	export PERL5LIB=""

	## not loaded ...
	# bam-readcount
	# igvtools
	# SeqPrep
	# pigz
	# tabix
	# vcfsorter
	# vcflib/1.0.0 vcftools/0.1.16 rsem/1.3.1 
fi 

# checkQC=/group/bioinformatics/Pipelines/Development/CRI-GMD-Pipeline/gmd-EC2-v1.0.0/BDS-ExScaliburGMD-test/util/Check_QC.pl
# pear=/group/bioinformatics/software/pear/0.9.6/pear-0.9.6-bin-64
# flash=/group/bioinformatics/software/FLASH/1.2.11/flash
# ADAPTERS=/apps/software/java-jdk-1.8.0_92/trimmomatic/0.36/adapters

## --------------------------------------
## Settings
## --------------------------------------

threads=12 ## 12 for star, change mem to 36Gb; 4 for other parts, change mem to 10Gb 
adpFasta="/ix/rbao/Software/bbmap/38.96/resources/adapters.fa" ## TruSeq3-SE.fa NexteraPE-PE.fa
qcMatrics="QUAL:30:50,NFD:0:10,GC:40:60,N:0:2"
trimEncodingParam=" -phred33"
clipEncodingParam=""
pearEncodingParam=" -b 33"
flashEncodingParam=" -p 33"
if [ $adpFasta == 'TruSeq3-PE.fa' ]; then 
	adp1="TACACTCTTTCCCTACACGACGCTCTTCCGATCT" ## TruSeq3-PE.fa R1 
	adp2="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT" ## TruSeq3-PE.fa R2
fi 
if [ $adpFasta == 'TruSeq3-SE.fa' ]; then 
	adp1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  
	adp2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA" 
fi 
if [ $adpFasta == 'NexteraPE-PE.fa' ]; then 
	adp1="AGATGTGTATAAGAGACAG" 
	adp2="AGATGTGTATAAGAGACAG" 
fi 

insertSizeMean=250
insertSizeStd=50

## 10/09/2018: changed dbnsfp33a to dbnsfp35a; clinvar_20170905 to clinvar_20180603
## 10/22/2018: updated db and added more fields....
if [ ! -z "$annotator" ] && [ $annotator == 'annovar' ]; then 
	genomeAnno=hg38
	protocols=( \
		'ensGene' \
		'avsnp150' \
		'cosmic83' \
		'gnomad_exome' \
		'gnomad_genome' \
		'esp6500siv2_all,esp6500siv2_ea,esp6500siv2_aa' \
		'exac03nontcga' \
		'1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas' \
		'kaviar_20150923' \
		'hrcr1' \
		'abraom' \
		'gme' \
		'nci60' \
		'clinvar_20180603' \
		'dbscsnv11' \
		'mcap' \
		'revel' \
		'dbnsfp33a' \
		'dbnsfp31a_interpro' \
		)
	protocol=`(IFS=,; echo "${protocols[*]}")`
	operation='g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f'
fi


## --------------------------------------
## Paths and directories
## --------------------------------------

projPath=/ix/jluke/Projects/sunny/TIIL-017-RoseProj
outputPath="$projPath/results/exome/${project}_samples_${genome}/$sample"
dirlist=" qc_reports clean_reads alignment alignment_stats alignment_stats/picard  alignment_stats/bedtools alignment_stats/mosdepth alignment_stats/somalier somatic_calls somatic_calls/$caller germline_calls germline_calls/gatkhc annotation annotation/$annotator tmp "

if [ $project == 'TCGAUVM' ]; then 
	outputPath="/scratch/rbao/projects/TCGAUVM_tmp/results/exome/${project}_samples_$genome/$sample"
fi 

tmpDir=$outputPath/tmp


## --------------------------------------
## Files
## --------------------------------------

refGenome=
refGenomeDict=
chromSize=
metadata=$projPath/sampleinfo/${project}.exome.metadata.20241231.txt
tumornormalpairs=$projPath/sampleinfo/${project}.exome.tumor_normal.20241231.list
# somalierGroups=$projPath/sampleinfo/$project.exome.somalier_group.txt

if [ $caller == 'mutect2_2' ]; then tumornormalpairs=$tumornormalpairs.2; fi 
if [ $caller == 'mutect2multi' ]; then tumornormalpairs=$tumornormalpairs.multi; fi 

ponInstring=''
ponOutstring=''
ponVCF=$outputPath/../$project.$aligner.mutect2pon${softclipOutstring}.pon.vcf.gz

if [ $project == 'JLUVM' ]; then ponVCF='/scratch/rbao/projects/TCGAUVM_tmp/results/exome/TCGAUVM_samples_grch38/TCGAUVM.bwamem.mutect2pon.nosc.pon.vcf.gz'; fi

if [ $project == 'HCC14' ]; then ponVCF='/ix/rbao/Projects/HCC-CBS-005-Hillman-BLFerris/results/exome/ECOG12_samples_grch38/ECOG12.bwamem.mutect2pon.pon.vcf.gz'; fi 

## not ideal for MLS since not the same sequencing center (CGF vs UGC), baits (V6 vs V6UTR), or technician
## but might be better than none?! MLS does not have matched normal samples, at all....
if [ $project == 'MLS' ] || [ $project == 'AMMelanoma' ] || [ $project == 'MCC' ] | [ $project == 'PDO_sm30' ] | [ $project == 'PDO_sm31B2' ]; then ponVCF='/ix/rbao/Projects/HCC-CBS-005-Hillman-BLFerris/results/exome/ECOG12_samples_grch38/ECOG12.bwamem.mutect2pon.pon.vcf.gz'; fi 

if [ $usepon -eq 1 ]; then ponInstring=" -pon $ponVCF"; ponOutstring='.wpon';fi 

## --------------------------------------
## Commands
## --------------------------------------

for module in $modules; do module load $module; done 
module list 

echo "## --------------------------------------"

cd $projPath
pwd=$PWD
echo "Current Path = $pwd"

for dir in $dirlist; do if [ ! -d $outputPath/$dir ]; then mkdir -p $outputPath/$dir; fi; done 

echo "## --------------------------------------"

echo "Aligner = $aligner"

echo "## --------------------------------------"

## get attributes from the metadata table
rgTotal=`awk -F"\t" -v s=$sample '$1==s' $metadata | cut -f 1- | wc -l` 
echo "readgroup total = $rgTotal"

## =========================== PART I ================================ ## 

## loop through readgroups per sample 
bamString=""
readgroups=""
r1files=""
r2files=""
for (( rgline=1; rgline<=$rgTotal; rgline++ )); do 

	echo "## --------------------------------------"
	echo "readgroup line = $rgline"

	values=( `awk -F"\t" -v s=$sample '$1==s' $metadata | cut -f 1- | head -$rgline | tail -1 ` )
	header=`head -1 $metadata`

	# echo $header
	smIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Sample" {print $1-1}'`
	lbIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Library" {print $1-1}'`
	rgIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="ReadGroup" {print $1-1}'`
	plIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Platform" {print $1-1}'`
	cnIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="SequencingCenter" {print $1-1}'`
	dtIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Date" {print $1-1}'`
	lnIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Lane" {print $1-1}'`
	puIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Unit" {print $1-1}'`
	flIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Flavor" {print $1-1}'`
	ecIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Encoding" {print $1-1}'`
	locIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Location" {print $1-1}'`
	r1Index=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Seqfile1" {print $1-1}'`
	r2Index=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Seqfile2" {print $1-1}'`

	## note shell array index starts with offset 0
	sample_id=${values[smIndex]}
	library=${values[lbIndex]}
	readgroup=${values[rgIndex]}
	platform=${values[plIndex]}
	center=${values[cnIndex]}
	date=${values[dtIndex]}
	lane=${values[lnIndex]}
	unit=${values[puIndex]}
	flavor=${values[flIndex]}
	encoding=${values[ecIndex]}
	loc=${values[locIndex]}
	r1=${values[r1Index]}
	# echo "r2Index = $r2Index"
	if [ ! -z $r2Index ] && [ $r2Index != "" ]; then r2=${values[r2Index]}; fi 
	# echo "r2 = $r2"

	## this should be set for each readgroup NOT over all readgroups
	libraryType='SE'
	if [ ! -z $r2 ] && [ $r2 != "" ];then libraryType='PE'; fi 
	if [ -z $libraryType ] || [ $libraryType == "" ]; then echo "libraryType is not valid. Program terminates."; exit; fi 

	## set adapter sequence accordingly
	if [ $adpFasta == "NexteraPE-PE.fa" ]; then 
		if [ $libraryType != 'PE' ]; then 
			echo "adpFasta is wrong!"
			echo $sample_id $library $readgroup $platform $center $date $lane $unit $flavor $encoding $loc $r1 $r2 $libraryType $adpFasta
			exit
		fi 
	else 
		if [ $libraryType == 'SE' ]; then adpFasta="/ix/rbao/Software/bbmap/38.96/resources/adapters.fa";fi 
		if [ $libraryType == 'PE' ]; then adpFasta="/ix/rbao/Software/bbmap/38.96/resources/adapters.fa";fi 
	fi 

	if [ $adpFasta == 'TruSeq3-PE.fa' ]; then 
		adp1="TACACTCTTTCCCTACACGACGCTCTTCCGATCT" 
		adp2="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT" 
	fi 
	if [ $adpFasta == 'TruSeq3-SE.fa' ]; then 
		adp1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" 
		adp2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA" 
	fi 
	if [ $adpFasta == 'NexteraPE-PE.fa' ]; then 
		adp1="AGATGTGTATAAGAGACAG" 
		adp2="AGATGTGTATAAGAGACAG" 
	fi 

	echo $sample_id $library $readgroup $platform $center $date $lane $unit $flavor $encoding $loc $r1 $r2 $libraryType $adpFasta

	## set reference files 
	if [ $genome == 'grch38' ]; then 
		refGenome=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/GDC/GRCh38.d1.vd1/GRCh38.d1.vd1.fa
		refGenomeDict=$refGenome.dict
		chromSize=$refGenome.chromsize
		bwaIndex=$refGenome
		# GATKbundle=/group/referenceFiles/Homo_sapiens/GDC/GATKbundle/hg38/v0
		broadref=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/GDC/broad-references/hg38/v0
		dbSNP=$broadref/dbsnp_146.hg38.vcf.gz
		G1000snp=$broadref/1000G_phase1.snps.high_confidence.hg38.vcf.gz
		hapmap=$broadref/hapmap_3.3.hg38.vcf.gz
		mills=$broadref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
		omni=$broadref/1000G_omni2.5.hg38.vcf.gz
		axiom=$broadref/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
		knownIndels=$broadref/Homo_sapiens_assembly38.known_indels.vcf.gz
		cosmic=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/cosmic/v90/cosmic_coding_and_noncoding_chr_M_sorted.vcf
		gnomad=$broadref/af-only-gnomad.hg38.vcf.gz
		# ponVCF=
		variantsForContamination=$broadref/small_exac_common_3.hg38.vcf.gz
		## 05/19/2018: changed from S07604624_V6r2_Covered to S07604624_V6UTRr2_Covered
		# targetBED=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/Agilent/S07604624_V6UTRr2_Covered.hg38.merged.srt.bed
		## 07/07/2018: changed from liftover bed to picard bed
		# targetBED=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/Agilent/S07604514_V6r2_Covered.picard.hg38.merged.srt.exRandom_Alt.bed
		## 04/30/2020: changed to UGC target bed 
		# targetBED=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/UGCxgen/ugc-xgen-exome-research-panel-targets-liftOver-to-hg38.bed.exRandom_Alt.bed
		## 11/05/2020: change to ampliseq
		# targetBED=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/Agilent/ampliseq.ComprehensiveCancer.dna_manifest.20180509.picard.hg38.exRandom_Alt.bed
		## 0/13/2024: update to V8 hg38 bed
		targetBED=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/Agilent/S33266340_hg38/S33266340_Covered.exRandom_Alt.bed
		targetInterval=$targetBED.interval_list

		if [ $project == 'NBL' ]; then 
			targetBED=/gpfs/data/bao-lab/projects/20180421-NBL/data/exome/targetBED/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.hg38.bed
			targetInterval=/group/bao-lab/projects/20180421-NBL/data/exome/targetBED/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.hg38.list
		fi 
		if [ $project == 'TCGAUVM' ]; then 
			targetBED=/group/bioinformatics/Projects/CRI-BIO-560-MED-JLuke-rbao/target/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.picard.hg38.exRandom_Alt.bed
			targetInterval=/group/bioinformatics/Projects/CRI-BIO-560-MED-JLuke-rbao/target/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.picard.hg38.exRandom_Alt.bed.interval_list
		fi 
		if [ $project == 'EnRich_Control_Trial' ]; then 
			targetBED=/ix/rbao/Projects/HCC-CBS-100-Hillman-ALee/target/EnRich_Control_Trial.exome.hg38.pad100.srt.merged.bed
			targetInterval=/ix/rbao/Projects/HCC-CBS-100-Hillman-ALee/target/EnRich_Control_Trial.exome.hg38.pad100.srt.merged.bed.interval_list
		fi 

		somalierSites=/ix/rbao/Software/somalier_2/0.2.15/sites_files/sites.hg38.vcf.gz
	fi 

	echo -e "refGenome = $refGenome\nbwaIndex = $bwaIndex\ntargetBED = $targetBED"

	## sample sanity check 
	if [ $sample != $sample_id ]; then 
		echo "Sample is not consistent with metadata table. $sample vs $sample_id"
		exit
	fi 

	## check fastqc format (default assuming it is 33)
	if [ $encoding == "64" ] || [ $encoding == "phred64" ] || [ $encoding == "Phred64" ]; then
		echo "Fastq Encoding = $encoding"
		trimEncodingParam=" -phred64"
		clipEncodingParam=" -6"
		pearEncodingParam=" -b 64"
		flashEncodingParam=" -p 64"

		# ## convert 64 to 33
		# echo "[ " `date` " ] Convert FastQ Phred33 to 64: Trimmomatic"
		# for in in $loc/$r1 $loc/$r2; do 
		# 	echo $in 
		# 	out=$in.TOPHRED33.fq.gz
		# 	java -Xmx12G -jar ${TRIMMOMATIC} SE -threads $threads $trimEncodingParam $in $out TOPHRED33

		# 	mv $in $in.bak
		# 	mv $out $in
		# done 
	fi 

	## readgroup header string (for aligners)
	# RGstring="@RG\tID:$readgroup\tPL:$platform\tLB:$library\tSM:$sample\tCN:$center\tDT:$date\tPU:$unit.$lane"
	## 05/19/2018: do not attach $lane to $unit; $unit already has lane string concatenated
	RGstring="@RG\tID:$readgroup\tPL:$platform\tLB:$library\tSM:$sample\tCN:$center\tDT:$date\tPU:$unit"
	echo "RGstring = $RGstring"

	# RGstringSTAR="ID:$readgroup PL:$platform LB:$library SM:$sample CN:$center DT:$date PU:$unit.$lane"
	RGstringBWA=$RGstring
	
	## =========================== PART II ================================ ## 

	if [ $part == 'I' ]; then 

		## do not set as we want bwa to fail on rc=1
		# set -eo pipefail ## important! otherwise script will keep going despite failed

		echo "## --------------------------------------"

		## fastqc (raw reads)
		echo "[ " `date` " ] Collect raw reads qc metrics: fastqc"
		inString="$loc/$r1"
		if [ $libraryType == 'PE' ]; then inString="$inString $loc/$r2"; fi 
		echo "$inString"
		fastqc --extract -o $outputPath/qc_reports -t $threads --nogroup $inString

		echo "## --------------------------------------"

		## trimmomatic: trim low quality bases, clip adapters
		if [ $alignRawReads -eq 0 ]; then 
			if [ $libraryType == 'SE' ]; then 
				echo "[ " `date` " ] Trim low quality bases (and clip adapters): Trimmomatic"
				in1=$loc/$r1
				out1=$outputPath/clean_reads/$readgroup.R1.trim.fq.gz
				java -Xmx12G -jar ${TRIMMOMATIC} SE -threads $threads $trimEncodingParam $in1 $out1 ILLUMINACLIP:${adpFasta}:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36

				# ## join all SE reads
				# echo "Double check all adapters were removed ... [adapter = $adp1 $adp2]"
				# for file in $out1 ; do 
				# 	echo $file
				# 	zcat $file | grep $adp1 | wc -l
				# 	zcat $file | grep $adp2 | wc -l
				# done 
			fi 
			if [ $libraryType == 'PE' ]; then 
				echo "[ " `date` " ] Trim low quality bases (and clip adapters): Trimmomatic"
				in1=$loc/$r1
				in2=$loc/$r2
				out1=$outputPath/clean_reads/$readgroup.R1.trim.pe.fq.gz
				out2=$outputPath/clean_reads/$readgroup.R1.trim.unpaired.fq.gz
				out3=$outputPath/clean_reads/$readgroup.R2.trim.pe.fq.gz
				out4=$outputPath/clean_reads/$readgroup.R2.trim.unpaired.fq.gz
				out5=$outputPath/clean_reads/$readgroup.R1R2.trim.unpaired.fq.gz
				java -Xmx12G -jar ${TRIMMOMATIC} PE -threads $threads $trimEncodingParam $in1 $in2 $out1 $out2 $out3 $out4 ILLUMINACLIP:${ADAPTERS}/$adpFasta:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36

				## join all SE reads
				tmpFile=$outputPath/clean_reads/fq.tmp
				if [ -e $tmpFile ]; then rm $tmpFile; fi 
				for file in $out2 $out4; do echo $file; zcat $file >> $tmpFile; done 
				gzip -c $outputPath/clean_reads/fq.tmp > $out5
				rm $tmpFile
				# echo "Double check all adapters were removed ... [adapter = $adp1 $adp2]"
				# for file in $out1 $out2 $out3 $out4 $out5; do 
				# 	echo $file
				# 	zcat $file | grep $adp1 | wc -l
				# 	zcat $file | grep $adp2 | wc -l
				# done 
			fi 
		fi 

		echo "## --------------------------------------"

		## flash2: merge 3' overlapping mates
		if [ $alignRawReads -eq 1 ] && [ $skipMergeReads -eq 0 ]; then 
			if [ $libraryType == 'PE' ]; then 
				echo "[ " `date` " ] Merge 3' overlapping mates: flash2"
				in1=$loc/$r1
				in2=$loc/$r2
				outdir=$outputPath/clean_reads
				out=$readgroup.flash
				out1=$out.notCombined_1.fq.gz
				out2=$out.notCombined_2.fq.gz
				out3=$out.extendedFrags.fq.gz
				out4=$out.single.joint.fq.gz
				flash2 -m 10 -M 151 -e 35 --allow-outies --compress-prog=pigz --compress-prog-args='-p8' -t $threads -o $out -d $outdir $in1 $in2 
				rename .fastq.pigz .fq.gz $outdir/$out.*.fastq.pigz

				## join all SE reads
				tmpFile=$outputPath/clean_reads/fq.tmp
				if [ -e $tmpFile ]; then rm $tmpFile; fi 
				for file in $outdir/$out3; do echo $file; zcat $file >> $tmpFile; done 
				gzip -c $outputPath/clean_reads/fq.tmp > $outdir/$out4
				rm $tmpFile
			fi 
		fi 
		if [ $alignRawReads -eq 0 ] && [ $skipMergeReads -eq 0 ]; then 
			if [ $libraryType == 'PE' ]; then 
				echo "[ " `date` " ] Merge 3' overlapping mates: flash2"
				in1=$outputPath/clean_reads/$readgroup.R1.trim.pe.fq.gz
				in2=$outputPath/clean_reads/$readgroup.R2.trim.pe.fq.gz
				in3=$outputPath/clean_reads/$readgroup.R1R2.trim.unpaired.fq.gz
				outdir=$outputPath/clean_reads
				out=$readgroup.trim.flash
				out1=$out.notCombined_1.fq.gz
				out2=$out.notCombined_2.fq.gz
				out3=$out.extendedFrags.fq.gz
				out4=$out.single.joint.fq.gz
				flash2 -m 10 -M 151 -e 35 --allow-outies --compress-prog=pigz --compress-prog-args='-p8' -t $threads -o $out -d $outdir $in1 $in2 
	
				rename .fastq.pigz .fq.gz $outdir/$out.*.fastq.pigz

				## join all SE reads
				tmpFile=$outputPath/clean_reads/fq.tmp
				if [ -e $tmpFile ]; then rm $tmpFile; fi 
				for file in $in3 $outdir/$out3; do echo $file; zcat $file >> $tmpFile; done 
				gzip -c $outputPath/clean_reads/fq.tmp > $outdir/$out4
				rm $tmpFile
			fi 
		fi

		echo "## --------------------------------------"

		## fastqc (clean reads )
		if [ $alignRawReads -eq 0 ] || ( [ $alignRawReads -eq 1 ] && [ $skipMergeReads -eq 0 ]  ); then 

			echo "[ " `date` " ] Collect clean reads qc metrics: fastqc"
			inString=""
			for file in $outputPath/clean_reads/$readgroup.*.gz; do 
				echo $file
				inString="$inString $file"
			done 
			fastqc --extract -o $outputPath/qc_reports -t $threads --nogroup $inString
		fi 

		## =========================== PART III ================================ ## 

		echo "## --------------------------------------"

		if [ $aligner == 'bwamem' ] && [ 1 -eq 1 ]; then 
			## bwa mem
			echo "[ " `date` " ] Map to reference genome: BWA mem"
			if [ $alignRawReads -eq 1 ]; then 
				inString=""
				if [ $libraryType == 'SE' ]; then 
					in1=$loc/$r1
					inString=" $in1"
				fi 
				if [ $libraryType == 'PE' ]; then 
					in1=$loc/$r1
					in2=$loc/$r2
					inString=" $in1 $in2"
					
				fi 
				out1=$outputPath/alignment/$readgroup.$aligner.bam

				echo "Align $libraryType reads .... "
				echo "$inString"
				bwa mem -t $threads -T 30 -R $RGstringBWA $bwaIndex $inString | grep -v "^@PG" | sambamba view -f bam -S -t $threads -o $out1 /dev/stdin
				sambamba sort --tmpdir $tmpDir -o $out1.srt.bam -t $threads $out1

				mv $out1.srt.bam $out1
			fi 
			if [ $alignRawReads -eq 0 ]; then 
				if [ $libraryType == 'SE' ]; then 
					in1=$outputPath/clean_reads/$readgroup.R1.trim.fq.gz
					out1=$outputPath/alignment/$readgroup.$aligner.bam

					echo "Align SE reads .... "
					bwa mem -t $threads -T 30 -R $RGstringBWA $bwaIndex $in1 | grep -v "^@PG" | sambamba view -f bam -S -t $threads -o $out1 /dev/stdin
					sambamba sort --tmpdir $tmpDir -o $out1.srt.bam -t $threads $out1

					mv $out1.srt.bam $out1
				fi 
				if [ $libraryType == 'PE' ]; then 
					in1=$outputPath/clean_reads/$readgroup.trim.flash.notCombined_1.fq.gz
					in2=$outputPath/clean_reads/$readgroup.trim.flash.notCombined_2.fq.gz
					in3=$outputPath/clean_reads/$readgroup.trim.flash.single.joint.fq.gz
					out1=$outputPath/alignment/$readgroup.$aligner.pe.bam
					out2=$outputPath/alignment/$readgroup.$aligner.se.bam
					out3=$outputPath/alignment/$readgroup.$aligner.bam

					echo "Align PE reads .... "
					bwa mem -t $threads -T 30 -R $RGstringBWA $bwaIndex $in1 $in2 | grep -v "^@PG" | sambamba view -f bam -S -t $threads -o $out1 /dev/stdin
					sambamba sort --tmpdir $tmpDir -o $out1.srt.bam -t $threads $out1

					echo "Align SE reads .... "
					bwa mem -t $threads -T 30 -R $RGstringBWA $bwaIndex $in3 | grep -v "^@PG" | sambamba view -f bam -S -t $threads -o $out2 /dev/stdin
					sambamba sort --tmpdir $tmpDir -o $out2.srt.bam -t $threads $out2

					echo "Merge PE and SE read alignment ..."
					sambamba merge -t $threads $out3 $out1.srt.bam $out2.srt.bam
					sambamba index $out3

					echo "Remove intermediate bams..."
					rm $out1 $out1.srt.bam* $out2 $out2.srt.bam*
				fi 
			fi

		fi 

	fi 

	echo "## --------------------------------------"

	## record which bam file output 
	bamString="$bamString $outputPath/alignment/$readgroup.$aligner.bam"

	## record readgroups
	readgroups="$readgroups $readgroup"

	## record R1 fastq files
	r1files="$r1files $loc/$r1"
	r2files="$r2files $loc/$r2"

done

## =========================== PART IV ================================ ## 

if [ $part == 'I' ]; then 

	set -eo pipefail ## important! otherwise script will keep going despite failed

	echo "## --------------------------------------"

	## those steps do not require target bed , run them first 5/22/2022 sunny
	if [  1 -eq 1 ]; then

		## combined readgroup BAM files per sample 
		echo "[ " `date` " ] Merge readgroup alignment into sample level: sambamba"
		echo $bamString
		echo $readgroups
		out=$outputPath/alignment/$sample.$aligner.merged.bam
		if [ $rgTotal -eq 1 ]; then 
			echo "readgroup total = $rgTotal. No need to merge!"
			link=${bamString/*\/}
			if [ ! -e $out ]; then ln -s $link $out; fi 
			if [ ! -e $out.bai ]; then sambamba index -t $threads $out; fi 
		fi 
		if [ $rgTotal -gt 1 ]; then 
			echo "readgroup total = $rgTotal. Merge bamfiles!"
			sambamba merge -t $threads $out $bamString
		fi

		echo "## --------------------------------------"

		if [ $skipmarkdup -eq 0 ]; then 

			## sambamba: mark dups
			echo "[ " `date` " ] Mark duplicates: sambamba "
			in=$outputPath/alignment/$sample.$aligner.merged.bam
			out=$outputPath/alignment/$sample.$aligner.merged.dedup.bam
			sambamba markdup -t $threads --tmpdir=$tmpDir $in $out 

			echo "## --------------------------------------"

			## sambamba: remove dups (for excavator2 copy number detection)
			echo "[ " `date` " ] Remove duplicates: sambamba "
			in=$outputPath/alignment/$sample.$aligner.merged.dedup.bam
			out=$outputPath/alignment/$sample.$aligner.merged.rmdup.bam
			sambamba markdup -t $threads --tmpdir=$tmpDir -r $in $out 
		fi 

		echo "## --------------------------------------"

		## to save space ....
		if [ $rgTotal -gt 1 ] && [ $rmBAM -eq 1 ]; then 
			echo "[ " `date` " ] Removing bam files to save space "
			in=$outputPath/alignment/$sample.$aligner.merged.bam
			echo -e "unlink $in\nrm $in.bai\nrm $bamString"
			unlink $in 
			rm $in.bai
			rm $bamString
		fi 

	fi 
	echo "## --------------------------------------"

	## hmm, might still need it, as some spike in carry mutations (some carry fusions)
	## but first I need to figure out the target bed ..... run these steps second 5/22/2022 sunny
	## at 11 pm - finally I built custom target BED myself by checking previous emails and materials sent by Megan (script build_targetBED.sh)
	if [  1 -eq 1 ]; then 

		## 5/22: I need to quickly check read orientation, manta failed on all samples???
		if [ 0 -eq 1 ]; then 

			echo "## --------------------------------------"

			## picard: collect matrics
			echo "[ " `date` " ] Collect alignment metrics: picard CollectMultipleMetrics"

			echo "Collect multiple metrics ..." ## PROGRAM=RnaSeqMetrics will throw nullpoints in java! 
			## cannot process in mixed PE & SE bam files ... need to run on each
			in=$outputPath/alignment/$sample.$aligner.merged.rmdup.bam
			out=$outputPath/alignment_stats/picard/$sample.$aligner.merged.rmdup.bam
			java -XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir -jar ${PICARD} CollectMultipleMetrics I=$in O=$out R=$refGenome PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000
		fi 
		
		echo "## --------------------------------------"

		## gatk: BQSR
		echo "[ " `date` " ] Base-score quality recalibration: GATK "
		in=$outputPath/alignment/$sample.$aligner.merged${markdupString}.bam
		out1=$outputPath/alignment/$sample.$aligner.merged${markdupString}.recal_table.csv
		out2=$outputPath/alignment/$sample.$aligner.merged${markdupString}.recal.bam

		## if GDC bam files ....
		# in=$loc/${readgroup}_gdc_realn.bam
		# out1=$outputPath/alignment/$sample.$aligner.merged${markdupString}.recal_table.csv
		# out2=$outputPath/alignment/$sample.$aligner.merged${markdupString}.recal.bam

		echo -e "Running BaseRecalibrator ..."
		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" BaseRecalibrator -R $refGenome -I $in -O $out1 -L $targetInterval --known-sites $mills --known-sites $knownIndels --known-sites $dbSNP --use-original-qualities --interval-padding 100

		echo -e "Running ApplyBQSR ..."
		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" ApplyBQSR -R $refGenome -I $in -O $out2 -L $targetInterval -bqsr $out1 --use-original-qualities --interval-padding 100
		sambamba index $out2 

		echo "## --------------------------------------"

		## sambamba: remove dups (for bam-readcount)
		echo "[ " `date` " ] Remove duplicates from recal bam: sambamba "
		in=$outputPath/alignment/$sample.$aligner.merged${markdupString}.recal.bam
		out=$outputPath/alignment/$sample.$aligner.merged${markdupString}.recal.rmdup.bam
		sambamba markdup -t $threads --tmpdir=$tmpDir -r $in $out 

		echo "## --------------------------------------"

		## GATK CollectFragmentCounts (coverage)
		## replaced with CollectReadCounts after GATK 4.0.3.0....
		echo "[ " `date` " ] Collecting fragment counts from bam: GATK"
		in=$outputPath/alignment/$sample.$aligner.merged${markdupString}.recal.bam
		out=$outputPath/alignment/$sample.$aligner.merged${markdupString}.recal.bam.counts

		echo -e "Running CollectReadCounts ..."
		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" CollectReadCounts -I $in -L $targetInterval -R $refGenome --interval-merging-rule OVERLAPPING_ONLY -O $out.hd5 --format HDF5 
		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" CollectReadCounts -I $in -L $targetInterval -R $refGenome --interval-merging-rule OVERLAPPING_ONLY -O $out.tsv --format TSV 

		echo "## --------------------------------------"

		## GATK HaplotypeCaller (germline calls)
		echo "[ " `date` " ] Germline variant calling: HaplotypeCaller"
		in=$outputPath/alignment/$sample.$aligner.merged${markdupString}.recal.bam
		out=$outputPath/germline_calls/gatkhc/$sample.$aligner.gatkhc${softclipOutstring}.g.vcf
		out2=${out/.g.vcf/.bam}

		bamoutString=''
		if [ $HaplotypeCallerBAM -eq 1 ]; then bamoutString=" --bam-output $out2"; fi 

		## do not use -nct ....
		echo -e "skipsoftclip = $skipsoftclip"
		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" HaplotypeCaller -I $in -L $targetInterval -R $refGenome -ERC GVCF --dbsnp $dbSNP -O $out -mbq 20 --interval-padding 100 $softclipInstring $bamoutString

		## Many users have reported issues running HaplotypeCaller with the -nct argument, so we recommend using Queue to parallelize HaplotypeCaller instead of multithreading
		## Note that when HaplotypeCaller is used in GVCF mode (using either -ERC GVCF or -ERC BP_RESOLUTION) the call threshold is automatically set to zero. Call confidence thresholding will then be performed in the subsequent GenotypeGVCFs command. 

		echo "## --------------------------------------"

		## picard: collect matrics
		echo "[ " `date` " ] Collect alignment metrics: picard CollectMultipleMetrics"

		echo "Collect multiple metrics ..." ## PROGRAM=RnaSeqMetrics will throw nullpoints in java! 
		## cannot process in mixed PE & SE bam files ... need to run on each
		in=$outputPath/alignment/$sample.$aligner.merged${markdupString}.recal.bam
		out=$outputPath/alignment_stats/picard/$sample.$aligner.merged${markdupString}.recal.bam
		java -XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir -jar ${PICARD} CollectMultipleMetrics I=$in O=$out R=$refGenome PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000

		echo "Estimate library complexity ..." 
		java -XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir -jar ${PICARD} EstimateLibraryComplexity I=$in O=$out.lib_complex_metrics TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000 


		# ## 03/11/2019
		# ## for large BAM may need large mem
		# echo "Estimate Hs metrics ..." 
		# for in in $outputPath/alignment/$sample.$aligner.merged.bam $outputPath/alignment/$sample.$aligner.merged${markdupString}.bam $outputPath/alignment/$sample.$aligner.merged${markdupString}.recal.bam; do 
			
		# 	echo $in
		# 	out=${in/alignment\//alignment_stats/picard\/}
		# 	echo $out

		# 	java -XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir -jar ${PICARD} CollectHsMetrics I=$in O=$out.hs_metrics PER_TARGET_COVERAGE=$out.per_target_cov R=$refGenome BAIT_INTERVALS=$targetInterval TARGET_INTERVALS=$targetInterval TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000 
		# done 

		# echo "## --------------------------------------"

		# ## output MapQ distribution to make sure featureCounts behave as expected 
		# ## (eg. counts only primary & unique alignment)
		# echo "[ " `date` " ] Count MapQ distribution of the bam file ..."
		# in=$outputPath/alignment/$sample.$aligner.merged${markdupString}.recal.bam
		# sambamba view -t $threads $in | cut -f 5 | sort | uniq -c > $in.mapq_hist

	fi 

	echo "## --------------------------------------"

	## those steps require target bed , run them second 5/22/2022 sunny
	if [ 1 -eq 1 ]; then

		## mosdepth: coverage
		echo "[ " `date` " ] Calculating read coverage: mosdepth, with dups"
		in=$outputPath/alignment/$sample.$aligner.merged${markdupString}.bam
		out=$outputPath/alignment_stats/mosdepth/$sample.$aligner.merged${markdupString}.bam

		## https://broadinstitute.github.io/picard/explain-flags.html
		## sam flag 1796 = read unmapped (0x4) + not primary alignment (0x100) + read fails platform/vendor quality checks (0x200) + read is PCR or optical duplicate (0x400)
		## developer: use no more than 4 threads 
		## include dup reads
		mosdepth -t 4 -b $targetBED -n -F 772 -T 1,5,10,20,30 $out.wdup $in
		## exclude dup reads
		mosdepth -t 4 -b $targetBED -n -F 1796 -T 1,5,10,20,30 $out.wodup $in

		# echo "## --------------------------------------"

		# # decided to swap bedtools out with mosdepth; bedtools takes too much
		# # mem and not multithreaded!

		# bedtools coverage -b $in -a $targetBED -hist | grep ^all > $out

		# scaled=50000000
		# scaleFactor=`awk -v s=$scaled '$2==255 {print s/$1}' $outputPath/alignment/$sample.$aligner.merged.bam.mapq_hist`
		# echo -e "scaleFactor = $scaleFactor"

		# echo "Calculating genome coverage ... "
		# bedtools genomecov -bga -split -ibam $in -scale $scaleFactor > $out.bdg

		# echo "Converting bedgraph to bigwig ... "
		# bedGraphToBigWig $out.bdg $chromSize $out.bw

		# echo "## --------------------------------------"

		# ## somalier: check sample swaps
		# echo "[ " `date` " ] Extracting informative sites for relateness eval: somalier"
		# in=$outputPath/alignment/$sample.$aligner.merged${markdupString}.recal.bam
		# out=$outputPath/alignment_stats/somalier

		# somalier extract -d $out/ --sites $somalierSites -f $refGenome $in
	fi 

fi 

# =========================== PART V ================================ ## 

if [ $part == 'II' ]; then 

	echo -e "tumor = $tumor; normal = $normal"

	set -eo pipefail ## important! otherwise script will keep going despite failed

	echo "## --------------------------------------"

	if [ $caller == 'mutect2' ] || [ $caller == 'mutect2_2' ] || [ $caller == 'mutect2multi' ] || [ $caller == 'mutect2pon' ] || 
		[ $caller == 'lancet' ]; then 

		echo "caller = $caller. No need to run bam cocleaning!"
	else 

		## mrsa: indel cocleaning 
		echo "[ " `date` " ] indel realignment + cocleaning: abra "
		tumorBAM=$outputPath/alignment/$tumor.$aligner.merged${markdupString}.recal
		normalBAM=$outputPath/../$normal/alignment/$normal.$aligner.merged${markdupString}.recal
		knownIndelsABRA=${knownIndels/.gz}

		echo -e "tumorBAM = $tumorBAM.bam; normalBAM = $normalBAM.bam"
		java -Xmx12G -jar ${ABRA} --in $normalBAM.bam,$tumorBAM.bam --out $normalBAM.coclean.bam,$tumorBAM.coclean.bam --ref $refGenome --threads $threads --targets $targetBED --dist 1000 --tmpdir $tmpDir --in-vcf $knownIndelsABRA --log warn
		sambamba index $normalBAM.coclean.bam
		sambamba index $tumorBAM.coclean.bam
	fi 

	echo "## --------------------------------------"

	## mutect2: somatic calling 
	if [ $caller == 'mutect2' ] || [ $caller == 'mutect2_2' ] ; then 

		echo "[ " `date` " ] Somatic variant calling: $caller "
		tumorBAM=$outputPath/alignment/$tumor.$aligner.merged${markdupString}.recal.bam
		normalBAM=$outputPath/../$normal/alignment/$normal.$aligner.merged${markdupString}.recal.bam
		out=$outputPath/somatic_calls/$caller/$tumor.$aligner.$caller${softclipOutstring}${ponOutstring}${emitgermOutstring}${givenvcfOutstring}.vcf.gz
		out2=$outputPath/somatic_calls/$caller/$tumor.$aligner.$caller${softclipOutstring}${ponOutstring}${emitgermOutstring}${givenvcfOutstring}.f1r2.tar.gz
		out3=${out/.vcf.gz/.bam}

		bamoutString=''
		if [ $mutect2BAM -eq 1 ]; then bamoutString=" --bam-output $out3"; fi 

		gatk GetSampleName -I $tumorBAM -O $tumorBAM.SampleName.txt
		tumorName=`awk '{print $1}' $tumorBAM.SampleName.txt | perl -wpe 's/\s+//g'`

		normalInstring=''
		normalName=''
		if [ $normal == 'NA' ] || [ $normal == '' ] || [ $normal == '-' ]; then normalInstring=''; 
		else 
			gatk GetSampleName -I $normalBAM -O $normalBAM.SampleName.txt
			normalName=`awk '{print $1}' $normalBAM.SampleName.txt | perl -wpe 's/\s+//g'`

			normalInstring=" -I $normalBAM -normal $normalName" ## fix 10/09/2018
		fi 

		echo "tumorName = $tumorName; normalName = $normalName"

		echo -e "skipsoftclip = $skipsoftclip"
		echo -e "emitgerm = $emitgerm"
		# gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" Mutect2 -R $refGenome -I $tumorBAM -tumor $tumorName $normalInstring -O $out -L $targetInterval --germline-resource $gnomad -mbq 20 -VS LENIENT --interval-padding 100 $softclipInstring $ponInstring $emitgermInstring $givenvcfInstring

		## 11/16/2019: error: Cannot construct fragment from more than two reads
		## temp solution .... https://gatkforums.broadinstitute.org/gatk/discussion/24492/mutect2-java-lang-illegalargumentexception-cannot-construct-fragment-from-more-than-two-reads#latest
		## 5/4/2020: seems this issue was fixed in gatk 4.1.5.0
		## 5/4/2020: remove --independent-mates in this run (gatk 4.1.7.0)
		## 5/5/2020: new strategy to filter mutect2 calls!
		## https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
		## need to run Mutect2 with the --f1r2-tar-gz argument in order to use the read orientation fitler.... have to rerun mutect2 jobs!!
		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" Mutect2 -R $refGenome -I $tumorBAM -tumor $tumorName $normalInstring -O $out -L $targetInterval --germline-resource $gnomad -mbq 20 -VS LENIENT --interval-padding 100 --f1r2-tar-gz $out2 $softclipInstring $ponInstring $emitgermInstring $givenvcfInstring $bamoutString

		echo "## --------------------------------------"

		## mutect2: filter variants 
		## 5/5/2020: new strategy to filter mutect2 calls!
		## https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
		## need to run Mutect2 with the --f1r2-tar-gz argument in order to use the read orientation fitler.... have to rerun mutect2 jobs!!
		echo "[ " `date` " ] Somatic variant filtering: $caller "
		tumorBAM=$outputPath/alignment/$tumor.$aligner.merged${markdupString}.recal.bam
		normalBAM=$outputPath/../$normal/alignment/$normal.$aligner.merged${markdupString}.recal.bam
		fnVCF=$outputPath/somatic_calls/$caller/$tumor.$aligner.$caller${softclipOutstring}${ponOutstring}${emitgermOutstring}${givenvcfOutstring}
		in=$fnVCF.vcf.gz
		in2=$fnVCF.bam
		out1=$fnVCF.f1r2.tar.gz
		out2=$fnVCF.read-orientation-model.tar.gz
		out3=$fnVCF.pileups.table
		out3Normal=$outputPath/somatic_calls/$caller/$normal.$aligner.$caller${softclipOutstring}${ponOutstring}${emitgermOutstring}${givenvcfOutstring}.pileups.table
		out4=$fnVCF.segments.table
		out5=$fnVCF.contamination.table
		out6=$fnVCF.flt.vcf.gz
		out7=$fnVCF.flt.flt.vcf.gz

		echo -e "Passing f1r2.tar.gz from MuTect2 run to generate LearnReadOrientationModel...\n+++++++++++++++++"
 		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" LearnReadOrientationModel -I $out1 -O $out2

 		nomalContamString=''
 		nomalContamString="--tumor-segmentation $out4 --contamination-table $out5"
		if [ $normal == 'NA' ] || [ $normal == '' ] || [ $normal == '-' ]; then 
			echo -e "normal is $normal"

			nomalContamString="--tumor-segmentation $out4 --contamination-table $out5"

	 		echo -e "Running GetPileupSummaries to summarize read support for a set number of known variant sites...\n+++++++++++++++++"
			gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" GetPileupSummaries -I $tumorBAM -L $targetInterval -V $variantsForContamination -O $out3 --interval-padding 100

			echo -e "Estimating contamination with CalculateContamination...\n+++++++++++++++++"
		    gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" CalculateContamination -I $out3 --tumor-segmentation $out4 -O $out5

		else 


	 		echo -e "Running GetPileupSummaries to summarize read support for a set number of known variant sites...\n+++++++++++++++++"
			gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" GetPileupSummaries -I $tumorBAM -L $targetInterval -V $variantsForContamination -O $out3 --interval-padding 100
			gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" GetPileupSummaries -I $normalBAM -L $targetInterval -V $variantsForContamination -O $out3Normal --interval-padding 100

			echo -e "Estimating contamination with CalculateContamination...\n+++++++++++++++++"
		    gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" CalculateContamination -I $out3 -matched $out3Normal --tumor-segmentation $out4 -O $out5
		fi 

	    echo -e "Passing the learned read orientation model to FilterMutectCalls with the -ob-priors argument...\n+++++++++++++++++"
		echo -e "skipsoftclip = $skipsoftclip"
		## 03/09/2019
		## A USER ERROR has occurred: dont-use-soft-clipped-bases is not a recognized option
		## changed to the following command...
		## 08/09/2019: from gatk 4.1.0.0 to 4.1.3.0, FilterMutectCalls requires -R param now... added
		## 05/05/2020: big changes regarding mutect2 variant filtering .... steps all changed. look above for notes!!
		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" FilterMutectCalls -R $refGenome -V $in $nomalContamString --ob-priors $out2 -O $out6 -L $targetInterval --interval-padding 100 --read-filter OverclippedReadFilter

		## 09/16/2021
		## add step after filtermutectcalls (Kyle): filteralnartifacts
		echo -e "Identify and filter alignment artifacts...\n+++++++++++++++++"
		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" FilterAlignmentArtifacts -R $refGenome -V $out6 -I $in2 --bwa-mem-index-image $refGenome.img -O $out7

		if [ $usegivenvcf -eq 1 ]; then 

			echo -e "Selecting variants that overlap with userVCF...\n+++++++++++++++++"
			gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" SelectVariants -R $refGenome -V $out7 --concordance $userVCF --output $out7.select.vcf.gz

			echo -e "Making 1-sample VCF for downstream analysis...\n+++++++++++++++++"
			gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" SelectVariants -R $refGenome -V $out7.select.vcf.gz --output $out7.select.1sm.vcf.gz -sn $tumor
		fi 

	fi 

	echo "## --------------------------------------"

	## 11/24/2020
	## mutect2multi: somatic calling multiple tumors
	## https://gatk.broadinstitute.org/hc/en-us/community/posts/360071839192-Mutect2-multi-sample-pipeline
	## https://gatk.broadinstitute.org/hc/en-us/community/posts/360062528691-mutect2-multi-sample-
	if [ $caller == 'mutect2multi' ] ; then 

		## even tho this is multi sample somatic mut calling, I still used the first
		## tumor sample as the output file name to be consistent with the pipeline design
		## this could be confusing but I just need to keep this in my mind! :)
		echo "[ " `date` " ] Somatic variant calling: $caller "
		out=$outputPath/somatic_calls/$caller/$tumor.$aligner.$caller${softclipOutstring}${ponOutstring}${emitgermOutstring}${givenvcfOutstring}.vcf.gz
		out2=$outputPath/somatic_calls/$caller/$tumor.$aligner.$caller${softclipOutstring}${ponOutstring}${emitgermOutstring}${givenvcfOutstring}.f1r2.tar.gz
		out3=${out/.vcf.gz/.bam}

		echo -e "Gathering input BAM and Names for $caller ...\n+++++++++++++++++"

		patient=`awk -F"\t" -v t=$tumor '$2==t {print $1}' $tumornormalpairs`
		tumorList=`awk -F"\t" -v p=$patient '$1==p && $3=="Tumor" {print $2}' $tumornormalpairs | perl -wpe 's/\s+/ /g'`
		normalList=`awk -F"\t" -v p=$patient '$1==p && $3=="Normal" {print $2}' $tumornormalpairs | perl -wpe 's/\s+/ /g'`
		echo -e "tumor = $tumor\tpatient = $patient\ttumorList = $tumorList\tnormalList = $normalList"

		## removed "-tumor $tumorName" as log returned "tumor name cannot be specified more than once"
		## As of v4.1, there is no longer a need to specify the tumor sample name with -tumor. You need only specify the normal sample name with -normal, if you include a normal.
		tumorInstring='' 
		for myt in $tumorList; do echo $myt; tumorBAM=''; tumorName=''; tumorBAM=$outputPath/../$myt/alignment/$myt.$aligner.merged${markdupString}.recal.bam; gatk GetSampleName -I $tumorBAM -O $tumorBAM.SampleName.txt; tumorName=`awk '{print $1}' $tumorBAM.SampleName.txt | perl -wpe 's/\s+//g'`; tumorInstring="$tumorInstring -I $tumorBAM"  ; done 

		normalInstring=''
		if [ $normalList == 'NA' ] || [ $normalList == '' ] || [ $normalList == '-' ]; then 
			normalInstring=''
		else 
			for myn in $normalList; do echo $myn; normalBAM=''; normalName=''; normalBAM=$outputPath/../$myn/alignment/$myn.$aligner.merged${markdupString}.recal.bam; gatk GetSampleName -I $normalBAM -O $normalBAM.SampleName.txt; normalName=`awk '{print $1}' $normalBAM.SampleName.txt | perl -wpe 's/\s+//g'`; normalInstring="$normalInstring -I $normalBAM -normal $normalName" ; done 
		fi 

		bamoutString=''
		if [ $mutect2BAM -eq 1 ]; then bamoutString=" --bam-output $out3"; fi 

		echo -e "tumorInstring = $tumorInstring\nnormalInstring = $normalInstring"

		echo -e "skipsoftclip = $skipsoftclip"
		echo -e "emitgerm = $emitgerm"

		echo -e "Starting $caller ...\n+++++++++++++++++"
		# gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" Mutect2 -R $refGenome -I $tumorBAM -tumor $tumorName $normalInstring -O $out -L $targetInterval --germline-resource $gnomad -mbq 20 -VS LENIENT --interval-padding 100 $softclipInstring $ponInstring $emitgermInstring $givenvcfInstring

		## 11/16/2019: error: Cannot construct fragment from more than two reads
		## temp solution .... https://gatkforums.broadinstitute.org/gatk/discussion/24492/mutect2-java-lang-illegalargumentexception-cannot-construct-fragment-from-more-than-two-reads#latest
		## 5/4/2020: seems this issue was fixed in gatk 4.1.5.0
		## 5/4/2020: remove --independent-mates in this run (gatk 4.1.7.0)
		## 5/5/2020: new strategy to filter mutect2 calls!
		## https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
		## need to run Mutect2 with the --f1r2-tar-gz argument in order to use the read orientation fitler.... have to rerun mutect2 jobs!!
		## 11/24/2020: implement mutect2 multi-tumor mode
		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" Mutect2 -R $refGenome $tumorInstring $normalInstring -O $out -L $targetInterval --germline-resource $gnomad -mbq 20 -VS LENIENT --interval-padding 100 --f1r2-tar-gz $out2 $softclipInstring $ponInstring $emitgermInstring $givenvcfInstring $bamoutString

		echo "## --------------------------------------"

		## mutect2: filter variants 
		## 5/5/2020: new strategy to filter mutect2 calls!
		## https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
		## need to run Mutect2 with the --f1r2-tar-gz argument in order to use the read orientation fitler.... have to rerun mutect2 jobs!!
		echo "[ " `date` " ] Somatic variant filtering: $caller "
		
		# normalBAM=$outputPath/../$normal/alignment/$normal.$aligner.merged${markdupString}.recal.bam
		fnVCF=$outputPath/somatic_calls/$caller/$tumor.$aligner.$caller${softclipOutstring}${ponOutstring}${emitgermOutstring}${givenvcfOutstring}
		in=$fnVCF.vcf.gz
		in2=$fnVCF.bam
		out1=$fnVCF.f1r2.tar.gz
		out2=$fnVCF.read-orientation-model.tar.gz
		# out3=$fnVCF.pileups.table
		# out4=$fnVCF.segments.table
		# out5=$fnVCF.contamination.table
		out6=$fnVCF.flt.vcf.gz
		out7=$fnVCF.flt.flt.vcf.gz

		echo -e "Passing f1r2.tar.gz from MuTect2 run to generate LearnReadOrientationModel...\n+++++++++++++++++"
 		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" LearnReadOrientationModel -I $out1 -O $out2

 		echo -e "Running GetPileupSummaries to summarize read support for a set number of known variant sites...\n+++++++++++++++++"

 		for myt in $tumorList; do echo $mys; tumorBAM=''; fnVCF=''; out3=''; tumorBAM=$outputPath/../$myt/alignment/$myt.$aligner.merged${markdupString}.recal.bam; fnVCF=$outputPath/somatic_calls/$caller/$myt.$aligner.$caller${softclipOutstring}${ponOutstring}${emitgermOutstring}${givenvcfOutstring}; out3=$fnVCF.pileups.table; gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" GetPileupSummaries -I $tumorBAM -L $targetInterval -V $variantsForContamination -O $out3 --interval-padding 100; done 

 		## Argument 'matched/matched-normal' cannot be specified more than once for the next step ... "CalculateContamination"
		matchInstring=''
 		if [ $normalList == 'NA' ] || [ $normalList == '' ] || [ $normalList == '-' ]; then 
			normalInstring=''
		else 

			for myn in $normalList; do echo $myn; normalBAM=''; fnVCF=''; out3=''; normalBAM=$outputPath/../$myn/alignment/$myn.$aligner.merged${markdupString}.recal.bam; fnVCF=$outputPath/somatic_calls/$caller/$myn.$aligner.$caller${softclipOutstring}${ponOutstring}${emitgermOutstring}${givenvcfOutstring}; out3=$fnVCF.pileups.table; gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" GetPileupSummaries -I $normalBAM -L $targetInterval -V $variantsForContamination -O $out3 --interval-padding 100; done

			## can only pick one normal sample ...
			myn=$normal
			fnVCF=''
			out3='' 
			if [ $patient == 'PHS13_31013' ]; then myn=S20_0024; fi 
			fnVCF=$outputPath/somatic_calls/$caller/$myn.$aligner.$caller${softclipOutstring}${ponOutstring}${emitgermOutstring}${givenvcfOutstring}
			out3=$fnVCF.pileups.table 
			matchInstring="$matchInstring -matched $out3"; 
		fi

		echo -e "Estimating contamination with CalculateContamination...\n+++++++++++++++++"

		tumorSegInstring=''
		contamTableInstring=''
 		for myt in $tumorList; do echo $mys; tumorBAM=''; fnVCF=''; out3=''; out4=''; out5=''; tumorBAM=$outputPath/../$myt/alignment/$myt.$aligner.merged${markdupString}.recal.bam; fnVCF=$outputPath/somatic_calls/$caller/$myt.$aligner.$caller${softclipOutstring}${ponOutstring}${emitgermOutstring}${givenvcfOutstring}; out3=$fnVCF.pileups.table; out4=$fnVCF.segments.table; out5=$fnVCF.contamination.table; gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" CalculateContamination -I $out3 $matchInstring --tumor-segmentation $out4 -O $out5; tumorSegInstring="$tumorSegInstring --tumor-segmentation $out4"; contamTableInstring="$contamTableInstring --contamination-table $out5"; done 

	    echo -e "Passing the learned read orientation model to FilterMutectCalls with the -ob-priors argument...\n+++++++++++++++++"
		echo -e "skipsoftclip = $skipsoftclip"
		## 03/09/2019
		## A USER ERROR has occurred: dont-use-soft-clipped-bases is not a recognized option
		## changed to the following command...
		## 08/09/2019: from gatk 4.1.0.0 to 4.1.3.0, FilterMutectCalls requires -R param now... added
		## 05/05/2020: big changes regarding mutect2 variant filtering .... steps all changed. look above for notes!!
		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" FilterMutectCalls -R $refGenome -V $in $tumorSegInstring $contamTableInstring --ob-priors $out2 -O $out6 -L $targetInterval --interval-padding 100 --read-filter OverclippedReadFilter

		## 09/16/2021
		## add step after filtermutectcalls (Kyle): filteralnartifacts
		echo -e "Identify and filter alignment artifacts...\n+++++++++++++++++"
		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" FilterAlignmentArtifacts -R $refGenome -V $out6 -I $in2 --bwa-mem-index-image $refGenome.img -O $out7

		if [ $usegivenvcf -eq 1 ]; then 

			echo -e "Selecting variants that overlap with userVCF...\n+++++++++++++++++"
			gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" SelectVariants -R $refGenome -V $out7 --concordance $userVCF --output $out7.select.vcf.gz

			echo -e "Making 1-sample VCF for downstream analysis...\n+++++++++++++++++"
			gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" SelectVariants -R $refGenome -V $out7.select.vcf.gz --output $out7.select.1sm.vcf.gz -sn $tumor
		fi 

	fi 

	echo "## --------------------------------------"

	## mutect2 PON: normal calling 06/23/2018
	## 05/03/2020: update commands to use GenomicsDBImport
	## https://gatk.broadinstitute.org/hc/en-us/articles/360042479112-CreateSomaticPanelOfNormals-BETA-
	## Note that as of May, 2019 -max-mnp-distance must be set to zero to avoid a known bug in GenomicsDBImport (Bad input: GenomicsDBImport does not support GVCFs with MNPs)
	## also need gatk 4.1.6.0 or higher; before this version even tho you set the param to zero mutect2 still emits mnp and GenomicsDBImport still fails!
	if [ $caller == 'mutect2pon' ] ; then 

		echo "[ " `date` " ] Making PON, tumor-only normal calling: $caller "
		normalBAM=$outputPath/../$normal/alignment/$normal.$aligner.merged${markdupString}.recal.bam
		out=$outputPath/somatic_calls/$caller/$normal.$aligner.$caller${softclipOutstring}.vcf.gz

		gatk GetSampleName -I $normalBAM -O $normalBAM.SampleName.txt
		normalName=`awk '{print $1}' $normalBAM.SampleName.txt | perl -wpe 's/\s+//g'`

		echo "normalName = $normalName"

		echo -e "skipsoftclip = $skipsoftclip"
		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" Mutect2 -R $refGenome -I $normalBAM -tumor $normalName -O $out -L $targetInterval -mbq 20 --interval-padding 100 $softclipInstring -max-mnp-distance 0

	fi 

	echo "## --------------------------------------"

	## lancet: somatic calling
	if [ $caller == 'lancet' ]; then 

		echo "[ " `date` " ] Somatic variant calling: $caller "
		tumorBAM=$outputPath/alignment/$tumor.$aligner.merged${markdupString}.recal.bam
		normalBAM=$outputPath/../$normal/alignment/$normal.$aligner.merged${markdupString}.recal.bam
		out=$outputPath/somatic_calls/$caller/$tumor.$aligner.$caller.$chr.vcf

		## by default the program outputs 'config.txt' in current running directory
		## to avoid overlapping of this file while multiple jobs are running....
		pwd=$PWD 
		cd $outputPath/somatic_calls/$caller 
		lancet --tumor $tumorBAM --normal $normalBAM --ref $refGenome --bed $targetBED.$chr.bed -X $threads -q 10 -C 17 -b 30 -c 5 -P 300 -a 3 -m 0 -e 0.04 -i 0 -o 4 -z 10 -v > $out 2> $out.log 
		bgzip -f $out
		tabix -p vcf -f $out.gz
		cd $pwd
	fi 
fi 

# =========================== PART VI ================================ ## 

if [ $part == 'III' ]; then 

	echo -e "tumor = $tumor; normal = $normal; annotator = $annotator"

	set -eo pipefail ## important! otherwise script will keep going despite failed

	echo "## --------------------------------------"

	if [ $annotator == 'annovar' ]; then 
		echo "[ " `date` " ] Running $annotator"

		for mystring in .flt .flt.flt; do 
			echo $mystring
			in=''
			out=''

			in=$outputPath/somatic_calls/$caller/$tumor.$aligner.$caller${softclipOutstring}${ponOutstring}${emitgermOutstring}${givenvcfOutstring}${mystring}.vcf.gz
			out=$outputPath/annotation/$annotator/$tumor.$aligner.$caller${softclipOutstring}${ponOutstring}${emitgermOutstring}${givenvcfOutstring}${givenvcfOutstring}${mystring}.avinput

			echo -e "Converting VCF to annovar input ..."
			## from annovar: use vcf4old format if VCF has more than 1 sample column!
			convert2annovar.pl -format vcf4old -includeinfo $in -outfile $out
			
			echo -e "Annotating variants ..."
			table_annovar.pl -thread $threads -maxgenethread $threads -otherinfo -buildver $genomeAnno -protocol $protocol -operation $operation -outfile $out $out ${HUMANDB} -remove

			echo -e "Replacing header ..."
			header=`zcat $in | grep '^#CHROM' `
			sed -i "s|\tOtherinfo$|\t$header|" $out.hg38_multianno.txt
		done 

	fi 

	echo "## --------------------------------------"

	if [ $annotator == 'funcotator' ]; then 

		# dataSourcesFolder=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/GDC/GATKbundle/funcotator/dataSources.v1.6.20190124s
		## 11/26/2020
		# dataSourcesFolder=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/GDC/GATKbundle/funcotator/funcotator_dataSources.v1.7.20200521s
		## 9/14/2024
		dataSourcesFolder=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/GDC/GATKbundle/funcotator/funcotator_dataSources.v1.8.hg38.20230908s

		genomeString=''
		if [ $genome == 'grch38' ]; then genomeString=hg38; fi 

		echo -e "dataSourcesFolder=$dataSourcesFolder\ngenomeString=$genomeString\ncenter=$center\ntumor=$tumor\nnormal=$normal\n"

		echo "[ " `date` " ] Running $annotator"
		## 9/16/2021: after adding filteralnartifact step after filtermutect2 step above, the .flt.flt.vcf.gz has way less calls than .flt.vcf.gz because only the PASS variants from filtermutect2 output were kept in filteralnartifact output with new flag [alignment] added
		for mystring in .flt .flt.flt; do 
			echo $mystring
			in=''
			out=''

			in=$outputPath/somatic_calls/$caller/$tumor.$aligner.$caller${softclipOutstring}${ponOutstring}${emitgermOutstring}${givenvcfOutstring}${mystring}.vcf.gz
			out=$outputPath/annotation/$annotator/$tumor.$aligner.$caller${softclipOutstring}${ponOutstring}${emitgermOutstring}${givenvcfOutstring}${givenvcfOutstring}${mystring}.$annotator

			tumor0=$tumor
			normal0=$normal
			if [ $caller == 'mutect2multi' ]; then 
				patient=`awk -F"\t" -v t=$tumor '$2==t {print $1}' $tumornormalpairs`
				tumor0=`awk -F"\t" -v p=$patient '$1==p && $3=="Tumor" {print $2}' $tumornormalpairs | perl -wpe 's/\s+/,/g'`
				normal0=`awk -F"\t" -v p=$patient '$1==p && $3=="Normal" {print $2}' $tumornormalpairs | perl -wpe 's/\s+/,/g'`
			echo -e "tumor = $tumor\tpatient = $patient\ttumor0 = $tumor0\tnormal0 = $normal0"
			fi 

			echo -e "Annotating variants ..."
			gatk Funcotator \
				--data-sources-path ${dataSourcesFolder}/ \
				-O $out.maf \
				--output-file-format MAF \
				--ref-version $genomeString \
				-R $refGenome \
				-V $in \
				--remove-filtered-variants true \
				--sequence-dictionary $refGenomeDict \
				--tmp-dir $tmpDir \
				--transcript-selection-mode CANONICAL \
				--annotation-default Center:$center \
				--annotation-default Tumor_Sample_Barcode:$tumor0 \
				--annotation-default Matched_Norm_Sample_Barcode:$normal0

			echo -e "Replacing double quotes ..."
			sed 's/"//g' $out.maf > $out.fixed.maf

		done 

	fi 

fi 

# =========================== PART V ================================ ## 

if [ $part == 'IV' ]; then 

	echo -e "project = $project"

	set -eo pipefail ## important! otherwise script will keep going despite failed

	echo "## --------------------------------------"

	## calculate relateness 
	echo "[ " `date` " ] relateness calculation: somalier "
	in=`ls $outputPath/../*/alignment_stats/somalier/*.somalier`
	out=$outputPath/../$project.somalier.relateness
	somalier relate --groups $somalierGroups --min-depth 7 --min-ab 0.3 --infer -o $out $in

	echo "## --------------------------------------"

	if [ $caller == 'mutect2pon' ]; then 		

		echo "[ " `date` " ] Creating PON VCF list file as input for next step"
		in=$tumornormalpairs
		out=$ponVCF.inlist

		echo -e "tumornormalpairs = $tumornormalpairs"
		if [ -e $out ]; then rm $out; fi 
		for sample in `awk -F"\t" 'NR>1 && $5=="Normal" && $2!="NA" {print $2}' $in`; do echo $sample; echo -e "$outputPath/../$sample/somatic_calls/$caller/$sample.$aligner.$caller.vcf.gz" >> $out;done 

		echo "## --------------------------------------"

		## GenomeImportDB extremely slow (running for 10 hours and not finish)
		## roll back to gatk 4.1.0.0, up until this version CreateSomaticPanelOfNormals has --vcfs option and  GenomeImportDB is NOT required
		## for the same list of vcfs GenomeImportDB took 20s to finish!
	    echo "[ " `date` " ] Combine MuTect2 normal calls to generate PON: CreateSomaticPanelOfNormals"
		in=$ponVCF.inlist
		out=$ponVCF.2.vcf.gz

		echo -e "out = $out; "`wc -l $in`" input normal vcfs"
		vcfileString=''
		for file in `cat $in`; do echo $file; vcfileString="$vcfileString --vcfs $file";done 
	    gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" CreateSomaticPanelOfNormals $vcfileString -O $out

	fi 

fi 

# =========================== PART VI ================================ ## 

if [ $part == 'V' ] && [ $caller == 'gatkhc' ]; then 

	echo -e "project = $project"

	set -eo pipefail ## important! otherwise script will keep going despite failed

	## still testing 
	# echo "## --------------------------------------"

	# echo "[ " `date` " ] Consolidating germline GVCFs from all normal samples: GenomicsDBImport"
	# in=$gVCF.inlist
	# out=$gVCF

	# echo -e "out = $out; "`wc -l $in`" input normal vcfs"
	# vcfileString=''
	# for file in `cat $in`; do echo $file; vcfileString="$vcfileString -vcfs $file";done 

	# gatk CreateSomaticPanelOfNormals $vcfileString -O $out 

	# echo "## --------------------------------------"

	# echo "[ " `date` " ] Joint-calling germline variants: GenotypeGVCFs"

	# java –jar GenomeAnalysisTK.jar –T GenotypeGVCFs\
	# 	–R human.fasta \
	# 	–V sample1.g.vcf \
	# 	–V sample2.g.vcf \
	# 	–V sampleN.g.vcf \
	# 	–o output.vcf

	# echo "## --------------------------------------"

	# echo "[ " `date` " ] Splitting germline vcf into snp and indel each: VariantRecalibrator, ApplyRecalibration"

	# echo "## --------------------------------------"

	# echo "[ " `date` " ] Filtering Variants by Variant (Quality Score) Recalibration: VariantRecalibrator, ApplyRecalibration"

	# echo -e "Recalibrate SNPs..."
	# java –jar GenomeAnalysisTK.jar –T VariantRecalibrator \
	# 	–R human.fasta \
	# 	–input raw.SNPs.vcf \
	# 	-resource:hapmap,known=false,training=true,truth=true,prior=15.0
	# 	hapmap_3.3.b37.sites.vcf \
	# 	-resource:omni,known=false,training=true,truth=false,prior=12.0
	# 	omni2.5.b37.sites.vcf \
	# 	-resource:1000G,known=false,training=true,truth=false,prior=10.0
	# 	1000G.b37.sites.vcf \
	# 	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0
	# 	dbsnp_137.b37.vcf \
	# 	–an DP –an QD –an FS –an MQRankSum {...} \
	# 	–mode SNP \
	# 	–recalFile raw.SNPs.recal \
	# 	–tranchesFile raw.SNPs.tranches \
	# 	–rscriptFile recal.plots.R

	# java –jar GenomeAnalysisTK.jar –T ApplyRecalibra8on \
	# 	–R human.fasta \
	# 	–input raw.SNPs.vcf \
	# 	–mode SNP \
	# 	–recalFile raw.SNPs.recal \
	# 	–tranchesFile raw.SNPs.tranches \
	# 	–o recal.SNPs.vcf \
	# 	–ts_filter_level 99.0


fi 

echo "## --------------------------------------"

echo "[ " `date` " ] Program finished!"




