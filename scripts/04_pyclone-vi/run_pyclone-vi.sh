
## You cannot load the module twice. This is the nature of modules using biocontainers. You can run "module purge" and then "module load pyclone-vi/0.1.6" in your slurm script to avoid this problem. 

echo "hostname = $HOSTNAME"
hostname=htc 
now=$(date +"%m-%d-%Y_%H:%M:%S")

## ------------------------------------------

## global settings
project=PDO_sm61 ## ATM_PAAD MSLNpilot
genome=grch38
cnvcaller=facets ## gatkscnv, facets
part=$1 ## I II III IV IV 
sample=$2 ## for some parts of analysis below: part II

if [ -z $part ]; then echo -e "$0 <part> [<sample>]"; echo "Param missing. Program exit"; exit; fi 

## ------------------------------------------

## tools
modules='gcc/8.2.0 bcftools/1.15.1 bedtools/2.29.0 bedops/2.4.35 bzip2/1.0.6 xz/5.2.3 r/4.1.0 sambamba/0.6.8 samtools/1.12 igv/2.4.16 java/21.0.2-openjdk'
export PATH=/ix/rbao/Software/UCSCtools:/ix/rbao/Software/bamreadcount/0.8.0/bin:/ix/rbao/Software/mosdepth/0.2.6:/ix/rbao/Software/gatk/4.6.0.0:/ix/rbao/Software/somalier_2/0.2.15:/ix/rbao/Software/pigz/2.4:/ix/rbao/Software/tabix/0.2.6/:$PATH
PICARD=/ihome/crc/install/picard/2.18.12/picard.jar

## ------------------------------------------

## parameters
threads=16

## ------------------------------------------

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

	somalierSites=/ix/rbao/Software/somalier_2/0.2.15/sites_files/sites.hg38.vcf.gz
fi

## ------------------------------------------

## path
projPath=/ix/jluke/Projects/sunny/TIIL-017-RoseProj
outputPath="$projPath/results/exome/${project}_samples_${genome}"

## ------------------------------------------

## files 
metadata=$projPath/sampleinfo/${project}.exome.metadata.20250102.txt
tumornormallist=$projPath/sampleinfo/${project}.exome.tumor_normal.20250102.list

## ------------------------------------------

echo START `date`
echo -e "Part = $part"

echo -e "\n===STEP 0 - Prepare environment===\n"
if [ 1 -eq 1 ] ; then 

	module purge ## fangping
	
	module load $modules 
	module list 

	cd $projPath
	pwd=$PWD
	echo "Current Path = $pwd"

	patients=`awk -F"\t" 'NR>1 {print $4}' $tumornormallist | sort | uniq | perl -wpe 's/\s+/ /g'`
	echo -e "Patients = $patients"

fi 

## ------------------------------------------

echo -e "\n===STEP 1 - create united set of known VCFs from mutect2 calls across samples===\n"

echo -e "STEP 1.1 - filter mutect2 VCF to keep PASS variants: gatk SelectVariants"
if [ $part == 'I' ] & [ 0 -eq 1 ]; then 

	for pt in $patients; do 
		echo $pt

		tumors=''
		tumors=`awk -F"\t" -v pt=$pt '$4==pt{print $1}' $tumornormallist | perl -wpe 's/\s+/ /g'`
		normals=''
		normals=`awk -F"\t" -v pt=$pt '$4==pt{print $2}' $tumornormallist | sort | uniq | perl -wpe 's/\s+/ /g'`

		echo -e "Tumors = $tumors\nNormals = $normals"

		vcflist=''
		for sample in $tumors; do 
			echo $sample

			## select only PASS variants ...
			in=${outputPath}/${sample}/somatic_calls/mutect2/${sample}.bwamem.mutect2.wpon.flt.flt.vcf.gz
			out=${outputPath}/${sample}/somatic_calls/mutect2/${sample}.bwamem.mutect2.wpon.flt.flt.pass.vcf.gz
			tmpDir=${outputPath}/${sample}/tmp
			if [ ! -d $tmpDir ]; then mkdir $tmpDir; fi 

			gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" SelectVariants -R $refGenome -V $in --exclude-filtered true -O $out

		done 

	done 

fi 

echo -e "STEP 1.2 - merge all PASS variants into one united set [not intersect!!] across all samples: bcftools merge"
if [ $part == 'I' ] & [ 0 -eq 1 ]; then 

	vcflistAll='' ## all patients
	for pt in $patients; do 
		echo $pt

		## I had manually sorted the sample list to always list source tumor (NOT PDO) as 1st sample of the same patient in the tumor normal list file
		tumorSource=`awk -F"\t" -v pt=$pt '$4==pt{print $1}' $tumornormallist | head -1`

		echo -e "Source Tumor = $tumorSource"

		tumors=''
		tumors=`awk -F"\t" -v pt=$pt '$4==pt{print $1}' $tumornormallist | perl -wpe 's/\s+/ /g'`
		normals=`awk -F"\t" -v pt=$pt '$4==pt{print $2}' $tumornormallist | sort | uniq | perl -wpe 's/\s+/ /g'`

		echo -e "Tumors = $tumors\nNormals = $normals"

		vcflist='' ## per patient 
		for sample in $tumors; do 
			echo $sample

			in=${outputPath}/${sample}/somatic_calls/mutect2/${sample}.bwamem.mutect2.wpon.flt.flt.pass.vcf.gz

			# vcflist="$vcflist $in"
			# vcflistAll="$vcflistAll $in"

			## 1/20: okay so ... downstream you can see I had to exclude multiallic and INDELs only keep biallic and SNPs for gatk ASE run. however if I merge VCFs first then filter on merged VCF, sites with overlapping INDELs and SNPs will show as multiallic in merged VCF, and excluded by the biallic filter. As a result, lots sites were excluded.
			## so I decided to exclude INDELS on per sample VCF first before merging to see if I can include more sites ...

			out=${outputPath}/${sample}/somatic_calls/mutect2/${sample}.bwamem.mutect2.wpon.flt.flt.pass.snp.vcf.gz
			tmpDir=${outputPath}/${sample}/tmp
			if [ ! -d $tmpDir ]; then mkdir $tmpDir; fi 

			echo -e "Filter out INDELs sites before merging: $pt, $sample"
			gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" SelectVariants -R $refGenome -V $in --select-type-to-include SNP -O $out
			tabix -f -p vcf $out

			vcflist="$vcflist $out"
			vcflistAll="$vcflistAll $out"

		done 

		echo -e "VCFlist = $vcflist"

		## running bcftools ...
		echo -e "Save output in the source tumor folder per patient: $pt, $tumorSource"
		out=${outputPath}/${tumorSource}/somatic_calls/mutect2/${pt}.bwamem.mutect2.wpon.flt.flt.pass.merge.00.vcf.gz
		bcftools merge --force-samples $vcflist | bgzip > $out 
		# bcftools merge --force-samples $vcflist -m none | bgzip > $out 
		tabix -f -p vcf $out

		## 1/20 update: the merged VCF did not work for gatk ASE run downsteam, got truncated output 
		## error msg: "gatk ASE run: A USER ERROR has occurred: More then one variant context at position: 
		## solution: supply with [-m none] in bcftools merge above: "no new multiallelics, output multiple records instead"; so I re-ran this step, and manually confirmed that there is no ',' (if multiallelic it would be something like [A,T] etc. ) in ALT field in the merged VCF
		## more update: NOPE that did not fix it!! still got truncated gatk ASE output with same error lol ## "gatk ASE run: A USER ERROR has occurred: More then one variant context at position: chr3:142843144"
		## this is why: https://github.com/broadinstitute/gatk/issues/7249
		## It turns out it's caused either by the input VCF having 2 SNPs on separate rows but with the same position, or by an INDEL that overlaps a SNP
		## The solution was to filter out INDELs and multiallelic sites before running ASEReadCounter
		## so I reverted back to not using [-m none], so that bcftools will merge by default
		## Grrrrrrrrrrrr >=(

		echo -e "Filter out INDELs and multiallelic sites before running ASEReadCounter: $pt, $tumorSource"
		in=${outputPath}/${tumorSource}/somatic_calls/mutect2/${pt}.bwamem.mutect2.wpon.flt.flt.pass.merge.00.vcf.gz
		# out1=${outputPath}/${tumorSource}/somatic_calls/mutect2/${pt}.bwamem.mutect2.wpon.flt.flt.pass.merge.vcf.0.gz
		out2=${outputPath}/${tumorSource}/somatic_calls/mutect2/${pt}.bwamem.mutect2.wpon.flt.flt.pass.merge.vcf.gz
		tmpDir=${outputPath}/${sample}/tmp
		if [ ! -d $tmpDir ]; then mkdir $tmpDir; fi 

		# gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" SelectVariants -R $refGenome -V $in --restrict-alleles-to BIALLELIC --select-type-to-include SNP -O $out1
		# tabix -f -p vcf $out1
		# ## note bcftools rm dup behavior...
		# ## https://github.com/samtools/bcftools/issues/1089
		# bcftools norm --rm-dup exact $out1 | bgzip > $out2
		# tabix -f -p vcf $out2

		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" SelectVariants -R $refGenome -V $in --restrict-alleles-to BIALLELIC --select-type-to-include SNP -O $out2
		tabix -f -p vcf $out2

	done 

	echo -e "vcflistAll = $vcflistAll"

	## running bcftools ...
	echo -e "Save output in the cohort folder across all pts: $outputPath"
	out=${outputPath}/${project}.bwamem.mutect2.wpon.flt.flt.pass.merge.00.vcf.gz
	bcftools merge --force-samples $vcflistAll | bgzip > $out 
	# bcftools merge --force-samples $vcflistAll -m none | bgzip > $out 
	tabix -f -p vcf $out

	in=${outputPath}/${project}.bwamem.mutect2.wpon.flt.flt.pass.merge.00.vcf.gz
	out2=${outputPath}/${project}.bwamem.mutect2.wpon.flt.flt.pass.merge.vcf.gz
	tmpDir=${outputPath}/${sample}/../tmp
	if [ ! -d $tmpDir ]; then mkdir $tmpDir; fi 
	gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" SelectVariants -R $refGenome -V $in --restrict-alleles-to BIALLELIC --select-type-to-include SNP -O $out2
	tabix -f -p vcf $out2

	## ==============
	## test purpose only : how many variants were excluded by each filter?
	gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" SelectVariants -R $refGenome -V $in --restrict-alleles-to BIALLELIC -O ${outputPath}/${project}.bwamem.mutect2.wpon.flt.flt.pass.merge.tmp1.vcf.gz

	gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" SelectVariants -R $refGenome -V $in --select-type-to-include SNP -O ${outputPath}/${project}.bwamem.mutect2.wpon.flt.flt.pass.merge.tmp2.vcf.gz

	# PDO_sm61.bwamem.mutect2.wpon.flt.flt.pass.merge.00.vcf.gz
	# 12684
	# PDO_sm61.bwamem.mutect2.wpon.flt.flt.pass.merge.tmp1.vcf.gz
	# 12684
	# PDO_sm61.bwamem.mutect2.wpon.flt.flt.pass.merge.tmp2.vcf.gz
	# 12684
	# PDO_sm61.bwamem.mutect2.wpon.flt.flt.pass.merge.vcf.gz
	# 12684

	## note that the rows stay the same before and after filters above... it is because I already excluded indels before merging :) and it turned out that all SNPs are uniq sites and no multiallelic sites were generated after merging per sample VCFs
	## I have manually checked the filters above and they DID WORK if I did not filter indels before merging first.
	## ==============


fi

## ------------------------------------------

echo -e "\n===STEP 2 - collect ref and alt allelic specific counts on known VCF===\n"
## submit jobs per sample: loop too slow!
if [ $part == 'II' ] & [ 0 -eq 1 ]; then 

	echo -e "Sample = $sample"

	if [ ! -z $sample ]; then 

		pt=''
		pt=`awk -F"\t" -v sm=$sample '$1==sm{print $4}' $tumornormallist | head -1`

		echo -e "Patient = $pt"
		if [ -z $pt ]; then echo "Patient is missing. Program exit!"; exit; fi 

		## ------------------------------------------

		## I had manually sorted the sample list to always list source tumor (NOT PDO) as 1st sample of the same patient in the tumor normal list file
		tumorSource=`awk -F"\t" -v pt=$pt '$4==pt{print $1}' $tumornormallist | head -1`

		echo -e "Source Tumor = $tumorSource"

		tumors=''
		tumors=`awk -F"\t" -v pt=$pt '$4==pt{print $1}' $tumornormallist | perl -wpe 's/\s+/ /g'`
		normals=`awk -F"\t" -v pt=$pt '$4==pt{print $2}' $tumornormallist | sort | uniq | perl -wpe 's/\s+/ /g'`

		echo -e "Tumors = $tumors\nNormals = $normals"

		knownVCF1=${outputPath}/${tumorSource}/somatic_calls/mutect2/${pt}.bwamem.mutect2.wpon.flt.flt.pass.merge.vcf.gz
		knownVCF2=${outputPath}/${project}.bwamem.mutect2.wpon.flt.flt.pass.merge.vcf.gz

		## ------------------------------------------

		in=${outputPath}/${sample}/alignment/${sample}.bwamem.merged.rmdup.bam
		out1=${outputPath}/${sample}/alignment/${sample}.bwamem.merged.rmdup.bam.knownSitesPt.txt
		out2=${outputPath}/${sample}/alignment/${sample}.bwamem.merged.rmdup.bam.knownSitesAll.txt
		tmpDir=${outputPath}/${sample}/tmp
		if [ ! -d $tmpDir ]; then mkdir $tmpDir; fi 

		echo -e "Count ref and alt allele counts: gatk ASEReadCounter, known sites=per pt"
		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" ASEReadCounter -R $refGenome -I $in -V $knownVCF1 -O $out1 -L $targetInterval -mbq 20 --lenient --interval-padding 100 

		echo -e "Count ref and alt allele counts: gatk ASEReadCounter, known sites=all pts"
		gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir" ASEReadCounter -R $refGenome -I $in -V $knownVCF2 -O $out2 -L $targetInterval -mbq 20 --lenient --interval-padding 100 

	fi 


fi 
## ------------------------------------------

echo -e "\n===STEP 3 - collect copy number on known VCF sites===\n"
if [ $part == 'III' ] & [ 1 -eq 1 ]; then 

	knownVCF=${outputPath}/${project}.bwamem.mutect2.wpon.flt.flt.pass.merge.vcf.gz

	## ------------------------------------------

	if [ $cnvcaller == 'facets' ]; then 

		# [rib37@htc-n36 PDO_sm61_multisample_grch38]$ head -1 PDO_sm61.exome.facetsSuite.txt | sed 's/\t/\n/g' | awk '{print NR"\t"$0}'
		# 1       Tumor <<<<<<<<<<<<<<<
		# 2       purity
		# 3       ploidy
		# 4       dipLogR
		# 5       flags
		# 6       em_flags
		# 7       chrom <<<<<<<<<<<<<<<
		# 8       seg
		# 9       num.mark
		# 10      nhet
		# 11      cnlr.median
		# 12      mafR
		# 13      segclust
		# 14      cnlr.median.clust
		# 15      mafR.clust
		# 16      start <<<<<<<<<<<<<<<
		# 17      end <<<<<<<<<<<<<<<
		# 18      cf.em
		# 19      tcn.em 
		# 20      lcn.em 
		# 21      cf
		# 22      tcn
		# 23      lcn

		# [rib37@htc-n36 PDO_sm61_multisample_grch38]$ head -1 PDO_sm61.exome.facets.fit.wType.txt | sed 's/\t/\n/g' | awk '{print NR"\t"$0}'
		# 1       Tumor
		# 2       purity
		# 3       ploidy
		# 4       dipLogR
		# 5       chrom
		# 6       seg
		# 7       num.mark
		# 8       nhet
		# 9       cnlr.median
		# 10      mafR
		# 11      segclust
		# 12      cnlr.median.clust
		# 13      mafR.clust
		# 14      start
		# 15      end
		# 16      cf.em
		# 17      tcn.em
		# 18      lcn.em
		# 19      type


		## facets documentation in R:
		## In the output, cf, tcn, lcn are the initial estimates of cellular fraction, total and minor copy number estimates, and cf.em, tcn.em, lcn.em are the estimates by the mixture model optimized using the EM-algorithm. cf is used as initial values for the EM algorithm. For diploid normal segments (total copy=2, minor copy=1), we report cellular fraction as 1 (100% normal). The logOR data for a segment are summarized using the square of expected log-odds-ratio (mafR column). Estimated tumor sample purity and ploidy are reported.

		## Is there a way to create outputs for Pyclone??
		## https://github.com/mskcc/facets/issues/57
		echo -e "Convert tsv segs into bed files and sort: bedtools"
		# in=$projPath/results/exome/${project}_multisample_${genome}/${project}.exome.facetsSuite.txt
		# out=$projPath/results/exome/${project}_multisample_${genome}/${project}.exome.facetsSuite.srt.bed
		# out2=$projPath/results/exome/${project}_multisample_${genome}/${project}.exome.facetsSuite.srt.bed.wknownVCF.bed

		in=$projPath/results/exome/${project}_multisample_${genome}/${project}.exome.facets.fit.wType.txt
		out=$projPath/results/exome/${project}_multisample_${genome}/${project}.exome.facets.fit.wType.bed
		out2=$projPath/results/exome/${project}_multisample_${genome}/${project}.exome.facets.fit.wType.bed.wknownVCF.bed

		## facets is 1-based and bed is 0-based
		awk -F"\t" 'NR>1{print "chr"$5"\t"$14-1"\t"$15"\t"$1":"$5":"$14":"$15"\t.\t+"}' $in | bedtools sort -i stdin > $out

		echo -e "Collect copy number from segs: bedtools intersect, known sites=all pts"
		bedtools intersect -a $knownVCF -b $out -wao | bedtools sort -i stdin > $out2
		cut -f 1-7,114-120 $out2 > $out2.slim
		# cut -f 1-7,114-120 $out2 | sed 's/:/\t/g' > $out2.slim

		echo -e "Convert tsv segs into bed files and sort: bedtools"
		in=$projPath/results/exome/${project}_multisample_${genome}/${project}.exome.facets.fit.wType.LOH.txt
		out=$projPath/results/exome/${project}_multisample_${genome}/${project}.exome.facets.fit.wType.LOH.bed
		out2=$projPath/results/exome/${project}_multisample_${genome}/${project}.exome.facets.fit.wType.LOH.bed.wknownVCF.bed

		## facets is 1-based and bed is 0-based
		awk -F"\t" 'NR>1{print "chr"$5"\t"$14-1"\t"$15"\t"$1":"$5":"$14":"$15"\t.\t+"}' $in | bedtools sort -i stdin > $out

		echo -e "Collect copy number from LOH segs: bedtools intersect, known sites=all pts"
		bedtools intersect -a $knownVCF -b $out -wao | bedtools sort -i stdin > $out2
		cut -f 1-7,114-120 $out2 > $out2.slim
		# cut -f 1-7,114-120 $out2 | sed 's/:/\t/g' > $out2.slim

	fi 


fi 

## ------------------------------------------

echo -e "\n===STEP 4 - create input for pyclone-vi by merging sommut and scnv results===\n"

## run this on my local computer: prepare_pyclone-vi_input.R

## ------------------------------------------

echo -e "\n===STEP 5 - run pyclone-vi===\n"

if [ $part == 'V' ] & [  0 -eq 1 ] ; then 

	patientsTumorPDO=" HN23-10730 HN24-10810 HN24-10856 HN24-10900 HN24-10909 HN24-10936 "

	echo -e "patientTumorPDO = $patientsTumorPDO"

	## ------------------------------------------

	module load pyclone-vi/0.1.6

	# dataInPath='pyclone-vi_altDP0'
	dataInPath='pyclone-vi_altDP3'
	# dataInPath='pyclone-vi_altDP8'
	echo -e "$dataInPath"

	for pt in $patientsTumorPDO; do 
		echo $pt

		myTag='Pt'
		in=${projPath}/results/exome/${project}_multisample_${genome}/${dataInPath}/${myTag}/${project}.gatkASE.pyclone.list.${pt}.${myTag}.txt
		out=$in.out.h5
		out2=$in.out.txt
		
		echo -e "Run pyclone-vi fit: $pt, $myTag"
		# pyclone-vi fit -i $in -o $out -c 40 -d beta-binomial -r 10
		pyclone-vi fit -i $in -o $out -c 40 -d beta-binomial -r 50

		echo -e "Run pyclone-vi write-results-file: $pt, $myTag"
		pyclone-vi write-results-file -i $out -o $out2

	done 

	

fi 

## ------------------------------------------

echo END `date`
