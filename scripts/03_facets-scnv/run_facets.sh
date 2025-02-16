
echo "hostname = $HOSTNAME"
hostname=htc ## stats or hpc; if hpc for some reason samtools cannot be found; only use stats even for hpc jobs

now=$(date +"%m-%d-%Y_%H:%M:%S")
	
runATM=0 ## only use ATM polymorphisms....

## --------------------------------------
## command line options 
## --------------------------------------

project=PDO_sm61 ## ATM_PAAD MSLNpilot
sample=$1 ## MTS-1
genome=$2 ## grch38 grch37 hg19
part=$3 ## I or II
aligner=$4 ## bwamem
caller=$5 ## mutect2
tumor=$6
normal=$7

echo "Sample = $sample; part = $part"
if [ -z $sample ] || [ -z $genome ] || [ -z $aligner ]; then echo -e "Input missing. \n$0 <sample> <genome>\n"; exit; fi 

## --------------------------------------
## Tools
## --------------------------------------

if [ $hostname == 'htc' ]; then 
	modules='gcc/8.2.0 bedtools/2.29.0 samtools/1.12 r/4.1.0'
	export PATH=$PATH:/ix/rbao/Software/facets/0.6.2/inst/extcode
fi 

## --------------------------------------
## Settings
## --------------------------------------

threads=16

## --------------------------------------
## Paths and directories
## --------------------------------------

projPath=/ix/jluke/Projects/sunny/TIIL-017-RoseProj
outputPath="$projPath/results/exome/${project}_samples_${genome}/$sample"
dirlist=" scnv_calls scnv_calls/gatkscnv"

## --------------------------------------
## Commands
## --------------------------------------

for module in $modules; do module load $module; done 
module list 

echo "## --------------------------------------"

cd $projPath
pwd=$PWD
echo "Current Path = $pwd"

echo -e "Running FACETS: tumor = $tumor; normal = $normal"

if [ ! -d $outputPath/scnv_calls/facets ]; then mkdir -p $outputPath/scnv_calls/facets; fi 

which samtools
which bedtools
which R

echo -e "[" `date` "] Generate SNP pileup: snp-pileup"

tumorBAM=$outputPath/alignment/$tumor.bwamem.merged.rmdup.bam
normalBAM=$outputPath/../$normal/alignment/$normal.bwamem.merged.rmdup.bam
out=$outputPath/scnv_calls/facets/$tumor.bwamem.merged.rmdup.bam.snp_pileup.csv.gz

## https://github.com/dariober/cnv_facets
## A VCF file of common, polymorphic SNPs. For human samples, a good source is the dbSNP file common_all.vcf.gz (https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/). See also NCBI human variation sets in VCF Format.
## however from my prev notes (2021) I used gnomad as below, so I am gonna continue using this one
knownVCF=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/GDC/broad-references/hg38/v0/af-only-gnomad.hg38.vcf.gz.WGSprimary_AF0.02.vcf.gz

## use all variants within ATM loci regardless of AF...
if [ $runATM -eq 1 ]; then 
	knownVCF=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/GDC/broad-references/hg38/v0/af-only-gnomad.hg38.vcf.gz.atm.vcf.gz
	out=$outputPath/scnv_calls/facets/$tumor.bwamem.merged.rmdup.bam.snp_pileup.atm.csv.gz
fi 

## if not removing it then will receive error msg as 'Output file /ix/jluke/Projects/sunny/TIIL-017-RoseProj/results/exome/PDO_sm61_samples_grch38/TIIL_0165/scnv_calls/facets/TIIL_0165.bwamem.merged.rmdup.bam.snp_pileup.csv.gz already exists!'
rm $out

snp-pileup -d4000 -g -q30 -Q20 -P100 -r25,0 -v $knownVCF $out $normalBAM $tumorBAM

## ------------------------------------------

## downstream analysis done in R, on my local computer
## facets, and facets-suite are both R packages
## https://github.com/mskcc/facets
## https://github.com/mskcc/facets-suite
## https://github.com/dariober/cnv_facets (seems most recent?)

echo "## --------------------------------------"

echo "[ " `date` " ] Program finished!"




