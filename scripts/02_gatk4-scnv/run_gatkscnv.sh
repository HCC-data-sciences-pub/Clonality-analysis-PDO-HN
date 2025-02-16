
echo "hostname = $HOSTNAME"
hostname=htc ## stats or hpc

now=$(date +"%m-%d-%Y_%H:%M:%S")

## -------------------------------------

set -eo pipefail ## important! otherwise script will keep going despite failed

## -------------------------------------

## remember to change several places where I labeled [PDO_sm30]
project=PDO_sm61 ## PDO_sm30
sample=$1 ## MTS-1
aligner=$2 ## bwamem
caller=$3 ## mutect2
datatype=$4 ## exome genome 

echo "Sample = $sample; Datatype = $datatype"
if [ -z $sample ] || [ -z $aligner ] || [ -z $caller ] || [ -z $datatype ]; then echo -e "Input missing. \n$0 <sample> <aligner> <caller> <datatype>\n"; exit; fi 

echo "## --------------------------------------"

if [ $hostname == 'htc' ]; then 
    ## 9/13/2024: gatk/4.6.0.0 requires higher java version; load java/21.0.2-openjdk resolves this
    modules='gcc/8.2.0 htslib/1.9 bwa/0.7.17 fastqc/0.11.7 bedtools/2.29.0 bedops/2.4.35 bzip2/1.0.6 xz/5.2.3 r/3.6.0 sambamba/0.6.8 samtools/1.9  igv/2.4.16 bcftools/1.11 java/21.0.2-openjdk'

    ## annovar/20180416 lancet/1.0.6 abra/2.12 flash/1.2.11 

    ## 9/13/2024: update gatk to 4.6.0.0; also update funcotator bundle 
    export PATH=/ix/rbao/Software/UCSCtools:/ix/rbao/Software/star/2.7.3a/bin/Linux_x86_64_static:/ix/rbao/Software/kallisto/0.46.1:/ix/rbao/Software/salmon/1.0.0/bin:/ix/rbao/Software/bamreadcount/0.8.0/bin:/ix/rbao/Software/mosdepth/0.2.6:/ix/rbao/Software/gatk/4.6.0.0:/ix/rbao/Software/annovar/20180416:/ix/rbao/Software/somalier_2/0.2.10:/ix/rbao/Software/tabix/0.2.6:$PATH
    PICARD=/ihome/crc/install/picard/2.18.12/picard.jar
    TRIMMOMATIC=/ihome/crc/install/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar
    HUMANDB=/ix/rbao/Software/annovar/humandb

    ## fastqc returned error:
    ## /usr/bin/perl: symbol lookup error: /ihome/crc/install/gcc-8.2.0/perl/5.28.0/lib/5.28.0/x86_64-linux-thread-multi/auto/Cwd/Cwd.so: undefined symbol: Perl_xs_handshake
    ## solution: https://github.com/bioconda/bioconda-recipes/issues/4390
    export PERL5LIB=""
fi 

echo "## --------------------------------------"

threads=16
genome=grch38
minContigLength=1000000 ## Threshold length (in bp) for contigs to be plotted. Contigs with lengths less than this threshold will not be plotted. This can be used to filter out mitochondrial contigs, unlocalized contigs, etc.
padding=0
binLength=0
minAF=0.02 ## keep common variants to reduce gnomad size ...
assayName=''
javamem=''

## https://gatk.broadinstitute.org/hc/en-us/articles/360035531092
## The default value of 250 bases was determined to work well empirically for TCGA targeted exome data. This argument is relevant for exome data, as binning without an intervals list does not allow for intervals expansion
if [ $datatype == 'exome' ]; then 
    #assayName=S07604624_V6UTRr2
    if [ $project == 'PDO_sm30' ] || [ $project == 'PDO_sm61' ]; then assayName=S33266340; fi
    padding=250; binLength=0; javamem="-Xmx12G"
fi 
if [ $datatype == 'genome' ]; then 
    assayName="WGSprimary_AF$minAF"
    padding=0; binLength=1000; javamem="-Xmx30G"
fi 

echo "## --------------------------------------"

projPath=/ix/jluke/Projects/sunny/TIIL-017-RoseProj
outputPath="$projPath/results/${datatype}/${project}_samples_$genome/$sample"
dirlist=" scnv_calls scnv_calls/gatkscnv"

tmpDir=$outputPath/tmp

if [ ! -d $tmpDir ]; then mkdir -p $tmpDir; fi
for dir in $dirlist; do echo $dir; if [ ! -d $outputPath/$dir ]; then mkdir -p $outputPath/$dir; fi  ;done 

echo "## --------------------------------------"

sampleinfo=''
if [ $datatype == 'exome' ]; then 
    if [ $project == 'PDO_sm30' ]; then 
        sampleinfo=$projPath/sampleinfo/$project.exome.tumor_normal.20241223.gatkscnv.tsv
        # sampleinfo=$projPath/sampleinfo/$project.exome.tumor_normal.20241223.gatkscnv.reorder.tsv ## for gistic2 heatmap

    fi
    if [ $project == 'PDO_sm61' ]; then 
        sampleinfo=$projPath/sampleinfo/$project.exome.tumor_normal.20250102.gatkscnv.reorder.tsv
        ## pooling both batches: sm61 (30+31)
    fi
fi 
if [ $datatype == 'genome' ]; then 
    sampleinfo=$projPath/sampleinfo/$project.genome.tumor_normal.20200526.delly.tsv
fi 

echo "## --------------------------------------"

refGenome=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/GDC/GRCh38.d1.vd1/GRCh38.d1.vd1.fa
refGenomeDict=$refGenome.dict
chromSize=$refGenome.chromsize
bwaIndex=$refGenome
# GATKbundle=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/GDC/GATKbundle/hg38/v0
broadref=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/GDC/broad-references/hg38/v0
dbSNP=$broadref/dbsnp_146.hg38.vcf.gz
G1000snp=$broadref/1000G_phase1.snps.high_confidence.hg38.vcf.gz
hapmap=$broadref/hapmap_3.3.hg38.vcf.gz
mills=$broadref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
omni=$broadref/1000G_omni2.5.hg38.vcf.gz
axiom=$broadref/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
knownIndels=$broadref/Homo_sapiens_assembly38.known_indels.vcf.gz
cosmic=/group/bioinformatics/ReferenceData/cosmic/v83/cosmic_coding_and_noncoding_chr_M_sorted.vcf
gnomad=$broadref/af-only-gnomad.hg38.vcf.gz
variantsForContamination=$broadref/small_exac_common_3.hg38.vcf.gz
## if it was exome ...
if [ $datatype == 'exome' ]; then 
    # targetBED=/ix/jluke/Projects/sunny/TIIL-001-Augustin-Proj/target/S07604514_V6r2_Covered.picard.hg38.merged.srt.exRandom_Alt.bed ## AMMelanoma
    if [ $project == 'PDO_sm30' ] || [ $project == 'PDO_sm61' ]; then 
        targetBED=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/Agilent/S33266340_hg38/S33266340_Covered.exRandom_Alt.bed
    fi 
fi 
## if it was genome ... pick the primary chrs only!!!!
if [ $datatype == 'genome' ]; then 
    targetBED=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/GDC/GRCh38.d1.vd1/GRCh38.d1.vd1.fa.manta.bed
fi 
targetInterval=$targetBED.interval_list

## blacklist, most recent
## paper: https://www.nature.com/articles/s41598-019-45839-z
## data: https://github.com/Boyle-Lab/Blacklist/tree/master/lists
## direct link: https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
blacklist=/ix/rbao/Projects/HCC-CBS-000-ReferenceFiles/Homo_sapiens/Blacklist/hg38-blacklist.v2.bed

echo "## --------------------------------------"

# echo -e "patient = $patient\ntumor = $tumor\nnormal = $normal"
echo -e "sample = $sample"

echo "## --------------------------------------"

for module in $modules; do module load $module;done 
module list

echo "## --------------------------------------"

##  for each target bed file...
echo -e "targetBED = $targetBED"

## only run once if target bed file stays the same ....
## if target bed changes, need to re-run this!!
if [ 0 -eq 1 ]; then 

    echo "[ " `date` " ] 1. Collect raw counts data with PreprocessIntervals and CollectFragmentCounts"
    gatk --java-options "-XX:ParallelGCThreads=2 $javamem -Djava.io.tmpdir=$tmpDir" PreprocessIntervals \
        -L $targetInterval \
        -R $refGenome \
        --bin-length $binLength \
        --padding $padding \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O $targetInterval.$caller.preprocessed.interval_list
fi 

echo "## --------------------------------------"

## only run once if target bed file stays the same ....
## if target bed changes, need to re-run this!!

## https://gatkforums.broadinstitute.org/gatk/discussion/23207/what-interval-should-i-give-to-collectalleliccounts-for-somatic-cnv-in-exome
## https://gatkforums.broadinstitute.org/gatk/discussion/11683#ref9
## https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants
## For targeted exomes, it may be convenient to subset these to the preprocessed intervals, e.g. with SelectVariants for use with CollectAllelicCounts. This is not necessary, however, as ModelSegments drops sites outside the target regions from its analysis in the joint-analysis approach.
## For whole genomes, depending on the desired resolution of the analysis, consider subsetting the gnomAD sites to those commonly variant, e.g. above an allele frequency threshold. Note that SelectVariants, as of this writing, can filter on AF allele frequency only for biallelic sites. Non-biallelic sites make up ~3% of the gnomAD SNPs-only resource.
if [ 0 -eq 1 ]; then 

    echo "[ " `date` " ] Generate a smaller VCF from gnomad as input for step 5, otherwise taking too much time to read the 3G vcf.gz"
    
    echo -e "subset variants ..."
    if [ $datatype == 'exome' ]; then 

        gatk --java-options "-XX:ParallelGCThreads=2 ${javamem} -Djava.io.tmpdir=$tmpDir" SelectVariants \
            -L $targetInterval \
            --interval-padding $padding \
            -V $gnomad \
            -R $refGenome \
            -O $gnomad.$assayName.vcf.gz

    fi 

    if [ $datatype == 'genome' ]; then

        ## need to filter for primary chrs ... 
        ## otherwise gatk complains "Badly formed genome unclippedLoc: Contig chr1_KI270766v1_alt given as location, but this contig isn't present in the Fasta sequence dictionary" in downstream  CollectAllelicCounts step!!!
        gatk --java-options "-XX:ParallelGCThreads=2 ${javamem} -Djava.io.tmpdir=$tmpDir" SelectVariants \
            -L $targetInterval \
            --interval-padding $padding \
            -V $gnomad \
            -R $refGenome \
            -O $gnomad.$assayName.tmp.vcf.gz

        ## okay, if I did not subset this the CollectAllelicCounts step will never finish even with 60G mem
        ## I could never get VariantFiltration in gatk to work
        ## I give up. Use bcftools instead 
        bcftools filter \
            --include "AF > $minAF" \
            -o $gnomad.$assayName.vcf.gz -O z \
            --threads $threads \
            $gnomad.$assayName.tmp.vcf.gz

        tabix -p vcf -f $gnomad.$assayName.vcf.gz
        rm $gnomad.$assayName.tmp.vcf.gz
    fi 

    echo -e "convert vcf to interval list ...(using -L with vcf files could be very slow"
    java -XX:ParallelGCThreads=2 ${javamem} -Djava.io.tmpdir=$tmpDir -jar ${PICARD} VcfToIntervalList \
        I=$gnomad.$assayName.vcf.gz \
        O=$gnomad.$assayName.vcf.gz.interval_list 

fi 

echo "## --------------------------------------"

## https://gatk.broadinstitute.org/hc/en-us/articles/360035531092?id=11682
echo -e "GATK-SCNV detection part I: per sample "

if [ 0 -eq 1 ]; then 

    ## run this for both tumor and control samples ...
    echo -e "sample = $sample"

    echo "[ " `date` " ] 1. Count ref and alt alleles at common germline variant sites using CollectAllelicCounts"
    in=$outputPath/alignment/${sample}.bwamem.merged.rmdup.bam
    out=$outputPath/scnv_calls/$caller/${sample}.$aligner.$caller.counts.hdf5
    gatk --java-options "-XX:ParallelGCThreads=2 ${javamem} -Djava.io.tmpdir=$tmpDir" CollectReadCounts \
        -I $in \
        -L $targetInterval.$caller.preprocessed.interval_list \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O $out

fi 

echo "## --------------------------------------"

echo -e "GATK-SCNV detection part I: panel of normal across all normals"

if [ 0 -eq 1 ]; then 

    # ## 12/26/2024: comment this part out!! otherwise job halts with error msg below and never proceed 
    # 
    # ## https://github.com/bioconda/bioconda-recipes/issues/12944
    # ## gatk scnv needs openblas 0.3.7 or higher but it is not bundled within gatk, so need to load it seperately! otherwise you will receive error such as  libopenblas.so cannot be found
    # ## in addition to load the module, you also need to explicitely export the path to libopenblas.so ....
    # module load gcc/4.8.5 openblas/0.3.10
    # export LD_PRELOAD=/ihome/crc/install/openblas/0.3.10/gcc-4.8.5/lib/libopenblas.so
    #
    # ## returned error ... 12/26/2024
    # ## Exception: java.lang.StackOverflowError thrown from the UncaughtExceptionHandler in thread "process reaper"
    # ## have to comment the above part out; now the script can run

    ## --------------------------------------------

    echo "[ " `date` " ] 2. Generate a CNV panel of normals with CreateReadCountPanelOfNormals"
    out=$outputPath/../${project}.$aligner.$caller.pon.hdf5
    inString=''
    normalList=`awk -F"\t" '$2=="control" {print $1}' $sampleinfo | perl -wpe 's/\s+/ /g'`
    # normalList="PNATM-101-normal PNATM-102-normal"
    # echo -e "tumor = $tumor\tpatient = $patient\tnormalList = $normalList"
    echo -e "sample = $sample\tnormalList = $normalList"

    for myn in $normalList; do echo $myn; in=''; in=$outputPath/../${myn}/scnv_calls/$caller/${myn}.$aligner.$caller.counts.hdf5; inString="$inString -I $in"; done 
    echo "inString = $inString"

    gatk --java-options "-XX:ParallelGCThreads=2 ${javamem} -Djava.io.tmpdir=$tmpDir"  CreateReadCountPanelOfNormals \
        $inString \
        --minimum-interval-median-percentile 5.0 \
        -O $out

fi 

echo "## --------------------------------------"

echo -e "GATK-SCNV detection part I and II: per tumor vs PoN"

if [ 1 -eq 1 ]; then

    ## only run this for tumor samples ...
    tumor=$sample
    normal=

    echo -e "tumor = $tumor\nnormal = $normal"

    echo "[ " `date` " ] 3. Standardize and denoise case read counts against the PoN with DenoiseReadCounts"
    in=$outputPath/scnv_calls/$caller/${tumor}.$aligner.$caller.counts.hdf5
    out1=$outputPath/scnv_calls/$caller/${tumor}.$aligner.$caller.standardizedCR.tsv
    out2=$outputPath/scnv_calls/$caller/${tumor}.$aligner.$caller.denoisedCR.tsv
    gatk --java-options "-XX:ParallelGCThreads=2 -Xmx12G -Djava.io.tmpdir=$tmpDir"  DenoiseReadCounts \
        -I $in \
        --count-panel-of-normals $outputPath/../${project}.$aligner.$caller.pon.hdf5 \
        --standardized-copy-ratios $out1 \
        --denoised-copy-ratios $out2

    echo "## --------------------------------------"

    echo "[ " `date` " ] 4. Plot standardized and denoised copy ratios with PlotDenoisedCopyRatios"

    in1=$outputPath/scnv_calls/$caller/${tumor}.$aligner.$caller.standardizedCR.tsv 
    in2=$outputPath/scnv_calls/$caller/${tumor}.$aligner.$caller.denoisedCR.tsv
    out1=$outputPath/scnv_calls/$caller
    out2=${tumor}.$aligner.$caller
    gatk --java-options "-XX:ParallelGCThreads=2 ${javamem} -Djava.io.tmpdir=$tmpDir" PlotDenoisedCopyRatios \
        --standardized-copy-ratios $in1 \
        --denoised-copy-ratios $in2 \
        --sequence-dictionary $refGenomeDict \
        --minimum-contig-length $minContigLength \
        --output $out1 \
        --output-prefix $out2

    echo "## --------------------------------------"

    echo -e "GATK-SCNV detection part II:"

    echo "## --------------------------------------"

    echo "[ " `date` " ] 5. Count ref and alt alleles at common germline variant sites using CollectAllelicCounts"
    echo "Collect counts at germline variant sites for the tumor ..."
    in=$outputPath/alignment/$tumor.bwamem.merged.rmdup.bam
    out=$outputPath/scnv_calls/$caller/${tumor}.$aligner.$caller.allelicCounts.tsv
    gatk --java-options "-XX:ParallelGCThreads=2 ${javamem} -Djava.io.tmpdir=$tmpDir" CollectAllelicCounts \
        -L $gnomad.$assayName.vcf.gz.interval_list \
        -I $in \
        -R $refGenome \
        -O $out

    if [ ! -z $normal ] && [ $normal != '' ]; then 

        echo "Collect counts at the same sites for the matched-control ..."
        in=$outputPath/../${normal}/alignment/$normal.bwamem.merged.rmdup.bam
        out=$outputPath/scnv_calls/$caller/${normal}.$aligner.$caller.allelicCounts.tsv
        gatk --java-options "-XX:ParallelGCThreads=2 ${javamem} -Djava.io.tmpdir=$tmpDir"  CollectAllelicCounts \
            -L $gnomad.$assayName.vcf.gz.interval_list \
            -I $in \
            -R $refGenome \
            -O $out

    fi 

    echo "## --------------------------------------"

    echo "[ " `date` " ] 6. Group contiguous copy ratios into segments with ModelSegments"
    in1=$outputPath/scnv_calls/$caller/${tumor}.$aligner.$caller.denoisedCR.tsv
    in2=$outputPath/scnv_calls/$caller/${tumor}.$aligner.$caller.allelicCounts.tsv
    in3=$outputPath/scnv_calls/$caller/${normal}.$aligner.$caller.allelicCounts.tsv
    out1=$outputPath/scnv_calls/$caller
    out2=${tumor}.$aligner.$caller
    normalString=''
    if [ ! -z $normal ] && [ $normal != '' ]; then normalString=" --normal-allelic-counts $in3"; fi 
    echo -e "normalString = $normalString"
    gatk --java-options "-XX:ParallelGCThreads=2 ${javamem} -Djava.io.tmpdir=$tmpDir"  ModelSegments \
        --denoised-copy-ratios $in1 \
        --allelic-counts $in2 $normalString \
        --output $out1 \
        --output-prefix $out2

    echo "## --------------------------------------"

    echo "[ " `date` " ] 7. Call copy-neutral, amplified and deleted segments with CallCopyRatioSegments"
    in=$outputPath/scnv_calls/$caller/${tumor}.$aligner.$caller.cr.seg 
    out=$outputPath/scnv_calls/$caller/${tumor}.$aligner.$caller.called.seg
    gatk --java-options "-XX:ParallelGCThreads=2 ${javamem} -Djava.io.tmpdir=$tmpDir" CallCopyRatioSegments \
        --input $in \
        --output $out

    echo "## --------------------------------------"

    echo "[ " `date` " ] 8. Plot modeled copy ratio and allelic fraction segments with PlotModeledSegments"
    in1=$outputPath/scnv_calls/$caller/${tumor}.$aligner.$caller.denoisedCR.tsv
    in2=$outputPath/scnv_calls/$caller/${tumor}.$aligner.$caller.hets.tsv
    in3=$outputPath/scnv_calls/$caller/${tumor}.$aligner.$caller.modelFinal.seg
    out1=$outputPath/scnv_calls/$caller
    out2=${tumor}.$aligner.$caller
    gatk --java-options "-XX:ParallelGCThreads=2 ${javamem} -Djava.io.tmpdir=$tmpDir" PlotModeledSegments \
        --denoised-copy-ratios $in1 \
        --allelic-counts $in2 \
        --segments $in3 \
        --sequence-dictionary $refGenomeDict \
        --minimum-contig-length $minContigLength \
        --output $out1 \
        --output-prefix $out2

fi 

echo "## --------------------------------------"

echo "[ " `date` " ] Program finished!"
