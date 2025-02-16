#!/bin/bash
#
#SBATCH -N 1 
#SBATCH -c 8
#SBATCH --mem 16g
#SBATCH -t 0-96:00:00

echo "hostname = $HOSTNAME"

now=$(date +"%m-%d-%Y_%H:%M:%S")

projPath=/ix/jluke/Projects/sunny/TIIL-017-RoseProj

project=PDO_sm31B2
sample=S024
genome=grch38
part=I
aligner=bwamem
caller=mutect2
tumor=
normal=
annotator=
script=scripts/run_exome.sh
log=
chr=

if [ $part == 'I' ]; then 
	log=logs/batches/run_exome.$sample.$genome.$part.$aligner.$now.log
fi 
if [ $part == 'II' ]; then 
	log=logs/batches/run_exome.$sample.$genome.$part.$aligner.$caller.$now.log

	if [ $caller == 'lancet' ]; then 
		log=logs/batches/run_exome.$sample.$genome.$part.$aligner.$caller.$chr.$now.log
	fi 
fi 
if [ $part == 'III' ]; then 
	log=logs/batches/run_exome.$sample.$genome.$part.$aligner.$caller.$annotator.$now.log

fi 
if [ $part == 'IV' ]; then 
	log=logs/batches/run_exome.$project.$genome.$part.$aligner.$caller.$now.log

fi 

echo START `date`
cd $projPath
./$script $sample $genome $part $aligner $caller $tumor $normal $annotator $chr >& $log 
echo END `date`

