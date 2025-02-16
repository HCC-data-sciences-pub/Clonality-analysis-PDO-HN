#!/bin/bash
#
#SBATCH -N 1 
#SBATCH -c 16
#SBATCH --mem 36g
#SBATCH -t 0-46:00:00

echo "hostname = $HOSTNAME"

now=$(date +"%m-%d-%Y_%H:%M:%S")

projPath=/ix/jluke/Projects/sunny/TIIL-017-RoseProj

project=PDO_sm61
sample=MTS12
genome=grch38
part=I
aligner=bwamem
caller=facets
tumor=
normal=
script=scripts/run_facets.sh
log=logs/batches/run_facets.$sample.$genome.$now.log

echo START `date`
cd $projPath
./$script $sample $genome $part $aligner $caller $tumor $normal >& $log 
echo END `date`

