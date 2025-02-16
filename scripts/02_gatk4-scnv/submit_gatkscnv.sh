#!/bin/bash
#
#SBATCH -N 1 
#SBATCH -c 4
#SBATCH --mem 60g
#SBATCH -t 0-22:00:00

echo "hostname = $HOSTNAME"

now=$(date +"%m-%d-%Y_%H:%M:%S")

projPath=/ix/jluke/Projects/sunny/TIIL-017-RoseProj

project=PDO_sm61
sample=PNATM-101-tumor
aligner=bwamem
caller=gatkscnv
datatype=exome
script=scripts/run_gatkscnv.sh
log=logs/batches/run_gatkscnv.$datatype.$sample.$now.log

echo START `date`
cd $projPath
./$script $sample $aligner $caller $datatype >& $log 
echo END `date`
