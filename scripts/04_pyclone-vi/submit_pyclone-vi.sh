#!/bin/bash
#
#SBATCH -N 1 
#SBATCH -c 8
#SBATCH --mem 24g
#SBATCH -t 0-46:00:00

echo "hostname = $HOSTNAME"

now=$(date +"%m-%d-%Y_%H:%M:%S")

projPath=/ix/jluke/Projects/sunny/TIIL-017-RoseProj

project=PDO_sm61
genome=grch38
sample=TIIL_0165
part=II
script=scripts/run_pyclone-vi.sh
log=logs/batches/run_pyclone-vi.$sample.$genome.$now.log

echo START `date`
cd $projPath
./$script $part $sample >& $log 
echo END `date`

