
projPath=/ix/jluke/Projects/sunny/TIIL-017-RoseProj

## exome
project=PDO_sm31B2 ## PDO_sm30
# metadata=$projPath/sampleinfo/${project}.exome.metadata.20240913.txt
# tumornormalpairs=$projPath/sampleinfo/${project}.exome.tumor_normal.20240913.list
metadata=$projPath/sampleinfo/${project}.exome.metadata.20241231.txt
tumornormalpairs=$projPath/sampleinfo/${project}.exome.tumor_normal.20241231.list
tumornormalpairs0=$tumornormalpairs
genomes=' grch38' ## grch37 grch38
aligners=' bwamem'
annotator=funcotator ## annovar funcotator
part=III ## I to IV
callers='mutect2' ## lancet mutect2 mutect2_2 mutect2pon mutect2multi
chrs='chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr2 chr20 chr21 chr22 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chrX chrY chrUn_KI270742v1 '

if [ $part == 'IV' ]; then callers='mutect2pon'; fi 

for sample in `grep -v '^Sample' $metadata | cut -f 1 | sort | uniq`; do 
	for genome in $genomes; do 
		for aligner in $aligners; do 
			for caller in $callers; do 
				# echo $sample $genome $aligner
				script=
				tumorFlag=
				tumor=
				normal=

				## mutect2_2 is only because one patient has 2 normal samples .. I could not decide which one to use. So I am running 2 analyses using each as normal sample... 
				if [ $caller == 'mutect2_2' ]; then tumornormalpairs=$tumornormalpairs0.2; fi 
				if [ $part == 'II' ] & [ $caller == 'mutect2_2' ]; then
					x=0 
					x=`awk -F"\t" -v s=$sample '$1==s' $tumornormalpairs | wc -l`;
					# echo "x=$x" 
					if [ $x -eq 0 ]; then 
						echo "skip"
						continue
					fi

				fi 

				## mutect2 multi 
				if [ $caller == 'mutect2multi' ]; then tumornormalpairs2=$tumornormalpairs0.multi; fi 
				## I designed the run script to anchor on the 1st tumor sample from each patient
				if [ $caller == 'mutect2multi' ]; then 
					x=0
					rm $tumornormalpairs2.tmp
					x=`awk -F"\t" -v s=$sample '$2==s{print $1}' $tumornormalpairs2`
					awk -F"\t" -v x=$x '$1==x' $tumornormalpairs2 > $tumornormalpairs2.tmp
					x=`awk -F"\t" -v s=$sample '$2==s {print NR}' $tumornormalpairs2.tmp`
					# echo "x=$x" 
					if [ $x -ne 1 ]; then 
						# echo "skip"
						continue
					fi

				fi
				if [ $part == 'I' ]; then 
					script=scripts/batches/submit_exome.$sample.$genome.$part.$aligner.sh
					# script=scripts/batches/submit_pureCN.$sample.$genome.$part.$aligner.sh
				fi 
				if [ $part == 'IV' ]&& [ $caller == 'mutect2pon' ]; then 
					scriptsNew=scripts/batches/submit_exome.$project.$genome.$part.$aligner.$caller.sh
				fi
				if [ $part == 'II' ] || [ $part == 'III' ]; then 
					tumorFlag=`awk -v s=$sample '$1==s {print $1}' $tumornormalpairs | sort | uniq | wc -l`
					normalFlag=`awk -v s=$sample '$2==s {print $2}' $tumornormalpairs | sort | uniq | wc -l`

					# echo -e "$tumorFlag\t$normalFlag"

					if [ $caller == 'mutect2pon' ]; then 
						if [ $normalFlag -eq 1 ]; then 
							## could be multiple tumors !!! only pick one
							tumor=`awk -v s=$sample '$2==s {print $1}' $tumornormalpairs | sort | uniq | head -1`
							normal=$sample

							# echo -e "this is tumor! $tumor, $normal"
							script=scripts/batches/submit_exome.$sample.$genome.$part.$aligner.$caller.sh

						else  
							continue 
						fi  
					else 
						if [ $tumorFlag -eq 1 ]; then 
							tumor=$sample
							normal=`awk -v s=$sample '$1==s {print $2}' $tumornormalpairs | sort | uniq`

							# echo -e "this is tumor! $tumor, $normal"
							script=scripts/batches/submit_exome.$sample.$genome.$part.$aligner.$caller.sh
							# script=scripts/batches/submit_bamreadcount.$tumor.$genome.sh
							# script=scripts/batches/submit_sequenza.$tumor.$genome.sh
							# script=scripts/batches/submit_freec.$tumor.$genome.sh

							if [ $part == 'III' ]; then 
								script=scripts/batches/submit_exome.$sample.$genome.$part.$aligner.$caller.$annotator.sh
								# script=scripts/batches/submit_pureCN.$sample.$genome.$part.$aligner.sh
							fi 

						else  
							continue 
						fi 
					fi 
				fi 

				if [ $part == 'I' ]; then if [ $sample == 'TIIL_0260' ]; then continue; fi; fi 
				if [ $part == 'II' ] ; then if [ $sample == 'TIIL_0260' ]; then continue; fi; fi
				if [ $part == 'II' ] && [ $caller == 'mutect2pon' ]  ; then if [ $sample == 'TIIL_0260' ]; then continue; fi; fi
				if [ $part == 'III' ] ; then if [ $sample == 'TIIL_0260' ]; then continue; fi; fi
					
				if [ $part == 'II' ] && [ $sample == 'lancet' ]; then 
					for chr in $chrs; do 
						script2=${script/.sh/.$chr.sh}

						echo $script2
						# sbatch $script2
						# sleep 2
					done 
				else 
					if [ ! -z "$script" ] && [ $part != 'IV' ]; then 
						echo $script
						sbatch $script
						sleep 2
					fi 
				fi 



				# if [ $aligner == 'bwamem' ]; then sed -i "s|mem=.*|mem=12G|;s|walltime=.*|walltime=24:00:00|" $scriptsNew;  && [ $sample == 'S001' ]fi 
			done 

		done 
	done 
done 

if [ $part == 'IV' ]; then 
	script=scripts/batches/submit_exome.$project.$genome.$part.$aligner.$caller.sh
	echo $script
# 	sbatch $script
fi 

## tcga raw rnaseq data files .... untar them first!!
# for sample in `cat sampleinfo/TCGA_UVM_exome.filelist`; do echo $file; uuid=`echo $sample | sed 's/\//\t/g' | cut -f 8`; echo $uuid; script=scripts/batches/submit_biobambam.$uuid.sh; echo $script; sbatch $script; sleep 2; done


chmod u+x scripts/batches/*.sh

