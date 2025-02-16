
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
				echo $sample $genome $aligner
				script=scripts/submit_exome.sh
				# script=scripts/submit_bamreadcount.sh ## only part II
				# script=scripts/submit_pureCN.sh ## only part I and III
				# script=scripts/submit_sequenza.sh ## only part II
				# script=scripts/submit_freec.sh ## only part II
				# config=scripts/freec.gMTS100.cfg
				gender=''
				sex=XY
				scriptsNew=
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
					echo "x=$x" 
					if [ $x -ne 1 ]; then 
						echo "skip"
						continue
					fi

				fi


				# gender=`awk -F"\t" -v s=$sample '$1==s {print $5}' /gpfs/data/bioinformatics/Projects/CRI-BIO-630-RO-MSpiotto-rbao/sampleinfo/STEP_2_LNM_seq-revised_for_Riyue.mike-1.2.gender.txt`
				# if [ $gender == 'F' ]; then sex=XX; fi 
				# echo $gender $sex

				if [ $part == 'I' ]; then 
					scriptsNew=scripts/batches/submit_exome.$sample.$genome.$part.$aligner.sh
					# scriptsNew=scripts/batches/submit_pureCN.$sample.$genome.$part.$aligner.sh
				fi 
				if [ $part == 'II' ] || [ $part == 'III' ]; then 
					tumorFlag=`awk -v s=$sample '$1==s {print $1}' $tumornormalpairs | sort | uniq | wc -l`
					normalFlag=`awk -v s=$sample '$2==s {print $2}' $tumornormalpairs | sort | uniq | wc -l`

					echo -e "$tumorFlag\t$normalFlag"

					if [ $caller == 'mutect2pon' ]; then 
						if [ $normalFlag -eq 1 ]; then 
							## could be multiple tumors !!! only pick one
							tumor=`awk -v s=$sample '$2==s {print $1}' $tumornormalpairs | sort | uniq | head -1`
							normal=$sample

							echo -e "this is normal! $tumor, $normal"
							scriptsNew=scripts/batches/submit_exome.$sample.$genome.$part.$aligner.$caller.sh

						else  
							continue 
						fi  
					else 
						if [ $tumorFlag -eq 1 ]; then 
							tumor=$sample
							normal=`awk -v s=$sample '$1==s {
								print $2}' $tumornormalpairs | sort | uniq`

							echo -e "this is tumor! $tumor, $normal"
							scriptsNew=scripts/batches/submit_exome.$sample.$genome.$part.$aligner.$caller.sh
							# scriptsNew=scripts/batches/submit_bamreadcount.$tumor.$genome.sh
							# scriptsNew=scripts/batches/submit_sequenza.$tumor.$genome.sh
							# scriptsNew=scripts/batches/submit_freec.$tumor.$genome.sh
							# configNew=scripts/batches/freec.$tumor.cfg

							if [ $part == 'III' ]; then 
								scriptsNew=scripts/batches/submit_exome.$sample.$genome.$part.$aligner.$caller.$annotator.sh
								# scriptsNew=scripts/batches/submit_pureCN.$sample.$genome.$part.$aligner.sh
							fi 

						else 
							continue 
						fi 
					fi 
					
				fi 
				if [ $part == 'IV' ]; then 
					scriptsNew=scripts/batches/submit_exome.$project.$genome.$part.$aligner.$caller.sh
				fi

				# echo $script $scriptsNew

				cp $script $scriptsNew
				sed -i "s|^sample=.*|sample=$sample|;s|^genome=.*|genome=$genome|;s|^aligner=.*|aligner=$aligner|;s|^caller=.*|caller=$caller|;s|^part=.*|part=$part|" $scriptsNew

				if [ $part == 'II' ] && [ $tumorFlag -eq 1 ]; then 
					sed -i "s|^tumor=.*|tumor=$tumor|;s|^normal=.*|normal=$normal|;s|^annotator=.*|annotator=$annotator|" $scriptsNew


					# cp $config $configNew
					# sed -i "s|gMTS100|$tumor|g;s|^sex =.*|sex = $sex|" $configNew 
					
					if [ $caller == 'lancet' ]; then 
						for chr in $chrs; do 
							scriptsNew2=${scriptsNew/.sh/.$chr.sh}
							# echo $chr
							# echo $scriptsNew2
							cp $scriptsNew $scriptsNew2
							sed -i "s|^chr=.*|chr=$chr|;s|^threads=.*|threads=4|;s|mem=.*|mem=8g|;s|ppn=.*|ppn=4|;s|walltime=.*|walltime=24:00:00|" $scriptsNew2
						done 
						rm $scriptsNew
					fi 
				fi 
				if [ $part == 'II' ] && [ $normalFlag -eq 1 ]; then 

					
					sed -i "s|^tumor=.*|tumor=$tumor|;s|^normal=.*|normal=$normal|;s|^annotator=.*|annotator=$annotator|" $scriptsNew

				fi 
				if [ $part == 'III' ] && [ $tumorFlag -eq 1 ]; then 
					sed -i "s|^tumor=.*|tumor=$tumor|;s|^normal=.*|normal=$normal|;s|^annotator=.*|annotator=$annotator|;s|-t 0-.*|-t 0-6:00:00|" $scriptsNew
				fi

				# if [ $aligner == 'bwamem' ]; then sed -i "s|mem=.*|mem=12G|;s|walltime=.*|walltime=24:00:00|" $scriptsNew; fi 
			done 

		done 
	done 
done 


chmod u+x scripts/batches/*.sh

## tcga raw wes data files .... convert them first!!
# for sample in `cat sampleinfo/TCGA_UVM_exome.filelist`; do echo $file; uuid=`echo $sample | sed 's/\//\t/g' | cut -f 8`; echo $uuid; script=scripts/submit_biobambam.sh; scriptNew=scripts/batches/submit_biobambam.$uuid.sh; cp  $script $scriptNew; sed -i "s|^sample=.*|sample=$sample|" $scriptNew ;done


## for tumor only freec condig file .... make these changes manually!!
## delete entire [control] section 
## set intercept = 0 (default is 1)
# for sample in gMTS112 gMTS17 gMTS19 gMTS24 gMTS39 gMTS6 ; do sed -i 's/^intercept = 1/intercept = 0/' scripts/batches/freec.$sample.cfg;done
