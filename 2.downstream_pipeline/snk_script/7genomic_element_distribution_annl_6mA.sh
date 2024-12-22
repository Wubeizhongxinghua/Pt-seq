#!/bin/bash
#bedadd=$1
workdir=$1
shift
samples=("$@")
#bed=$1
#bedname=$2
mkdir -p ${workdir}/distribution

for sample in ${samples[@]}
do
	echo "This sample is ${sample}"
	echo "First one is ${samples[0]}"
	if [[ "$sample" = "${samples[0]}" ]]; then #process for ctrl
		/home/chengqiyi_pkuhpc/profiles/limingyang/software/homer/bin/annotatePeaks.pl ${workdir}/${sample}_ctrl_stopsite.bed hg38 -annStats ${workdir}/distribution/ctrl_elementStat.txt > ${workdir}/distribution/ctrl_annotatePeaks_Stats.txt

	#awk 'BEGIN{k=0} {if(k==2 && $1=="Annotation"){exit;} if($1=="Annotation"){k+=2;} {print $0;}}' ./distribution/${bed}_elementStat_${con}.txt > ./distribution/${bed}_broadelement_${con}.txt
		#print detailed annotation
		awk 'BEGIN{k=0}
			{
				if($1=="Annotation"){
					k += 1;
				}
			
				if(k==2){
					print $0; 
				}
			}
		' ${workdir}/distribution/ctrl_elementStat.txt > ${workdir}/distribution/ctrl_broadelement.txt
	fi

	/home/chengqiyi_pkuhpc/profiles/limingyang/software/homer/bin/annotatePeaks.pl ${workdir}/${sample}_treat.final_stopsite.bed hg38 -annStats ${workdir}/distribution/${sample}_elementStat.txt > ${workdir}/distribution/${sample}_annotatePeaks_Stats.txt

#awk 'BEGIN{k=0} {if(k==2 && $1=="Annotation"){exit;} if($1=="Annotation"){k+=2;} {print $0;}}' ./distribution/${bed}_elementStat_${con}.txt > ./distribution/${bed}_broadelement_${con}.txt
	awk 'BEGIN{k=0}
		{
			if($1=="Annotation"){
				k += 1;
			}
		
			if(k==2){
				print $0;
			}
		}
	' ${workdir}/distribution/${sample}_elementStat.txt > ${workdir}/distribution/${sample}_broadelement.txt
#python3 1make_table_for_plot_by_con.py -i ./distribution -b ${bed} -o ./distribution/${bed}_accCon_forplot.txt
#	Rscript 1plot_con.r -i ./distribution/${bed}_accCon_forplot.txt -o ./distribution/${bed}_accCon.pdf
done
##### divide subtypes
#mkdir -p ./distribution/subtypes/
#
#for con in hyper hypo
#do
##cat ${bedadd} | grep ${subtype} > ${bedadd}_${subtype}.bed
#	annotatePeaks.pl ${bedadd}__${con}.bed hg38 -annStats ./distribution/subtypes/${bed}_${subtype}_${con}_elementStat.txt > ./distribution/subtypes/${bed}_${subtype}_${con}_annotatePeaks_Stats.txt
#	awk 'BEGIN{k=0} {if(k==1 && $1=="Annotation"){exit;} if($1=="Annotation"){k+=1;} {print $0;}}' ./distribution/subtypes/${bed}_${subtype}_${con}_elementStat.txt > ./distribution/subtypes/${bed}_${subtype}_${con}_broadelement.txt
#done

/home/chengqiyi_pkuhpc/profiles/limingyang/miniforge3/bin/Rscript snk_script/7plot_subtype_persample.r -i ${workdir}/distribution -o ${workdir}/distribution/annotation.pdf
