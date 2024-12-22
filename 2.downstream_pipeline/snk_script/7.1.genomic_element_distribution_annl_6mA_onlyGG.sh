#!/bin/bash
#bedadd=$1
workdir=$1
shift
samples=("$@")
#bed=$1
#bedname=$2
mkdir -p ${workdir}/distribution_onlyGG

for sample in ${samples[@]}
do
	echo "This sample is ${sample}"
	echo "First one is ${samples[0]}"
	if [[ "$sample" = "${samples[0]}" ]]; then #process for ctrl
		/home/chengqiyi_pkuhpc/profiles/limingyang/software/homer/bin/annotatePeaks.pl ${workdir}/${sample}_GG_ctrl_stopsite.bed hg38 -annStats ${workdir}/distribution_onlyGG/ctrl_elementStat.txt > ${workdir}/distribution_onlyGG/ctrl_annotatePeaks_Stats.txt

	#awk 'BEGIN{k=0} {if(k==2 && $1=="Annotation"){exit;} if($1=="Annotation"){k+=2;} {print $0;}}' ./distribution_onlyGG/${bed}_elementStat_${con}.txt > ./distribution_onlyGG/${bed}_broadelement_${con}.txt
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
		' ${workdir}/distribution_onlyGG/ctrl_elementStat.txt > ${workdir}/distribution_onlyGG/ctrl_broadelement.txt
	fi

	/home/chengqiyi_pkuhpc/profiles/limingyang/software/homer/bin/annotatePeaks.pl ${workdir}/${sample}_GG_treat.final_stopsite.bed hg38 -annStats ${workdir}/distribution_onlyGG/${sample}_elementStat.txt > ${workdir}/distribution_onlyGG/${sample}_annotatePeaks_Stats.txt

#awk 'BEGIN{k=0} {if(k==2 && $1=="Annotation"){exit;} if($1=="Annotation"){k+=2;} {print $0;}}' ./distribution_onlyGG/${bed}_elementStat_${con}.txt > ./distribution_onlyGG/${bed}_broadelement_${con}.txt
	awk 'BEGIN{k=0}
		{
			if($1=="Annotation"){
				k += 1;
			}
		
			if(k==2){
				print $0;
			}
		}
	' ${workdir}/distribution_onlyGG/${sample}_elementStat.txt > ${workdir}/distribution_onlyGG/${sample}_broadelement.txt
#python3 1make_table_for_plot_by_con.py -i ./distribution_onlyGG -b ${bed} -o ./distribution_onlyGG/${bed}_accCon_forplot.txt
#	Rscript 1plot_con.r -i ./distribution_onlyGG/${bed}_accCon_forplot.txt -o ./distribution_onlyGG/${bed}_accCon.pdf
done
##### divide subtypes
#mkdir -p ./distribution_onlyGG/subtypes/
#
#for con in hyper hypo
#do
##cat ${bedadd} | grep ${subtype} > ${bedadd}_${subtype}.bed
#	annotatePeaks.pl ${bedadd}__${con}.bed hg38 -annStats ./distribution_onlyGG/subtypes/${bed}_${subtype}_${con}_elementStat.txt > ./distribution_onlyGG/subtypes/${bed}_${subtype}_${con}_annotatePeaks_Stats.txt
#	awk 'BEGIN{k=0} {if(k==1 && $1=="Annotation"){exit;} if($1=="Annotation"){k+=1;} {print $0;}}' ./distribution_onlyGG/subtypes/${bed}_${subtype}_${con}_elementStat.txt > ./distribution_onlyGG/subtypes/${bed}_${subtype}_${con}_broadelement.txt
#done

/home/chengqiyi_pkuhpc/profiles/limingyang/miniforge3/bin/Rscript snk_script/7.1.plot_subtype_persample_onlyGG.r -i ${workdir}/distribution_onlyGG -o ${workdir}/distribution_onlyGG/annotation.pdf
