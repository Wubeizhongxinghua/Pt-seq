import argparse
import sys
import time
import subprocess
import re
from scipy.stats import binomtest
from multiprocessing import Pool

parser = argparse.ArgumentParser(description="parse the mpileup file")
parser.add_argument("-ctrl","--controlmpi",nargs="?",type=str,default=sys.stdin,help="controlbam")
parser.add_argument("-treat","--treatmpi",nargs="?",type=str,default=sys.stdin,help="bamfile")
parser.add_argument("-O","--output",nargs="?",type=str,default=sys.stdout,help="Output file prefix")
parser.add_argument("-stop_reads_cut","--stop_reads_cut",type=int,default=15, help="stop reads number greater than X")
parser.add_argument("-stop_rate_cut","--stop_rate_cut",type=int,default=20,help="stop rate greater than X")
parser.add_argument("-cov_cut","--cov_cut",type=int,default=15,help="stop rate greater than X")
parser.add_argument("-p_cut","--pvalue_cut",type=float,default=0.05,help="stop site with p-value of binominal test greater or equal than X (not suitable for ctrl sample)")
parser.add_argument("-proba","--proba_binom",type=float,default=0.05,help="probability of signal-noise ratio used for binom_test")

args = parser.parse_args()
#
#controlbam = args.controlmpi
#treatbam = args.treatmpi
#Output = args.output
#cov_cut = args.cov_cut
#stop_rate_cut = args.stop_rate_cut
#stop_reads_cut = args.stop_reads_cut
#p_cut = args.pvalue_cut
#proba = args.proba_binom
#
#t1 = time.time()

def stop_site(args):
	pileup, stop_reads_cut, stop_rate_cut, cov_cut, output, pcut, proba = args
	with open(output, 'w+') as out1:
		for line in open(pileup):
			chr_name, pos, base, cov, reads = line.strip().split("\t")
			cov = int(cov)
			reads = int(reads)
			if cov > cov_cut:
				stop_rate = (reads / float(cov)) * 100
				pvalue = binomtest(k=reads, n=cov, p=proba, alternative='greater').pvalue
				if reads > stop_reads_cut and stop_rate > stop_rate_cut and pvalue <= pcut:
					if re.search('forward', pileup) or re.search('fwd', pileup):
						pos = int(pos) - 1  # fwd zero-based
					else:
						pos = int(pos) + 1  # rvs zero-based
					out_list = [chr_name, pos, base, cov, reads, stop_rate, pvalue]
					out1.write("\t".join(map(str, out_list)) + '\n')

if __name__ == "__main__":
	t1 = time.time()
	pool = Pool(processes=2)  # Number of processes to create
	tasks = [
		(args.treatmpi, args.stop_reads_cut, args.stop_rate_cut, args.cov_cut, f'{args.output}.treat', args.pvalue_cut, args.proba_binom),
		(args.controlmpi, args.stop_reads_cut, args.stop_rate_cut, args.cov_cut, f'{args.output}.ctrl', 1, args.proba_binom) #ctrl不进行统计检验，保留尽可能多的背景信号进行剔除
		#(treatbam, stop_reads_cut, stop_rate_cut, cov_cut, f'{Output}.treat', p_cut, proba),
		#(controlbam, stop_reads_cut, stop_rate_cut, cov_cut, f'{Output}.ctrl', 1, proba)
	]
	pool.map(stop_site, tasks)
	pool.close()
	pool.join()
	print("Time taken:", time.time() - t1)







#def stop_site(pileup, stop_reads_cut, stop_rate_cut, cov_cut,output):
#	 print("**************2",pileup)
#	out1 = open(output,'w+')
#	first_list = []
#	if re.search('forward', pileup) or re.search('fwd', pileup):
#		 print("**************3", pileup)
#		for lines in open(pileup):
#			line = lines.strip().split("\t")
#			cov = int(line[3])
#			if cov > cov_cut:
#				print(line)
#				mp = str(line[4].upper())
#				stop = int(mp.count("^"))
#				stop = int(line[4])
#				stop_rate = (stop / float(cov)) * 100
#				pvalue = binomtest(k=stop, n=cov, p=proba, alternative='greater')
#				print(stop,stop_reads_cut,stop_rate_cut,line)
#				if stop > stop_reads_cut and stop_rate > stop_rate_cut and pvalue >= p_cut:
#					print(line)
#					chr_name = line[0]
#					base = line[2]
#					pos = int(line[1]) - 1 zero-based
#					out_list = [chr_name, pos, base, cov, stop, stop_rate, pvalue]
#					out1.write("\t".join(map(str, out_list))+'\n')
#						first_list.append(out_list)
#	else:
#		print("***************")
#		line_num = 0
#
#		for lines in open(pileup):
#			line_num += 1
#			line = lines.strip().split("\t")
#			try:
#				cov = int(line[3])
#			except:
#				print(line, line_num)
#			if cov > cov_cut:
#				mp = str(line[4].upper())
#				stop = int(mp.count("$"))
#				stop_rate = (stop / float(cov)) * 100
#				if stop > stop_reads_cut and stop_rate > stop_rate_cut:
#					chr_name = line[0]
#					base = line[2]
#					pos = int(line[1]) + 1
#					out_list = [chr_name, pos, base, cov, stop,stop_rate]
#					out1.write("\t".join(map(str, out_list))+'\n')
#
#	first_list2 = ["\t".join(map(str,i)) for i in first_list]
#	out1.writelines("\n".join(first_list2)+"\n")
#	return first_list
#
#
##def rm_ctrl(treat_list,ctrl_list,output):
##	print(len(treat_list),len(ctrl_list))
##	list_ctrl = []
##	list_treat = []
##	for i in ctrl_list:
##		p1 = "_".join(map(str,[i[0],i[1]]))
##		list_ctrl.append(p1)
##	print(len(list_ctrl))
##	for j in treat_list:
##		p2 = "_".join(map(str,[j[0],j[1]]))
##		if p2 not in list_ctrl:
##			# print(p1,p2)
##			list_treat.append(j)
##	print(len(list_treat))
##	first_list2 = ["\t".join(map(str,i)) for i in list_treat]
##	out1 = open(output, 'w+')
##	print("lalalala")
##	out1.writelines("\n".join(first_list2) + "\n")
#
#
#if __name__ == "__main__":
#	print("**************1")
#	stop_site(treatbam, stop_reads_cut, stop_rate_cut, cov_cut, Output+'.treat')
#	stop_site(controlbam, stop_reads_cut, stop_rate_cut, cov_cut, Output+'.ctrl')
#
#	# treat_list = [i.strip().split("\t") for i in open(Output+'.treat').readlines()]
#	# ctrl_list = [i.strip().split("\t") for i in open(Output+'.ctrl').readlines()]
#	#rm_ctrl(treat_list, ctrl_list, Output+'.treat.final')

