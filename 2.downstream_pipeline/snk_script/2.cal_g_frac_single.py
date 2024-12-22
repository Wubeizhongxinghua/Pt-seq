import click
import pandas as pd
import numpy as np
import re
from tqdm import tqdm
@click.command()
@click.option('-i','--inputfile',help='Sequence file list')
@click.option('-d','--workdir',help='Workdir')



def main(inputfile, workdir):
	def countg_fun(seq):
		count = np.zeros(21)
		for char, i in zip(seq.upper(), range(len(seq))):
			if char == 'G':
				count[i] += 1
		return count
	dffracall = pd.DataFrame(columns = ['sample','pos','fraction'])
	with open(f'{workdir}/{inputfile}') as il:
		inputfile = inputfile.strip()
		if len(re.findall('treat', inputfile)):
			con = 'treat_'
		else:
			con = 'ctrl_'

		samplename = '_'.join(inputfile.split('_')[0:2])
		samplename = con+samplename
		colseq = list(range(-10, 11))
		zero21 = np.zeros((1,21))
		dfg = pd.DataFrame(data = zero21, columns = colseq)
		dfall = pd.DataFrame(data = zero21, columns = colseq)
		with open(f'{workdir}/{inputfile}', 'r') as f:
			for line in tqdm(f, leave=False):
				seq = line.strip().split('\t')[1]
				countg = countg_fun(seq)
				dfg += countg
				dfall += np.ones(21)
		dffrac = dfg / dfall
		dffrac['sample'] = samplename
		dffrac = dffrac.melt(id_vars='sample', value_vars=colseq).rename(columns={'variable':'pos', 'value':'fraction'})
		dffracall = pd.concat([dffracall, dffrac])
	dffracall.to_csv(f'{workdir}/G_frac_{inputfile}.txt', header=True, index=False, sep='\t')

if __name__=='__main__':
	main()
