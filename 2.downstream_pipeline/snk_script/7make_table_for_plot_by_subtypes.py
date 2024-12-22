import pandas as pd
from pathlib import Path
import click

@click.command()
@click.option('-i','--input',help='Input broad element file dir without "/"')
@click.option('-o','--output',help='Output file for plot')
#@click.option('-o','--output',help='Outputfile name (no dir!!, the dir is the same as input)')

def main(input,output):


	elements = ['3UTR', 'Retroposon', 'RC?', 'RNA', 'miRNA', 'ncRNA', 'TTS', 'LINE', 'srpRNA', 'SINE', 'RC', 'tRNA', 'DNA?', 'pseudo', 'DNA', 'Exon', 'Intron', 'Intergenic', 'Promoter', '5UTR', 'snoRNA', 'LTR?', 'scRNA', 'CpG-Island', 'Low_complexity', 'LTR', 'Simple_repeat', 'snRNA', 'Unknown', 'SINE?', 'Satellite', 'rRNA']
	df = pd.DataFrame(columns=['Sample']+elements)
	i=0

	directory = Path(input)

	filenames = []
	samples = []

	for file in directory.glob('*_broadelement.txt'):
		filenames.append(file.name)
		sample = file.stem.rsplit('_broadelement', 1)[0]  # Extracting the prefix before '_broadelement'
		samples.append(sample)

	for sample in samples:
		dff = pd.read_table(f"{input}/{sample}_broadelement.txt",header=0)
		values = [sample]
		for element in elements:
			values.append(float(dff[dff['Annotation']==element][sample][0]))
		df.loc[i] = values
		i = i + 1

	df.to_csv(output,header=True,index=False,sep='\t')

if __name__=='__main__':
	main()
