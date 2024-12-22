'''
Copyright (c) 2024-04-20 by LiMingyang, YiLab, Peking University.

Author: Li Mingyang (limingyang200101@gmail.com)

Institute: AAIS, Peking University

File Name: /gpfs1/chengqiyi_pkuhpc/limingyang/cisplatin/tools/downstream_pipeline/snk_script/16.pattern_signal.py

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
'''
from rich import print, pretty
from rich.traceback import install
pretty.install()
install(show_locals=True)
import rich_click as click
import polars as pl

@click.command()
@click.option('-t','--trt',help="Sample treatment")
@click.option('-c','--cut',help="Cut cutoff")
@click.option('-r','--rate',help="Rate cutoff")

def main(trt, cut, rate):
	ptrn = pl.read_csv(f'2.flank_region/cut{cut}rate{rate}/{trt}_pattern.bed', has_header=False, new_columns = [
		'chr','start','end','seq','ptrn','strand','trt','cut','rate','con'
	], separator='\t')
	ptrnfinal = pl.DataFrame()

	for con in ['ctrl','treat']:
		for strand in ['fwd','rvs']:
			if strand == 'fwd':
				strandmarker = "+"
			else:
				strandmarker = "-"
			
			dftmp = pl.read_csv(f'2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.{con}', has_header=False, new_columns = [
				'chr','pos','base','cov','stop','ratio','pvalue'
			], separator='\t')
			dftmp = dftmp.with_columns(
				(pl.col('pos')-2).alias('start'),
				pl.lit(strandmarker).alias('strand'),
				pl.lit(con).alias('con')
			).rename({
				"pos":"end"   
			}).drop([
				"base"
			])
			ptrnfinal = ptrnfinal.vstack(
				ptrn.join(dftmp, how='inner', on=['chr','start','end','strand','con'])
			)
	ptrnfinal.rename({
		"chr": "#chr"
	}).write_csv(f'2.flank_region/cut{cut}rate{rate}/{trt}_pattern_signal.bed', include_header=True, separator='\t')


if __name__ == "__main__":
	main()
