'''
Copyright (c) 2024-04-11 by LiMingyang, YiLab, Peking University.

Author: Li Mingyang (limingyang200101@gmail.com)

Institute: AAIS, Peking University

File Name: /gpfs1/chengqiyi_pkuhpc/limingyang/cisplatin/tools/downstream_pipeline/snk_script/0.0.cal_sig_bg.py

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
import polars as pl
import rich_click as click

@click.command()
@click.option('-if','--inputsig_fwd', help="Input sig")
@click.option('-ir','--inputsig_rvs', help="Input sig")
@click.option('-o','--outputv', help="Input sig")
def main(inputsig_fwd, inputsig_rvs, outputv):
	dffwd = pl.read_csv(inputsig_fwd, has_header=False, separator='\t', new_columns=['chr','pos','base','cov','stop'])
	dfrvs = pl.read_csv(inputsig_rvs, has_header=False, separator='\t', new_columns=['chr','pos','base','cov','stop'])
	df = dffwd.vstack(dfrvs)
	bg_ratio = df['stop'].sum()/df['cov'].sum()
	with open(outputv, 'w') as f:
		f.write(f'{bg_ratio}\n')

if __name__ == "__main__":
	main()
