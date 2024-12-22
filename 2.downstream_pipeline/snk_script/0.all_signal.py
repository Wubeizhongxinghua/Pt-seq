'''
Copyright (c) 2024-04-10 by LiMingyang, YiLab, Peking University.

Author: Li Mingyang (limingyang200101@gmail.com)

Institute: AAIS, Peking University

File Name: /gpfs1/chengqiyi_pkuhpc/limingyang/cisplatin/output_q1to10_unique/batch7/align/genome/0.all_signal.py

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
import multiprocessing
from functools import partial
import rich_click as click

def process_chunk(chunk, strand):
	# 处理每个块的函数
	result = []
	for line in chunk:
		parts = line.strip().split('\t')
		if len(parts) < 5:
			raise ValueError(f"A line not complete: {parts}")
			#continue  # 跳过不完整的行
		if strand in ['fwd','+']:
			count_char = parts[4].count('^')  # 计算'^'的数量
		elif strand in ['rvs','-']:
			count_char = parts[4].count('$')  # 计算'$'的数量
		else:
			count_char = 0  # 如果没有正确的strand值，返回0
		# 构造新的行输出
		if count_char == 0:
			continue #跳过没有信号的位点
		new_line = '\t'.join(parts[:4] + [str(count_char)])
		result.append(new_line)
	return result

def read_in_chunks(file, chunk_size=1000):
	# 文件块读取生成器
	chunk = []
	for line in file:
		chunk.append(line)
		if len(chunk) >= chunk_size:
			yield chunk
			chunk = []
	if chunk:
		yield chunk



@click.command()
@click.option('-i','--filename',help='Input Mpileup')
@click.option('-o','--output_filename',help='Output bg and signal')
@click.option('-s','--strand',type=click.Choice(['fwd','+', 'rvs','-'], case_sensitive=False), help="Strand")
@click.option('-p','--threads', default=None, help='Number of threads.')
def main(filename, output_filename, strand, threads):
	# 主并行处理函数
	if threads is None:
		num_processes = multiprocessing.cpu_count()
	else:
		num_processes = int(threads)


	pool = multiprocessing.Pool(processes=num_processes)
	tasks = []

	with open(filename, 'r') as file:
		chunk_generator = read_in_chunks(file)
		for chunk in chunk_generator:
			task = pool.apply_async(process_chunk, (chunk, strand))
			tasks.append(task)

	# 等待所有进程处理完成，并收集结果
	results = []
	for task in tasks:
		results.extend(task.get())

	pool.close()
	pool.join()

	# 将结果写入到输出文件
	with open(output_filename, 'w') as f:
		for line in results:
			f.write(line + '\n')

if __name__ == "__main__":
	main()
