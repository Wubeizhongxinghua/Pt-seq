#!/usr/bin/env python
# encoding: utf-8

import click
import re
def judge_gg(aline, ptrn, treatment, cut, rate, treat):
    posall = aline[0]
    seq = aline[1]
    thechr = posall.split(':')[0]
    thestart = int(re.findall(r"(?<=:)[0-9]+(?=-)", posall)[0]) + 9
    theend = int(re.findall(r"(?<=-)[0-9]+", posall)[0]) - 10
    thestrand = re.findall(r"(?<=\()(\+|-)", posall)[0]

    if seq[9].upper() == 'G' and seq[10].upper() == 'G': #is GG
       ptrn.write(f"{thechr}\t{thestart}\t{theend}\t{seq}\tGG\t{thestrand}\t{treatment}\t{cut}\t{rate}\t{treat}\n") 
    elif seq[9].upper() == 'A' and seq[10].upper() == 'G': #is AG 
       ptrn.write(f"{thechr}\t{thestart}\t{theend}\t{seq}\tAG\t{thestrand}\t{treatment}\t{cut}\t{rate}\t{treat}\n") 
    else:
       ptrn.write(f"{thechr}\t{thestart}\t{theend}\t{seq}\t{seq[9:11].upper()}\t{thestrand}\t{treatment}\t{cut}\t{rate}\t{treat}\n") 




@click.command()
@click.option('-t','--treatment',help='Treatment name')
@click.option('-c','--cut',help='cuterage')
@click.option('-r','--rate',help='Rate')

def main(treatment, cut, rate):
    with open(f'cut{cut}rate{rate}/{treatment}_fwd_stop_articut_{cut}_{rate}.txt.ctrl_expand_flank10_cisplatin_site.txt') as trt_fwd_ctrl, \
    open(f'cut{cut}rate{rate}/{treatment}_rvs_stop_articut_{cut}_{rate}.txt.ctrl_expand_flank10_cisplatin_site.txt') as trt_rvs_ctrl, \
    open(f'cut{cut}rate{rate}/{treatment}_fwd_stop_articut_{cut}_{rate}.txt.treat_expand_flank10_cisplatin_site.txt') as trt_fwd_trt, \
    open(f'cut{cut}rate{rate}/{treatment}_rvs_stop_articut_{cut}_{rate}.txt.treat_expand_flank10_cisplatin_site.txt') as trt_rvs_trt, \
    open(f'cut{cut}rate{rate}/{treatment}_pattern.bed', 'w') as ptrn:

        for aline in trt_fwd_ctrl:
            theline = aline.strip().split('\t') 
            judge_gg(theline, ptrn, treatment, cut, rate, 'ctrl')

        for aline in trt_rvs_ctrl:
            theline = aline.strip().split('\t') 
            judge_gg(theline, ptrn, treatment, cut, rate, 'ctrl')

        for aline in trt_fwd_trt:
            theline = aline.strip().split('\t') 
            judge_gg(theline, ptrn, treatment, cut, rate, 'treat')

        for aline in trt_rvs_trt:
            theline = aline.strip().split('\t') 
            judge_gg(theline, ptrn, treatment, cut, rate, 'treat')

if __name__ == '__main__':
    main()
