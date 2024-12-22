import argparse
import sys
import subprocess
from Bio import SeqIO
import pysam
import re
import time
import gzip

def filtered_reads(read1,read2,outname_prx):
    final_1 = outname_prx + "_rmSbfI_1.fq.gz"
    final_2 = outname_prx + "_rmSbfI_2.fq.gz"
    subprocess.call("rm -f " + final_1, shell=True)
    subprocess.call("rm -f " + final_2, shell=True)
    # print('xxx')
    with pysam.FastxFile(read1) as fh1, pysam.FastxFile(read2) as fh2, gzip.open(final_1,  mode='a+') as final1, \
            gzip.open(final_2, mode='a+') as final2:
        for (entry1,entry2) in zip(fh1,fh2):
            # print(entry1,"************1")
            seq1=entry1.sequence
            seq2=entry2.sequence
            # print(seq1,seq2)
            # time.sleep(1000)
            if not re.search('(T)*CCTGCAGGGCTC', seq1):
                final1.write(str(entry1).encode() + b"\n")
                final2.write(str(entry2).encode() + b"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="get mutation rates for m6A_new")
    parser.add_argument("-read1", "--read1",nargs="?", type=str, default=sys.stdin, help="read1")
    parser.add_argument("-read2", "--read2",nargs="?", type=str, default=sys.stdin, help="read2")
    parser.add_argument("-outname_prx", "--outname_prx", nargs="?", type=str, default=sys.stdin, help="outname_prx")

    options = parser.parse_args()
    print("lalalalal")

    filtered_reads(options.read1,options.read2,options.outname_prx)
