from Bio.Seq import Seq
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

import sys
import os
from os import listdir
from os.path import isfile, isdir, join

import re
import edlib
import random

def conver2bed(tsvfile):
    outfile = tsvfile[:-len(".tsv")] + ".bed"
    mono = {}
    with open(outfile, "w") as fout:
        with open(tsvfile, "r") as fin:
            for ln in fin.readlines():
                lst = ln.strip().split("\t")
                tig_name = lst[0]  # .split(":")[0]
                # tstart, tend = lst[0].split(":")[1].split("-")
                tstart, tend = 0, 0  # int(tstart), int(tend)
                rev = "+"
                if "'" in lst[1]:
                    lst[1] = lst[1].replace("'", "")
                    rev = "-"
                if lst[1] in mono:
                    r, g, b = mono[lst[1]]
                else:
                    r, g, b = random.randint(0, 255), random.randint(0, 255), random.randint(0, 255)
                    mono[lst[1]] = [r, g, b]
                fout.write("\t".join([tig_name, str(tstart + int(lst[2])), str(tstart + int(lst[3])) \
                                         , lst[1], str(int(float(lst[4]))), rev, str(tstart + int(lst[2])),
                                      str(tstart + int(lst[3])) \
                                         , str(r) + "," + str(g) + "," + str(b)]) + "\n")
    print("Please find result in ", outfile)
    return outfile


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please specify tsv-file to convert: python3 convert2bed.py <filename>")
        exit(-1)
    conver2bed(sys.argv[1])