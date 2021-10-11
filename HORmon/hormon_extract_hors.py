############################################################################
# Copyright (c) 2021 Saint Petersburg State University
# Part of HORmon package. All Rights Reserved
# See file LICENSE for details.
############################################################################

from Bio.Seq import Seq
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

import sys
import os
from os import listdir
from os.path import isfile, isdir, join
import argparse
import copy
from pathlib import Path

MONOIDNT = 95

def parse_args():
    parser = argparse.ArgumentParser(description='HORmon HOR decomposition')
    parser.add_argument('-dec', '--decomposition', dest="monodecfile", help='tsv-file with SD monomer decomposition')
    parser.add_argument('-chor', '--canonical-hor', dest='horsfile', help='tsv-file with HOR from HORmon')
    parser.add_argument('-o', '--out-file', dest="outfilename", help="output file [default='./hor_decomposition.tsv']", default="./hor_decomposition.tsv",
                        required=False)
    return parser.parse_args()

def load_monodec(filename):
    dec = []
    monomers = set()
    rc_num = 0
    with open(filename, "r") as fin:
        for ln in fin.readlines():
           if len(ln.strip().split("\t")) < 5:
              continue
           ref, mon, start, end, idnt = ln.strip().split("\t")[:5]
           if float(idnt) > MONOIDNT:
               dec.append([ref, mon, start, end, idnt])
               monomers.add(mon)
               if mon.endswith("'"):
                   rc_num += 1
    if rc_num > 0.5*len(dec):
         revert = True
    else:
         revert = False
    return dec, monomers, revert

def shift(hor_lst):
    min_ind = 0
    for i in range(len(hor_lst)):
        if hor_lst[min_ind] > hor_lst[i]:
            min_ind = i
    return hor_lst[min_ind:] + hor_lst[:min_ind]

def build_cycle(hor_lst):
    hor = {}
    for i in range(len(hor_lst) - 1):
        hor[hor_lst[i]] = hor_lst[i + 1]
        hor[hor_lst[i + 1] + "'"] = hor_lst[i] + "'"
    hor[hor_lst[len(hor_lst) - 1]] = hor_lst[0]
    hor[hor_lst[0] + "'"] = hor_lst[len(hor_lst) - 1] + "'"
    return hor

def add_monomers(mono_mp, mono_cnt, hor_lst):
    for i in range(len(hor_lst)):
        mono_cnt += 1
        mono_mp[hor_lst[i]] = mono_cnt
        mono_mp[hor_lst[i]+"'"] = -mono_cnt
    mono_cnt += 1
    return mono_mp, mono_cnt

def load_horascycle(filename):
    hors, mono_mp = {}, {}
    mono_cnt, hor_cnt = 0, 0
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            if len(ln.split("\t")) < 2:
                continue
            hor_name, hor_seq = ln.strip().split("\t")[:2]
            hor_lst = shift(hor_seq.split(",")[:-1])
            hor = build_cycle(hor_lst)
            hors[hor_name] = [hor, len(hor_lst), hor_cnt]
            mono_mp, mono_cnt = add_monomers(mono_mp, mono_cnt, hor_lst)
            hor_cnt += 1
            print("Loaded HOR:", hor_lst)

    hors["Mono"] = ["Mono", 1, hor_cnt + 1]
    return hors, mono_mp

def decompose(monodec, hors, rev=False):
    hordec = []
    if len(monodec) == 0:
        return hordec

    fp = 0 if not rev else len(monodec) - 1

    inhor, cur_hor = 1, [monodec[fp]]
    cur_hor_name = "Mono"
    for h in hors:
        if monodec[fp][1] in hors[h][0]:
            cur_hor_name = h

    bgp, edp, stp = (1, len(monodec), 1) if not rev else (len(monodec) - 2, -1, -1)

    for i in range(bgp, edp, stp):
        prev, cur = (cur_hor[-1][1], monodec[i][1]) if stp == 1 else (monodec[i][1], cur_hor[0][1])
        mndist = abs(int(cur_hor[-1][3]) - int(monodec[i][2])) if stp == 1 else abs(int(monodec[i][3]) - int(cur_hor[0][2]))
        if float(monodec[i][4]) > MONOIDNT  and cur_hor_name in hors and prev in hors[cur_hor_name][0] and hors[cur_hor_name][0][prev] == cur and hors[cur_hor_name][1] > inhor and mndist  < 5:
            inhor += 1
            if stp == 1:
                cur_hor.append(monodec[i])
            else:
                cur_hor = [monodec[i]] + cur_hor
        else:
            idnt = sum([float(x[4]) for x in cur_hor])/len(cur_hor)
            hordec.append([cur_hor[0][0], ",".join([x[1] for x in cur_hor]), cur_hor[0][2], cur_hor[-1][3], "{:.2f}".format(idnt), cur_hor_name ])
            inhor = 1
            cur_hor = [monodec[i]]
            cur_hor_name = "Mono"
            for h in hors:
                if monodec[i][1] in hors[h][0]:
                    cur_hor_name = h
    idnt = sum([float(x[4]) for x in cur_hor])/len(cur_hor)
    hordec.append([cur_hor[0][0], ",".join([x[1] for x in cur_hor]), cur_hor[0][2], cur_hor[-1][3], "{:.2f}".format(idnt), cur_hor_name] )
    return hordec

def collapse_hor_name(mono_seq, mono_mp, hor):
    mono_lst = mono_seq.split(",")
    if len(mono_lst) == 1:
       return mono_lst[0]
       rev = "'" if mono_lst[0].endswith("'") else ""
       if mono_lst[0][1] in "0123456789X":
           return mono_lst[0][0] + rev
       else:
           return mono_lst[0][:2] + rev
    res = ""
    start, end = mono_mp[mono_lst[0]], mono_mp[mono_lst[0]]
    for i in range(1, len(mono_lst)):
        if mono_mp[mono_lst[i-1]] + 1 == mono_mp[mono_lst[i]]:
           end = mono_mp[mono_lst[i]]
        else:
           if end >= start or (start < 0 and end <= start):
              res += str(start) + "-" + str(end) + ","
           else:
              res += str(start) + ","
           start, end = mono_mp[mono_lst[i]], mono_mp[mono_lst[i]]
    if end >= start or (start < 0 and end <= start):
        res += str(start) + "-" + str(end)
    else:
        res += str(start)
    if len(mono_lst) == hor[1]:
        if start < 0:
            num = "-" + res[1:].split("-")[0]
        else:
            num = res.split("-")[0]
        return "c" + str(hor[2]) + "<sub>" + num + "</sub>"
    else:
        s, e = res.split("-")[0], res.split("-")[-1]
        if start < 0:
            s = "-" + res[1:].split("-")[0]
        if end < 0:
            e = "-" + res.split("-")[-1]
        return "p<sub>" + s + "-" + e + "</sub>"

def revert_hor(mon_seq):
    res = []
    for m in mon_seq.split(","):
      if m.endswith("'"):
        res.append(m[:-1])
      else:
        res.append(m+"'")
    return ",".join(res[::-1])

def collapse_hor(mono_seq, mono_mp, hor, revert):
    if revert:
      return collapse_hor_name(revert_hor(mono_seq), mono_mp, hor) + "'"
    else:
      return collapse_hor_name(mono_seq, mono_mp, hor)

def collapse_hordec(hordec, mono_mp, hors, revert):
    hordec_c = [hordec[0]]
    cnt, idnt = 1, float(hordec[0][4])
    for i in range(1, len(hordec)):
        if hordec[i][1] == hordec_c[-1][1]:
           cnt += 1
           idnt += float(hordec[i][4])
           hordec_c[-1][3] = hordec[i][3]
        else:
           cur_hor_name = hordec_c[-1][5]
           hordec_c[-1][1] = collapse_hor(hordec_c[-1][1], mono_mp, hors[cur_hor_name], revert)
           if cnt > 1:
               hordec_c[-1][1] += "<sup>" + str(cnt) + "</sup>"
           hordec_c[-1][4] = "{:.2f}".format(idnt/cnt)
           idnt, cnt = float(hordec[i][4]), 1
           hordec_c.append(hordec[i])
    hordec_c[-1][4] = "{:.2f}".format(idnt/cnt)
    cur_hor_name = hordec_c[-1][5]
    hordec_c[-1][1] = collapse_hor(hordec_c[-1][1], mono_mp, hors[cur_hor_name], revert)
    if cnt > 1:
        hordec_c[-1][1] += "<sup>" + str(cnt) + "</sup>"
    return hordec_c


def HORdecomposition(monodecfile, horsfile, outfile):
    monodec, monomers, revert = load_monodec(monodecfile)
    hors, mono_mp = load_horascycle(horsfile)
    outfilename = Path(outfile)
    hordec = decompose(monodec, hors)

    with open(outfilename, "w") as fout:
        for i in range(len(hordec)):
            if float(hordec[i][4]) > MONOIDNT:
                fout.write("\t".join(hordec[i][:-1]) + "\n")
    print("HOR decomposition saved to", outfilename)

    hordec_c = collapse_hordec(copy.deepcopy(hordec), mono_mp, hors, revert = False)
    outfilename_collapsedhors = outfilename.with_suffix(".collapsed.tsv")
    with open(outfilename_collapsedhors, "w") as fout:
        for i in range(len(hordec_c)):
           fout.write("\t".join(hordec_c[i]) + "\n")
    print("Collapsed HOR decomposition saved to", outfilename_collapsedhors)

    hordec_rc = hordec_c
    if revert:
      hordec_rc = collapse_hordec(hordec, mono_mp, hors, revert)
    outfilename_rcollapsedhors = outfilename.with_suffix(".reverted.collapsed.tsv")
    with open(outfilename_rcollapsedhors, "w") as fout:
        for i in range(len(hordec_c)):
           fout.write("\t".join(hordec_c[i]) + "\n")
    print("Collapsed HOR decomposition saved to", outfilename_rcollapsedhors)

def main():
    args = parse_args()
    HORdecomposition(args.monodecfile, args.horsfile, args.outfilename)


if __name__ == "__main__":
    main()
