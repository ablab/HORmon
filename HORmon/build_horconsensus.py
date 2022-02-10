############################################################################
# Copyright (c) 2021 Saint Petersburg State University
# Part of HORmon package. All Rights Reserved
# See file LICENSE for details.
############################################################################

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

from Bio import AlignIO

import sys
import os
from os import listdir
from os.path import isfile, isdir, join
import argparse
import pandas as pd
import numpy as np

import subprocess

import edlib
import random

import HORmon.HORmon_pipeline.utils as utils

MONOIDNT = 95
INF=10000000
sd_run = "stringdecomposer"

def load_fasta(filename, tp = "list"):
    if tp == "map":
        records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
        for r in records:
            records[r] = records[r].upper()
    else:
        records = list(SeqIO.parse(filename, "fasta"))
        for i in range(len(records)):
            records[i] = records[i].upper()
    return records

def make_record(seq, name, sid, d=""):
    return SeqRecord(seq, id=sid, name=name, description = d)

def save_fasta(filename, orfs):
    with open(filename, "w") as output_handle:
        SeqIO.write(orfs, output_handle, "fasta")
        #fasta_out = FastaIO.FastaWriter(output_handle)
        #fasta_out.write_file(orfs)

def load_monodec(cen):
    dec = []
    monomers = set()
    filename = os.path.join(PATH, "cen" + cen + "_dec.tsv")
    with open(filename, "r") as fin:
         for ln in fin.readlines():
             if len(ln.strip().split("\t")) < 5:
                 continue
             ref, mon, start, end, idnt = ln.strip().split("\t")[:5]
             dec.append([ref, mon, start, end, idnt])
             monomers.add(mon)
    #dec = revert_rcreads(dec)
    #dec = remove_badreads(dec)
    return dec, monomers

def load_bedfile(bedfile):
    dec = []
    monomers = set()
    with open(bedfile, "r") as fin:
         for ln in fin.readlines():
             if len(ln.strip().split("\t")) < 6:
                 continue
             ref, start, end, mon, idnt, rev = ln.strip().split("\t")[:6]
             dec.append([ref, mon, start, end, idnt, rev])
             monomers.add(mon)
    return dec, monomers

def shift(hor_lst):
    min_ind = 0
    for i in range(len(hor_lst)):
        if hor_lst[min_ind] > hor_lst[i]:
            min_ind = i
            break
    return hor_lst[min_ind:] + hor_lst[:min_ind]

def load_horascycle(filename, log=None):
    if log is not None:
        log.info("\n= Loading HORs cycles =", indent=1)
    hors = []
    hor_name= ""
    mono_mp = {}
    cnt, hor_cnt = 0, 0
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            if len(ln.split("\t")) < 2:
                continue
            hor_name, hor_seq = ln.strip().split("\t")[:2]
            hor_lst = hor_seq.split(",")[:-1]
            #hor_lst = shift(hor_lst)
            if log is None:
                print(hor_lst)
            else:
                log.info("Loaded HOR: " + str(hor_lst), indent=2)
            hors.append([hor_name, hor_lst])
    return hors


def run_clustal(mappings, clustal_dir, pair_name):
    from Bio.Align.Applications import ClustalwCommandline
    from Bio.Align.Applications import ClustalOmegaCommandline
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from Bio.Align import MultipleSeqAlignment

    cluster_seqs_path = os.path.join(clustal_dir, pair_name + "_seq.fasta")
    aln_file = os.path.join(clustal_dir, pair_name + "_seq.clu")
    if not os.path.isdir(clustal_dir):
       os.makedirs(clustal_dir)
    if len(mappings) == 1:
        save_fasta(cluster_seqs_path, mappings + mappings)
    else:
        save_fasta(cluster_seqs_path, mappings)

    cmd = ClustalOmegaCommandline(infile=cluster_seqs_path, outfile=aln_file, force=True, threads=10)
    stdout, stderr = cmd()
    align = AlignIO.read(aln_file, "fasta")

    summary_align = AlignInfo.SummaryInfo(align)
    consensus = summary_align.gap_consensus(threshold=0, ambiguous='N')
    consensus = str(consensus).replace('-', '')
    return consensus

def extract_consensus(total_alns):
    consensus = ""
    scores = []
    for i in range(len(total_alns[0])):
        score = {"A": 0, "C": 0, "G": 0, "T": 0, "-": 0}
        for j in range(len(total_alns)):
            score[total_alns[j][i]] += 1
        scores_lst = sorted([[it, score[it]] for it in score ], key = lambda x: -x[1])
        scores.append(scores_lst)
        max_score = scores_lst[0][0]
        if scores_lst[0][1] == scores_lst[1][1]:
            max_score = "N"
        if max_score != "-":
            consensus += max_score
    return consensus, scores

def align_mappings(mappings, clustal_dir, pair_name):
    pair_name = pair_name.replace("/", "_")
    pair_name = pair_name.replace("(", "_")
    pair_name = pair_name.replace(")", "_")
    pair_name = pair_name.replace(".", "_")
    pair_name = pair_name.replace("&", "_")
    consensus = run_clustal(mappings, clustal_dir, pair_name)
    return consensus

def edist(lst):
    if len(str(lst[0])) == 0:
        return INF, []
    if len(str(lst[1])) == 0:
        return INF, []
    result = edlib.align(str(lst[0]), str(lst[1]), mode="SHW", task="path", k=100)
    if result["editDistance"] == -1:
        return INF, []
    aln = edlib.getNiceAlignment(result, str(lst[0]), str(lst[1]))
    return result["editDistance"], aln

def glue_pairs(p1, p2, log=None):
    max_len = 200
    eds = []
    ed, aln = edist([p1[-max_len:], p2])
    longest, longest_ind = 0, -1
    cur_len = 0
    for i in range(len(aln["matched_aligned"])):
        if aln["matched_aligned"][i] == "|":
            cur_len += 1
        else:
            if cur_len > longest:
                longest, longest_ind = cur_len, i - cur_len
            cur_len = 0
    if cur_len > longest:
        longest, longest_ind = cur_len, len(aln["matched_aligned"]) - cur_len
    i, j = len(p1) - max_len + len(aln["query_aligned"][:longest_ind].replace("-", "")), len(aln["target_aligned"][:longest_ind].replace("-", ""))

    if log is None:
        print(i, j, ed, len(p1[:i] + p2[j:]))
        print("")
    else:
        log.info("i=" + str(i) + "; j=" + str(j) +"; edit dist=" + str(ed) + "; pair union len=" + str(len(p1[:i] + p2[j:])), indent=3)
    return i, j

def build_monoconsensus(monodec, ref, step, clustal_dir, log=None):
    mappings = {}
    for i in range(len(monodec)):
        m_name, m_s, m_e, m_idnt = monodec[i][1], int(monodec[i][2]), int(monodec[i][3]), float(monodec[i][4])
        if m_idnt > MONOIDNT:
            if m_name not in mappings:
                mappings[m_name] = []
            if monodec[i][5] == "+":
                mappings[m_name].append(make_record(ref[monodec[i][0]].seq[max(0, m_s - step): min(m_e + step, len(ref[monodec[i][0]].seq) )], m_name + str(m_s), m_name + str(m_s)))
            else:
                mappings[m_name].append(make_record(ref[monodec[i][0]].seq[max(0, m_s - step): min(m_e + step, len(ref[monodec[i][0]].seq) )].reverse_complement(), m_name + str(m_s) + "_rev", m_name + str(m_s) + "_rev"))
    m_consensus = []
    for m in mappings:
        if log != None:
            log.info("Detecting monoconsensus for " + m, indent=1)
        else:
            print(m)
        consensus = align_mappings(mappings[m], clustal_dir, m)
        consensus = consensus[step : len(consensus)-step]
        cur_consensus = make_record(Seq(consensus), m, m)
        m_consensus.append(cur_consensus)
    return m_consensus

def build_pairconsensus(hors, monodec, ref, clustal_dir, log=None):
    if log is not None:
        log.info("\n= Build Pair Consensus =", indent=1)

    hors_seq = []
    for hor_desc in hors:
        hor = hor_desc[1]
        pairs_consensus = []
        for i in range(len(hor)):
            m1, m2 = hor[i], hor[(i+1)%len(hor)]
            mappings = []
            for j in range(len(monodec)-1):
                if ((monodec[j][1] == m1 and monodec[j+1][1] == m2 and monodec[j][5] == monodec[j+1][5] == "+") or\
                    (monodec[j][1] == m2 and monodec[j+1][1] == m1 and monodec[j][5] == monodec[j+1][5] == "-"))  \
                    and float(monodec[j][4]) > MONOIDNT and float(monodec[j+1][4]) > MONOIDNT:
                    s1, e1, s2, e2 = int(monodec[j][2]), int(monodec[j][3]), int(monodec[j+1][2]), int(monodec[j+1][3])
                    if s2 - e1 < 10:
                        if monodec[j][5] == "+":
                            mappings.append(make_record(ref[monodec[j][0]].seq[s1: e2], m1+"_" + m2 + str(s1), m1 + "_" +m2 + str(s1) ))
                        else:
                            mappings.append(make_record(ref[monodec[j][0]].seq[s1: e2].reverse_complement(), m1+"_" + m2 + str(s1) + "_rev", m1 + "_" +m2 + str(s1) + "_rev" ))

            if log is None:
                print("pair: ", m1, m2, len(mappings))
            else:
                log.info("# pairs(" + m1 + ", " + m2 + ") = " + str(len(mappings)), indent=2)

            if len(mappings) > 0:
                pairs_consensus.append(align_mappings(mappings, clustal_dir, m1+"_" + m2))

        if len(pairs_consensus) > 0:
            cur_consensus = ""
            border = []
            for i in range(len(pairs_consensus)):
                if log is None:
                    print("Pair ", i)
                else:
                    log.info("Handle pair #" + str(i), indent=2)
                border.append(glue_pairs(pairs_consensus[i], pairs_consensus[(i+1)%len(pairs_consensus)], log))

            l, r = border[len(pairs_consensus)-1][1], 0
            for i in range(len(pairs_consensus)):
                r = border[i][0]
                if l > r:
                   if log is None:
                       print("Something went wrong!", i, l, r)
                       exit(-1)
                   else:
                       log.error("Something went wrong! " + str(i) + " " + str(l) + " " + str(r))

                cur_consensus += pairs_consensus[i][l: r]
                l = border[i][1]

            if log is None:
                print(len(cur_consensus))
            else:
                log.info("\n")
                log.info("HOR consensus len=" + str(len(cur_consensus)), indent=2)

            hors_seq.append(make_record(Seq(cur_consensus), hor_desc[0], hor_desc[0], str(len(cur_consensus)) + "bp " + ",".join(hor) ))
        else:
            hors_seq.append(make_record(Seq(""), hor_desc[0], hor_desc[0], str(0) + "bp " + ",".join(hor) ))
    return hors_seq


def run(sequences, monomers, num_threads, out_file, log=None):
    if log is not None:
        log.info("Run StringDecomposer", indent=2)
    utils.sys_call([sd_run, sequences, monomers, "-t", str(num_threads), "--out-file", out_file])
    with open(out_file + ".tsv", 'r') as f:
        out_decomposition = "".join(f.readlines())
    return out_decomposition


def divide_into_monomers(hors_lst, hors, monomers, horfile, monofile, outtsv, log=None):
    if log is not None:
        log.info("\n= Divide into monomers =", indent=1)
    hors_res, res, bed = [], [], []
    for i in range(len(hors)):
        h, start_mono, hor_len = hors[i], hors_lst[i][1][0], len(hors_lst[i][1])
        if len(h.seq) > 0:
          triple_hors = []
          triple_hors.append(make_record(h.seq + h.seq + h.seq, h.id, h.name, h.description))
          save_fasta(horfile, triple_hors)
          decomposition = run(horfile, monofile, "1", outtsv[:-len(".tsv")], log)
          inHOR, start_shift = False, 0
          newhor_seq = ""

          if log is None:
              print(start_mono, hor_len, len(h.seq))
          else:
              log.info("Start monomer=" + str(start_mono) + ", HOR len(#mon)=" + str(hor_len) + ", HOR len(bp)=" + str(len(h.seq)), indent=2)
          for ln in decomposition.split("\n")[1:-1]:
              ref, mono, start, end, score = ln.split("\t")[:5]
              if mono == start_mono:
                  inHOR, start_shift = True, int(start)
                  if log is None:
                      print("shift", start_shift)
                  else:
                      log.info("Start monomer shift=" + str(start_shift), indent=2)

              if inHOR and len(res) < hor_len:
                  if log is None:
                      print(mono, ref + ":" + str(int(start) - start_shift) + "-" + str(int(end) + 1 - start_shift), len(res) )
                  else:
                      log.info("Monomer# " + str(len(res)) + ": " + mono + " " + ref + ":" + str(int(start) - start_shift) + "-" + str(int(end) + 1 - start_shift), indent=2)

                  res.append(make_record(triple_hors[0].seq[int(start): int(end) + 1], mono, mono, ref + ":" + str(int(start) - start_shift) + "-" + str(int(end) + 1 - start_shift) ))
                  r, g, b = random.randint(0, 255), random.randint(0, 255), random.randint(0, 255)
                  bed.append("\t".join([ref, str(int(start) - start_shift), str(int(end) + 1 - start_shift), mono, str(int(float(score))), "+", str(int(start) - start_shift), str(int(end) + 1 - start_shift), ",".join([str(r), str(g), str(b)]) ]))
                  newhor_seq += triple_hors[0].seq[int(start): int(end) + 1]
          hors_res.append(make_record(newhor_seq, h.id, h.id, str(len(newhor_seq)) + "bp " + ",".join(hors_lst[i][1]) ))
    return hors_res, res, bed


def getHORconsensus(hors_tsv, monodec, ref, consensus, mono_outfilename, outdir, log=None):
    hors = load_horascycle(hors_tsv, log)
    hor_consensus = build_pairconsensus(hors, monodec, ref, os.path.join(outdir, "clustal_alns"), log=log)

    if len(hor_consensus) > 0:
        hor_outfilename = os.path.join(outdir, "hor_consensus.fasta")
        hor_consensus_shifted, pair_monomers, bed = divide_into_monomers(hors, hor_consensus, consensus,
                                                                         hor_outfilename, mono_outfilename,
                                                                         os.path.join(outdir, "sd.tsv"), log)
        hor_outfilename = os.path.join(outdir, "hor_consensus.fasta")
        save_fasta(hor_outfilename, hor_consensus_shifted)

        if log is not None:
            log.info("HOR consensus can be found" + hor_outfilename, indent=1)
        else:
            print("HOR consensus can be found", hor_outfilename)

        pmono_outfilename = os.path.join(outdir, "monomer_paired_consensus.fasta")
        save_fasta(pmono_outfilename, pair_monomers)

        if log is not None:
            log.info("Paired monomer consensus can be found " + pmono_outfilename, indent=1)
        else:
            print("Paired monomer consensus can be found", pmono_outfilename)


        outfilename = os.path.join(outdir, "monomer_paired_consensus.bed")
        with open(outfilename, "w") as fout:
            for ln in bed:
                fout.write(ln + "\n")

        if log is not None:
            log.info("Paired monomer consensus bed can be found " + outfilename, indent=1)
        else:
            print("Paired monomer consensus bed can be found", outfilename)
    else:
        print("Paired wasn't generated - no pairs found")


def build_horcons(sepath, bedfile, odir, horspath, extend=4, seed=123, log=None):
    random.seed(int(seed))
    ref = load_fasta(sepath, "map")

    if not os.path.exists(odir):
        os.makedirs(odir)

    mono_outfilename = os.path.join(odir, "monomer_consensus.fasta")
    monodec, monomers = load_bedfile(bedfile)
    consensus = build_monoconsensus(monodec, ref, int(extend), os.path.join(odir, "clustal_alns"), log=log)
    save_fasta(mono_outfilename, consensus)

    if log is not None:
        log.info("\nMonomer consensus can be found " + mono_outfilename, indent=1)
    else:
        print("Monomer consensus can be found ", mono_outfilename)

    if horspath != None:
        getHORconsensus(horspath, monodec, ref, consensus, mono_outfilename, odir, log=log)


def main():
    parser = argparse.ArgumentParser(description='Extracts monomer/HOR consensus from annotation')
    parser.add_argument('sequences', help='fasta-file with annotated sequences')
    parser.add_argument('annotation', help='bed-file with annotation')
    parser.add_argument('outdir', help='output directory')
    parser.add_argument('--hors',  help='tsv-file with HOR description', required=False)
    parser.add_argument('--extend',  help='number of bp to extend monomer alignment', default=4, required=False)
    parser.add_argument('--seed',  help='seed for colors generation in bed-files', default=123, required=False)
    args = parser.parse_args()

    build_horcons(args.sequences, args.annotation, args.outdir, args.hors, args.extend, args.seed)


if __name__ == "__main__":
    main()
