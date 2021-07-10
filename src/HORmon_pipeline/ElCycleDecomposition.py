#!/usr/bin/env python3

import os
import networkx as nx
from networkx.algorithms import bipartite
from networkx.drawing.nx_agraph import write_dot
from subprocess import check_call
import math
import os
import pandas as pd
from Bio import SeqIO

import HORmon_pipeline.DetectHOR as DetectHOR
import HORmon_pipeline.MergeAndSplitMonomers as splitMn
from HORmon_pipeline.utils import rc
import HORmon_pipeline.utils as utils
from HORmon_pipeline.utils import run_SD
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



def save_seqs(blocks, cluster_seqs_path):
    with open(cluster_seqs_path, "w") as fa:
        for i in range(len(blocks)):
            name = "block" + str(i)
            new_record = SeqRecord(Seq(blocks[i]), id=name, name=name, description="")
            SeqIO.write(new_record, fa, "fasta")


def get_consensus_seq(cluster_seqs_path, arg_threads):
    from Bio.Align.Applications import ClustalwCommandline
    from Bio.Align.Applications import ClustalOmegaCommandline
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from Bio.Align import MultipleSeqAlignment

    aln_file = '.'.join(cluster_seqs_path.split('.')[:-1]) + "_aln.fasta"
    cmd = ClustalOmegaCommandline(infile=cluster_seqs_path, outfile=aln_file, force=True, threads=arg_threads)
    stdout, stderr = cmd()
    align = AlignIO.read(aln_file, "fasta")

    summary_align = AlignInfo.SummaryInfo(align)
    consensus = summary_align.gap_consensus(threshold=0, ambiguous='N')
    consensus = str(consensus).replace('-', '')
    return consensus

def get_blocks(trpl, path_seq, tsv_res):
    blocks = []
    seqs_dict = {}
    for record in SeqIO.parse(path_seq, "fasta"):
        seqs_dict[record.id] = str(record.seq).upper()

    df_sd = pd.read_csv(tsv_res, "\t")
    print(df_sd.head())
    for i in range(1, len(df_sd) - 1):
        if df_sd.iloc[i,4] > 60:
            if df_sd.iloc[i, 1].rstrip("'") == trpl[1]:
                if df_sd.iloc[i - 1, 1].rstrip("'") == trpl[0] and df_sd.iloc[i + 1, 1].rstrip("'") == trpl[2]:
                    blocks.append(seqs_dict[df_sd.iloc[i,0]][df_sd.iloc[i,2]:(df_sd.iloc[i,3] + 1)])
                    if df_sd.iloc[i, 1][-1] == "'":
                        blocks[-1] = rc(blocks[-1])
    return blocks

def SplitMonomers(MnToSplit, mnpath,  sdtsv, path_seq, outd):
    mnlist = utils.load_fasta(mnpath)
    for mn in MnToSplit.keys():
        resmns = [mon for mon in mnlist if mon.id != mn]
        ci = 0
        for ctx in MnToSplit[mn]:
            blocks = get_blocks((ctx[0], mn, ctx[1]), path_seq, sdtsv)
            save_seqs(blocks, os.path.join(outd, "blseq.fa"))
            consensus = get_consensus_seq(os.path.join(outd, "blseq.fa"), 16)
            name = mn + "." + str(ci)
            ci += 1
            new_record = SeqRecord(Seq(consensus), id=name, name=name, description="")
            resmns.append(new_record)

        mnlist = resmns
    utils.savemn(os.path.join(outd, "mn.fa"), mnlist)


def ElCycleSplit(mn_path, seq_path, sd_tsv, outd, G, hybridSet, threads):
    cycles = DetectHOR.genAllCycles(G)
    mns_prm = [v for v in G.nodes() if v not in hybridSet]
    cl_all = []
    for cl in cycles:
        usedV = set([v for v in cl if v not in hybridSet])
        if usedV == set(mns_prm):
            cl_all = cl
            break

    if len(cl_all) == 0:
        return None

    cl_all = cl_all[:-1]
    if len(cl_all) == len(set(cl_all)):
        return None

    outElC = os.path.join(outd, "ElCycleSplit")
    if not os.path.exists(outElC):
        os.makedirs(outElC)

    MnSplit = {mn: [] for mn in cl_all if cl_all.count(mn) > 1}
    for i, v in enumerate(cl_all):
        if cl_all.count(v) > 1:
            MnSplit[v].append((cl_all[i - 1], cl_all[(i + 1)%len(cl_all)]))

    SplitMonomers(MnSplit, mn_path, sd_tsv, seq_path, outElC)
    tsv_res = run_SD(os.path.join(outElC, "mn.fa"), seq_path, outElC, threads)
    return outElC