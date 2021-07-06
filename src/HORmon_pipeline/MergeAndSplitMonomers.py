#!/usr/bin/env python3

import os
import shutil
import pandas as pd
import csv

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import HORmon_pipeline.utils as utils
from HORmon_pipeline.utils import rc
from HORmon_pipeline.utils import unique
from HORmon_pipeline.utils import load_fasta
from HORmon_pipeline.utils import run_SD
import HORmon_pipeline.TriplesMatrix as TriplesMatrix


def getMnSim(mon):
    sim = {}
    for mn1 in mon:
        for mn2 in mon:
            sim[(mn1.id, mn2.id)] = min(utils.seq_identity(mn1.seq, mn2.seq),
                                        utils.seq_identity(mn1.seq, rc(mn2.seq)))
    return sim


cntMerge = 0
cntSplit = 0
tsv_res = ""

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


def MergeMonomers(mn1, mn2, odir, mons, path_seq):
    print("====MERGE===" + mn1.id + "+" + mn2.id)
    global tsv_res
    resmns = [mn for mn in mons if mn.id != mn1.id and mn.id != mn2.id]

    blocks = utils.get_blocks(path_seq, tsv_res, [mn1.id, mn2.id])
    save_seqs(blocks, os.path.join(odir, "blseq.fa"))
    consensus = get_consensus_seq(os.path.join(odir, "blseq.fa"), 16)
    name = mn1.id + "+" + mn2.id
    new_record = SeqRecord(Seq(consensus), id=name, name=name, description="")
    resmns.append(new_record)
    return resmns


def get_blocks(trpl, path_seq, tsv_res):
    block1, block2 = [],[]
    seqs_dict = {}
    for record in SeqIO.parse(path_seq, "fasta"):
        seqs_dict[record.id] = str(record.seq).upper()

    df_sd = pd.read_csv(tsv_res, "\t")
    for i in range(1, len(df_sd) - 1):
        if df_sd.iloc[i,4] > 60:
            if df_sd.iloc[i, 1].rstrip("'") == trpl[1]:
                #print(df_sd.iloc[i, 1])
                #print(df_sd.iloc[i - 1, 1], df_sd.iloc[i, 1], df_sd.iloc[i + 1, 1])
                #print(trpl)
                if df_sd.iloc[i - 1, 1].rstrip("'") == trpl[0] and df_sd.iloc[i + 1, 1].rstrip("'") == trpl[2]:
                    block2.append(seqs_dict[df_sd.iloc[i,0]][df_sd.iloc[i,2]:(df_sd.iloc[i,3] + 1)])
                    if df_sd.iloc[i, 1][-1] == "'":
                        block2[-1] = rc(block2[-1])
                else:
                    block1.append(seqs_dict[df_sd.iloc[i, 0]][df_sd.iloc[i, 2]:(df_sd.iloc[i, 3] + 1)])
                    if df_sd.iloc[i, 1][-1] == "'":
                        block1[-1] = rc(block1[-1])
    return block1, block2


def SplitMn(trpl, nnm, odir, mons, path_seq):
    print("====SPLIT===", trpl, nnm)
    global tsv_res
    resmns = [mn for mn in mons if mn.id != trpl[1]]

    block1, block2 = get_blocks(trpl, path_seq, tsv_res)
    save_seqs(block1, os.path.join(odir, "blseq.fa"))
    consensus = get_consensus_seq(os.path.join(odir, "blseq.fa"), 16)
    name = trpl[1]
    new_record = SeqRecord(Seq(consensus), id=name, name=name, description="")
    resmns.append(new_record)

    save_seqs(block2, os.path.join(odir, "blseq.fa"))
    consensus = get_consensus_seq(os.path.join(odir, "blseq.fa"), 16)
    name = nnm
    new_record = SeqRecord(Seq(consensus), id=name, name=name, description="")
    resmns.append(new_record)

    return resmns


def Iteration(iterNum, seq_path, outdir, monsPath, thread=1):
    global tsv_res
    global cntMerge
    global cntSplit

    mons = unique(load_fasta(monsPath))
    sim = getMnSim(mons)

    odir = os.path.join(outdir, "i" + str(iterNum))
    if not os.path.exists(odir):
        os.makedirs(odir)

    tsv_res = run_SD(monsPath, seq_path, os.path.join(odir, "InitSD"), thread)
    k_cnt = TriplesMatrix.calc_mn_order_stat(tsv_res, maxk=3)
    posscore, cenvec = TriplesMatrix.handleAllMn(k_cnt[2], k_cnt[1], thr=0)

    print("Similarity: ", sim)

    def get_best(posscore):
        bstm = (-1, -1)
        bp = 0
        for i in range(len(mons)):
            for j in range(len(mons)):
                mn1 = mons[i]
                mn2 = mons[j]
                if i == j:
                    continue
                if sim[(mn1.id, mn2.id)] < 6:
                    scr = posscore.get((mn1.id, mn2.id), 0)
                    if scr > bp:
                        bp = scr
                        bstm = (i, j)
        return bstm, bp

    bstm, bp = get_best(posscore)

    if bp > 0.3:
        cntMerge += 1
        resmn = MergeMonomers(mons[bstm[0]], mons[bstm[1]], odir, mons, seq_path)
        utils.savemn(os.path.join(odir, "mn.fa"), resmn)
        return True

    posscore, cenvec = TriplesMatrix.PrefixPosScore(k_cnt[2], k_cnt[1], thr=0)
    bstm, bp = get_best(posscore)
    if bp > 0.3:
        cntMerge += 1
        resmn = MergeMonomers(mons[bstm[0]], mons[bstm[1]], odir, mons, seq_path)
        utils.savemn(os.path.join(odir, "mn.fa"), resmn)
        return True

    posscore, cenvec = TriplesMatrix.SuffixPosScore(k_cnt[2], k_cnt[1], thr=0)
    bstm, bp = get_best(posscore)
    if bp > 0.3:
        cntMerge += 1
        resmn = MergeMonomers(mons[bstm[0]], mons[bstm[1]], odir, mons, seq_path)
        utils.savemn(os.path.join(odir, "mn.fa"), resmn)
        return True

    exchTrp = TriplesMatrix.SplitAllMn(k_cnt[2], k_cnt[1])
    if len(exchTrp) > 0:
        for key in exchTrp.keys():
            cntSplit += 1
            resmn = SplitMn(key, exchTrp[key], odir, mons, seq_path)
            utils.savemn(os.path.join(odir, "mn.fa"), resmn)
            return True
    return False


def MergeSplitMonomers(mns_path, seq_path, outdir, thr):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    initmn = unique(load_fasta(mns_path))

    shutil.copyfile(mns_path, os.path.join(outdir, "mn.fa"))
    iterNum = 0
    while (Iteration(iterNum, seq_path, outdir, mns_path, thr)):
        iterNum += 1
        mns_path = os.path.join(outdir, "i" + str(iterNum - 1), "mn.fa")
        shutil.copyfile(mns_path, os.path.join(outdir, "mn.fa"))

    print(os.path.join(outdir, "i" + str(iterNum), "InitSD", "final_decomposition.tsv"))
    shutil.copyfile(os.path.join(outdir, "i" + str(iterNum), "InitSD", "final_decomposition.tsv"),
                    os.path.join(outdir, "fdec.tsv"))

    finalm = unique(load_fasta(mns_path))

    with open(os.path.join(outdir, "MergeSplitStat.csv"), "w") as f:
        csv_writer = csv.writer(f, delimiter=',')
        csv_writer.writerow(["#Init mon", "#Merge", "#Split", "#Final mon", "Final SqMean", "Final DBIndex"])
        csv_writer.writerow([str(len(initmn)), str(cntMerge), str(cntSplit), str(len(finalm)),
                             str(utils.blocks_sqmean(seq_path, tsv_res , finalm)),
                             str(utils.DaviesBouldinIndex(seq_path, tsv_res , finalm))])

    return finalm, mns_path