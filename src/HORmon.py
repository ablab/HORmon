#!/usr/bin/env python3

import os
import argparse

import HORmon_pipeline.utils as utils
import HORmon_pipeline.ExtractValuableMonomers as ValMon
import HORmon_pipeline.MergeAndSplitMonomers as MergeSplit
import HORmon_pipeline.DrawMonomerGraph as dmg
import HORmon_pipeline.DetectHOR as DetectHOR

def parse_args():
    parser = argparse.ArgumentParser(description="HORmon: updating monomers to make it consistent with CE postulate, and canonical HOR inferencing")
    parser.add_argument("--mon", dest="mon", help="path to initial monomers", required=True)
    parser.add_argument("--seq", dest="seq", help="path to centromere sequence", required=True)
    parser.add_argument("-t", dest="threads", help="number of threads(default=1)", default=1, type=int)
    parser.add_argument("-o", dest = "outdir", help="path to output directore", required=True)

    return parser.parse_args()


def getValuableMonomers(args):
    ValMonDir = os.path.join(args.outdir, "valuable_monomers")
    if (not os.path.exists(ValMonDir)):
        os.makedirs(ValMonDir)

    valDec = utils.run_SD(args.mon, args.seq, ValMonDir, args.threads)
    initMon = utils.load_fasta(args.mon)
    return ValMon.getValuableMonomers(args.seq, valDec, initMon, ValMonDir)

def main():
    args = parse_args()
    args.outdir = os.path.realpath(args.outdir)
    args.mon = os.path.realpath(args.mon)
    args.seq = os.path.realpath(args.seq)

    valMon, valMonPath = getValuableMonomers(args)
    valMonDir = os.path.dirname(valMonPath)
    dmg.BuildAndDrawMonomerGraph(valMon, os.path.join(valMonDir, "final_decomposition.tsv"), valMonDir, nodeThr=0, edgeThr=10)

    mergeSplDir = os.path.join(args.outdir, "merge_split")
    if (not os.path.exists(mergeSplDir)):
        os.makedirs(mergeSplDir)

    mon, mon_path = MergeSplit.MergeSplitMonomers(valMonPath, args.seq, mergeSplDir, args.threads)
    MonDir = os.path.dirname(mon_path)
    G = dmg.BuildAndDrawMonomerGraph(mon_path, os.path.join(MonDir, "fdec.tsv"), MonDir, nodeThr=0, edgeThr=10)

    DetectHOR.detectHORs(mon_path, os.path.join(MonDir, "fdec.tsv"), args.outdir, G)


if __name__ == "__main__":
    main()