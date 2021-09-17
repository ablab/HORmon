#!/usr/bin/env python3

import os
import argparse
from shutil import copy

import HORmon_pipeline.utils as utils
import HORmon_pipeline.ExtractValuableMonomers as ValMon
import HORmon_pipeline.MergeAndSplitMonomers as MergeSplit
import HORmon_pipeline.DrawMonomerGraph as dmg
import HORmon_pipeline.DetectHOR as DetectHOR
import HORmon_pipeline.Hybrid as hybrid
import HORmon_pipeline.ElCycleDecomposition as elCycl
import HORmon_pipeline.RenameMonomers as rename
import HORmon_pipeline.MonoRun as monorun
import HORmon_pipeline.BuildSimpleGraph as smpGr
import HORmon_pipeline.TriplesMatrix as tm

def parse_args():
    parser = argparse.ArgumentParser(description="HORmon: updating monomers to make it consistent with CE postulate, and canonical HOR inferencing")
    parser.add_argument("--mon", dest="mon", help="path to initial monomers", required=True)
    parser.add_argument("--seq", dest="seq", help="path to centromere sequence", required=True)
    parser.add_argument("--cen-id", dest="cenid", help="chromosome id", default="1")
    parser.add_argument("--monomer-thr", dest="vertThr", help="Minimum weight for valuable monomers", type=int, default=100)
    parser.add_argument("--min-edge-multiplicity",
                        "--edge-thr",
                        dest="edgeThr",
                        help="MinEdgeMultiplicite. "
                             "Remove edges in monomer-graph with multiplicite below min(MinEdgeMultiplicite, "
                             "minCountFraction * (min occurrence of valuable monomer))", type=int, default=100)
    parser.add_argument("--min-count-fraction",
                        "--edgeFr-thr",
                        dest="edgeFrThr",
                        help="minCountFraction. "
                             "Remove edges in monomer-graph with multiplicite below min(MinEdgeMultiplicite, "
                             "minCountFraction * (min occurrence of valuable monomer))", type=float, default=0.9)
    parser.add_argument("--min-traversals", dest="minTraversals",
                        help="minimum HOR(or monocycle) occurance",
                        type=int, default=10)
    parser.add_argument("--original_mn", dest = "IAmn", help="path to original monomer only for comparing", default="")
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


def getMonomerGraphEdgeThr(sdout, args):
    k_cnt = tm.calc_mn_order_stat(sdout, maxk=2)
    vcnt = k_cnt[0]
    minW = min(args.edgeThr, min([cnt for v, cnt in vcnt.items()]) * args.edgeFrThr)
    return minW

def main():
    args = parse_args()
    args.outdir = os.path.realpath(args.outdir)
    args.mon = os.path.realpath(args.mon)
    args.seq = os.path.realpath(args.seq)

    if args.IAmn != "":
        args.IAmn = os.path.realpath(args.IAmn)

    valMon, valMonPath = getValuableMonomers(args)
    valMonDir = os.path.dirname(valMonPath)

    if (not os.path.exists(os.path.join(args.outdir, "init"))):
        os.makedirs(os.path.join(args.outdir, "init"))

    sdout = os.path.join(valMonDir, "final_decomposition.tsv")
    dmg.BuildAndDrawMonomerGraph(args.mon, sdout,
                                 os.path.join(args.outdir, "init"),
                                 nodeThr=args.vertThr,
                                 edgeThr=getMonomerGraphEdgeThr(sdout, args))

    mergeSplDir = os.path.join(args.outdir, "merge_split")
    if (not os.path.exists(mergeSplDir)):
        os.makedirs(mergeSplDir)

    mon, mon_path = MergeSplit.MergeSplitMonomers(valMonPath, args.seq, mergeSplDir, args.threads, args.cenid)

    MonDir = os.path.dirname(mon_path)
    sdout = os.path.join(MonDir, "i0", "InitSD", "final_decomposition.tsv")
    dmg.BuildAndDrawMonomerGraph(valMon, sdout, valMonDir, nodeThr=args.vertThr,
                                 edgeThr=getMonomerGraphEdgeThr(sdout, args))

    sdout =  os.path.join(MonDir, "fdec.tsv")
    G = dmg.BuildAndDrawMonomerGraph(mon_path, sdout, MonDir,
                                     nodeThr=args.vertThr,
                                     edgeThr=getMonomerGraphEdgeThr(sdout, args), IAmn=args.IAmn)

    SG = smpGr.BuildSimpleGraph({}, args.seq, sdout, mon_path, edgeThr=getMonomerGraphEdgeThr(sdout, args), IAmn=args.IAmn)
    dmg.DrawMonomerGraph(SG, MonDir, "simpl_graph")

    mnrundir = os.path.join(args.outdir, "MonoRunRaw")
    if not os.path.exists(mnrundir):
        os.makedirs(mnrundir)

    fdec=os.path.join(MonDir, "fdec.tsv")
    monorun.BuildAndShowMonorunGraph(fdec, mnrundir, vLim=args.vertThr,  eLim=getMonomerGraphEdgeThr(fdec, args))

    fdec = os.path.join(MonDir, "fdec.tsv")
    hybridSet, hybridDict = hybrid.getHybridINFO(mon_path, fdec, getMonomerGraphEdgeThr(fdec, args))
    print("Hybrid: ", hybridSet)
    eDir = elCycl.ElCycleSplit(mon_path, args.seq, fdec, args.outdir, G, hybridSet, args.threads)
    if eDir is not None:
        mon_path = os.path.join(eDir, "mn.fa")
        fdec = os.path.join(eDir, "final_decomposition.tsv")

        mns = utils.load_fasta(mon_path)
        mns, mon_path, fdec = utils.updateMonomersByMonomerBlocks(mns, os.path.join(args.outdir, "finalMnUpdate"), fdec, args.seq, int(args.threads))

        G = dmg.BuildAndDrawMonomerGraph(mon_path, fdec, eDir,
                                         nodeThr=args.vertThr,
                                         edgeThr=getMonomerGraphEdgeThr(fdec, args), IAmn=args.IAmn)
        hybridSet, hybridDict = hybrid.getHybridINFO(mon_path, fdec, getMonomerGraphEdgeThr(fdec, args))
        SG = smpGr.BuildSimpleGraph({}, args.seq, fdec, mon_path, edgeThr=getMonomerGraphEdgeThr(fdec, args), IAmn=args.IAmn)
        dmg.DrawMonomerGraph(SG, os.path.join(args.outdir, "finalMnUpdate"), "simpl_graph")

    SG = smpGr.BuildSimpleGraph(hybridSet, args.seq, fdec, mon_path, edgeThr=getMonomerGraphEdgeThr(fdec, args))
    HORs = DetectHOR.detectHORs(mon_path, fdec, args.outdir, SG, hybridSet, args.minTraversals)
    if args.IAmn != "":
        finalOriginalDir = os.path.join(args.outdir, "finalOriginal")
        if not os.path.exists(finalOriginalDir):
            os.makedirs(finalOriginalDir)
        copy(mon_path, finalOriginalDir)
        fdec = utils.run_SD(mon_path, args.seq, finalOriginalDir, args.threads)
        G = dmg.BuildAndDrawMonomerGraph(mon_path, fdec, finalOriginalDir,
                                     nodeThr=args.vertThr,
                                     edgeThr=getMonomerGraphEdgeThr(fdec, args), IAmn=args.IAmn)
        DetectHOR.saveHOR(HORs, finalOriginalDir)

    newNames = rename.RenameMonomers(HORs, hybridDict)
    mon_path = rename.saveNewMn(mon_path, newNames, args.outdir)
    fdec = utils.run_SD(mon_path, args.seq, args.outdir, args.threads)
    G = dmg.BuildAndDrawMonomerGraph(mon_path, fdec, args.outdir,
                                     nodeThr=args.vertThr,
                                     edgeThr=getMonomerGraphEdgeThr(fdec, args))

    HORs = rename.updateHORs(HORs, newNames)
    DetectHOR.saveHOR(HORs, args.outdir)

    mnrundir = os.path.join(args.outdir, "MonoRun")
    if not os.path.exists(mnrundir):
        os.makedirs(mnrundir)
    monorun.BuildAndShowMonorunGraph(fdec, mnrundir, vLim=args.vertThr, eLim=getMonomerGraphEdgeThr(fdec, args))


if __name__ == "__main__":
    main()
