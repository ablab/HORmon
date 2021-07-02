#!/usr/bin/env python3

import os
import argparse

import HORmon_pipeline.utils as utils
import HORmon_pipeline.ExtractValuableMonomers as ValMon

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


if __name__ == "__main__":
    main()