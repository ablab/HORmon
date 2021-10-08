#!/usr/bin/env python3

import os
import csv

import HORmon.HORmon_pipeline.utils as utils


def genCycleInner(G, prefixCycle, usedEdges, cycleList):
    def samecc(c1, c2):
        if len(c1) != len(c2):
            return False
        for j in range(len(c1) - 1):
            if c1[j:-1] + c1[:j] == c2[:-1]:
                return True
        return False

    if len(prefixCycle) > 1 and prefixCycle[0] == prefixCycle[-1]:
        for cc in cycleList:
            if samecc(cc, prefixCycle):
                break
        else:
            cycleList.append(prefixCycle)

    v = prefixCycle[-1]
    for e in G.edges(v):
        if e not in usedEdges:
            genCycleInner(G, prefixCycle + [e[1]], usedEdges | {e}, cycleList)


def genAllCycles(G):
    cycleList = []
    for v in G.nodes():
        genCycleInner(G, [v], set(), cycleList)
    return cycleList


def getCycleCnt(cl, mncen):
    clcnt = 0
    for i in range(len(mncen) - len(cl)):
        if mncen[i:i+len(cl)] == cl:
            clcnt += 1

    ccl = [x + "'" for x in cl[len(cl)::-1]]

    for i in range(len(mncen) - len(cl)):
        if mncen[i:i+len(cl)] == ccl:
            clcnt += 1

    return clcnt

def filterCycles(cycles, hybridSet, mncen, minTrav):
    usedV = set()
    cycles.sort(key=lambda x: -len(x))
    res_cyc = []
    for cl in cycles:
        if len([v for v in cl if v in hybridSet]) > 0:
            continue

        if getCycleCnt(cl, mncen) < minTrav:
            continue

        for v in cl:
            if v not in usedV:
                res_cyc.append(cl)
                break
        usedV |= set(cl)
    return res_cyc


def detectHORs(mon_path, sd_path, outdir, G, hybridSet, minTrav):
    mncen = utils.get_monocent(sd_path)
    with open(os.path.join(outdir, "Monocen"), "w") as fw:
        fw.write(str(mncen))

    cycles = genAllCycles(G)
    cycles = filterCycles(cycles, hybridSet, mncen, minTrav)
    return cycles

def saveHOR(cycles, outdir):
    outfile = os.path.join(outdir, "HORs.tsv")
    with open(outfile, "w") as fw:
        csv_writer = csv.writer(fw, delimiter='\t')
        horid = 1
        for cycle in cycles:
            csv_writer.writerow(["H" + str(horid), ",".join(cycle)])
            horid += 1
    return outfile