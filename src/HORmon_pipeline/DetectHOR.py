#!/usr/bin/env python3

import os
import csv
import HORmon_pipeline.utils as utils

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


def filterCycles(cycles, hybridSet):
    usedV = set()
    cycles.sort(key=lambda x: -len(x))
    res_cyc = []
    for cl in cycles:
        if len([v for v in cl if v in hybridSet]) > 0:
            continue

        for v in cl:
            if v not in usedV:
                res_cyc.append(cl)
                break
        usedV |= set(cl)
    return res_cyc


def detectHORs(mon_path, sd_path, outdir, G, hybridSet):
    mncen = utils.get_monocent(sd_path)
    cycles = genAllCycles(G)
    cycles = filterCycles(cycles, hybridSet)
    with open(os.path.join(outdir, "HORs.tsv"), "w") as fw:
        csv_writer = csv.writer(fw, delimiter='\t')
        horid = 1
        for cycle in cycles:
            csv_writer.writerow(["H" + str(horid), ",".join(cycle)])
            horid += 1