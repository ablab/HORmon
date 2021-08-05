#!/usr/bin/env python3
import math

import os
import networkx as nx
from networkx.drawing.nx_agraph import write_dot
from subprocess import check_call
import numpy as np

import HORmon_pipeline.utils as utils
import HORmon_pipeline.TriplesMatrix as TriplesMatrix

class LongEdge:
    def __init__(self):
        self.epath = []
        self.weight = 0
        self.name = ""

def initLongEdgeWeight(outLongE, k2cnt):
    return


def initLongEdgeNames(outLongE):
    namecnt = {}
    for v, les in outLongE.items():
        for le in les:
            cname = "L" + str(len(le.epath) - 1)
            if le.epath[0] == le.epath[-1]:
                cname = "C" + str(len(le.epath) - 1)
            if cname not in namecnt:
                namecnt[cname] = 1
                le.name = cname
            else:
                le.name = cname + "-" + str(namecnt[cname])
                namecnt[cname] += 1


def canSplit(v, srunG, rG, epaths, k3cnt, handledV):
    oute = srunG.edges(v)
    ine = rG.edges(v)
    odeg = len(oute)
    ideg = len(ine)
    if ideg == 1 or odeg == 1:
        return True
    return False


def SplitV(v, mnrunG, rG, epaths, k3cnt, handledV, srunG):
    oute = list(srunG.edges(v))
    ine = list(rG.edges(v))
    odeg = len(oute)
    ideg = len(ine)

    alph = "".join([chr(ord('a') + i) for i in range(26)])
    vlist = []
    srunG.remove_node(v)
    rG.remove_node(v)
    for e0 in oute:
        for e1 in ine:
            vlist.append(v + alph[len(vlist)])
            srunG.add_node(vlist[-1])

            srunG.add_edge(e1[1], vlist[-1])
            print(v, e0[1])
            u = e0[1]
            if e0[1][-1] in alph:
                u = e0[1][:-1]

            srunG.add_edge(vlist[-1], e0[1])
            v1 = e1[1]
            if e1[1][-1] in alph:
                v1 = e1[1][:-1]

            u1 = v
            if len(ine) == 1:
                v1 = v
                u1 = u
            print(list(srunG.nodes()))
            print(list(srunG.edges()))
            print(srunG)
            print(e1)
            print(list(mnrunG.nodes()))
            srunG[e1[1]][vlist[-1]]["penwidth"] = mnrunG[v1][u1]["penwidth"]
            srunG[e1[1]][vlist[-1]]["label"] = mnrunG[v1][u1]["label"]
            srunG[vlist[-1]][e0[1]]["penwidth"] = mnrunG[v1][u1]["penwidth"]
            srunG[vlist[-1]][e0[1]]["label"] = mnrunG[v1][u1]["label"]

            rG.add_edge(vlist[-1], e1[1])
            rG.add_edge(e0[1], vlist[-1])
    return vlist



def SplitMnrunVert(mnrunG, epaths, k3cnt):
    srunG = nx.DiGraph()
    rG = nx.DiGraph()
    for v in mnrunG.nodes():
        rG.add_node(v)
        srunG.add_node(v)

    for e in mnrunG.edges():
        rG.add_edge(e[1], e[0])
        srunG.add_edge(e[0], e[1])
        srunG[e[0]][e[1]]["penwidth"] = mnrunG[e[0]][e[1]]["penwidth"]
        srunG[e[0]][e[1]]["label"] = mnrunG[e[0]][e[1]]["label"]

    handledV = {}

    for v in mnrunG.nodes():
        oute = mnrunG.edges(v)
        ine = rG.edges(v)
        odeg = len(srunG.edges(v))
        ideg = len(ine)

        if odeg < 2 and ideg < 2:
            handledV[v] = [v]
        elif canSplit(v, srunG, rG, epaths, k3cnt, handledV):
            handledV[v] = SplitV(v, mnrunG, rG, epaths, k3cnt, handledV, srunG)
        else:
            handledV[v] = [v]

    return srunG


def genCycleInner(mnrunG, prefixCycle, usedEdges, usedV, cycleList):
    def samecc(c1, c2):
        if len(c1) != len(c2):
            return False
        for j in range(len(c1) - 1):
            if c1[j:-1] + c1[:j] == c2[:-1]:
                return True
        return False

    if len(prefixCycle) > 2 and prefixCycle[-1] == prefixCycle[-2]:
        return

    if len(prefixCycle) > 1 and prefixCycle[0] == prefixCycle[-1]:
        for cc in cycleList:
            if samecc(cc, prefixCycle):
                break
        else:
            cycleList.append(prefixCycle)

    if len(prefixCycle) > 1 and prefixCycle[-1] == prefixCycle[-2]:
        return

    v = prefixCycle[-1]
    for e in mnrunG.edges(v):
        if e not in usedEdges and e[1] not in usedV:
            genCycleInner(mnrunG, prefixCycle + [e[1]], usedEdges | {e}, usedV, cycleList)


def genAllCycles(mnrunG):
    usedV = set()
    cycleList = []
    for v in mnrunG.nodes():
        if v not in usedV:
            genCycleInner(mnrunG, [v], set(), usedV, cycleList)
            usedV.add(v)

    return cycleList


def monomrunHOR2monomersHOR(cc, epaths):
    print(cc)
    mncc = []
    for i in range(len(cc) - 1):
        mncc += epaths[cc[i]][:-1]
    return mncc


def getHORcnt(mnHOR, monocen):
    cnt = 0
    for i in range(0, len(monocen) - len(mnHOR)):
        if mnHOR == monocen[i:i + len(mnHOR)]:
            cnt += 1
    return cnt


def detectNodes(mnrunG):
    inE = {v: [] for v in mnrunG.nodes()}
    outE = {v: [] for v in mnrunG.nodes()}

    for e in mnrunG.edges():
        inE[e[1]].append(e[0])
        outE[e[0]].append(e[1])

    usedV = set()
    Nodes = []
    for v in mnrunG.nodes():
        if v in usedV:
            continue

        if len(inE[v]) != 1 or len(outE[v]) != 1:
            usedV |= {v}
            for u in outE[v]:
                curE = [v]
                w = u
                while len(inE[w]) == 1 and len(outE[w]) == 1:
                    usedV |= {w}
                    curE.append(w)
                    w = outE[w][0]

                curE.append(w)
                Nodes.append(curE)

    for v in mnrunG.nodes():
        if v in usedV:
            continue

        usedV |= {v}
        for u in outE[v]:
            curE = [v]
            w = u
            while w != v:
                usedV |= {w}
                curE.append(w)
                w = outE[w]

            curE.append(w)
            Nodes.append(curE)
    return Nodes


def getWeight(NodesList, mnruns2mn, mncen):
    mnlist = []
    for nd in NodesList:
        mnlist += mnruns2mn[nd]
        mnlist = mnlist[:-1]

    cntIner = 0
    for i in range(len(mncen) - len(mnlist)):
        if mncen[i:i+len(mnlist)] == mnlist:
            cntIner += 1
        if [nm[:-1] for nm in mncen[i+len(mnlist)-1:i-1:-1]] == mnlist:
            cntIner += 1
    return cntIner

def BuildAndShorIterativeMonorunGraph(mnrunG, MonorunsToMonomers, Monocentromere, outdir, filenm="IterMnrun"):
    print("Iterative monorun graph:")

    iterMnNodes = detectNodes(mnrunG)
    print("Iterative Monorun Nodes:", iterMnNodes)

    mnIterG = nx.DiGraph()
    for nds in iterMnNodes:
        wgh = str(getWeight(nds, MonorunsToMonomers, Monocentromere))
        mnIterG.add_node(','.join(nds), label=''.join(nds[1:]) + "[" + wgh + "]")

    for nds1 in iterMnNodes:
        for nds2 in iterMnNodes:
            if nds1[-1] == nds2[0]:
                wgh = getWeight(nds1 + nds2[1:], MonorunsToMonomers, Monocentromere)
                if wgh > 0:
                    mnIterG.add_edge(','.join(nds1), ','.join(nds2), label=str(wgh), penwidth=str(math.log10(wgh)))

    ofile = os.path.join(outdir, filenm + ".dot")
    opng = os.path.join(outdir, filenm + ".png")
    write_dot(mnIterG, ofile)
    try:
        check_call(['dot', '-Tpng', ofile, '-o', opng])
    except Exception:
        return

    return iterMnNodes




def BuildAndShowMonorunGraph(tsv_res, outdir, vLim=100, eLim = 100):
    k_cnt = TriplesMatrix.calc_mn_order_stat(tsv_res, maxk=3)
    k2cnt = k_cnt[1]
    k3cnt = k_cnt[2]

    monocen = utils.get_monocent(tsv_res)

    vcnt = {v : 0 for v, u in k2cnt.keys()}
    eCnt = 0
    for vu, cnt in k2cnt.items():
        vcnt[vu[0]] += cnt
        if cnt > eLim:
            eCnt += 1


    #print(vcnt)

    ine = {v : [] for v, cnt in vcnt.items() if cnt >= vLim}
    oute = {v : [] for v, cnt in vcnt.items() if cnt >= vLim}

    for vu, cnt in k2cnt.items():
        if cnt < eLim:
            continue
        ine[vu[1]].append(vu[0])
        oute[vu[0]].append(vu[1])

    usedv = set()

    outLongE = {v: [] for v in oute.keys()}
    for v in ine.keys():
        if v not in usedv and (len(oute[v]) != 1 or len(ine[v]) != 1):
            usedv |= {v}
            for u in oute[v]:
                curE = LongEdge()
                outLongE[v].append(curE)
                curE.epath.append(v)
                cu = u
                while len(oute[cu]) == 1 and len(ine[cu]) == 1:
                    usedv.add(cu)
                    curE.epath.append(cu)
                    cu = oute[cu][0]
                curE.epath.append(cu)

    for v in ine.keys():
        if v not in usedv:
            usedv |= {v}
            curE = LongEdge()
            outLongE[v].append(curE)
            curE.epath.append(v)
            u = oute[v][0]
            while u != v:
                usedv |= {u}
                curE.epath.append(u)
                u = oute[u][0]
            curE.epath.append(v)

    initLongEdgeWeight(outLongE, k2cnt)
    initLongEdgeNames(outLongE)

    lesall = []
    mnrunG = nx.DiGraph()
    for v, les in outLongE.items():
        for le in les:
            lesall.append(le)
            mnrunG.add_node(le.name)

    def smpl(vr):
        if len(vr) > 1:
            vr = f'({vr})'
        return vr
        
    open(os.path.join(outdir, "L.csv"), "w").close()

    epaths = {}
    for le in lesall:
        with open(os.path.join(outdir, "L.csv"), "a") as fw:
            multipl = []
            for i in range(1, len(le.epath)):
                multipl.append(k2cnt[(le.epath[i - 1], le.epath[i])])
            multipl.sort()
            fw.write(le.name + "\t" + "".join([smpl(vr) for vr in le.epath]) + "\t" + str(min(multipl)) + "\n")
        epaths[le.name] = le.epath
        for le2 in lesall:
            if le.epath[-1] == le2.epath[0]:
                if (le.epath[-2], le2.epath[0], le2.epath[1]) in k3cnt and \
                        k3cnt[(le.epath[-2], le2.epath[0], le2.epath[1])] >= eLim:
                    mnrunG.add_edge(le.name, le2.name,
                                    penwidth=max(1, 1.5*(math.log(float(k3cnt[(le.epath[-2], le2.epath[0], le2.epath[1])])) - 4)),
                                    label=str(int(k3cnt[(le.epath[-2], le2.epath[0], le2.epath[1])])))

    ofile = os.path.join(outdir, "monorun.dot")
    write_dot(mnrunG, ofile)
    try:
        check_call(['dot', '-Tpng', ofile, '-o', ".".join(ofile.split(".")[:-1]) + ".png"])
    except Exception:
        return

    with open(os.path.join(outdir, "Monorunscnt.csv"), "w") as fw:
        fw.write("Edges #:\t" + str(eCnt) + "\n")
        fw.write("Monoruns #:\t" + str(len(mnrunG.nodes())) + "\n")

    srunG = SplitMnrunVert(mnrunG, epaths, k3cnt)
    sruno = ".".join(ofile.split(".")[:-1]) + "_splv.dot"
    srunopng = ".".join(ofile.split(".")[:-1]) + "_splv.png"
    write_dot(srunG, sruno)
    try:
        check_call(['dot', '-Tpng', sruno, '-o', srunopng])
    except Exception:
        return

    MonorunsToMonomers = {le.name: le.epath for le in lesall}
    for nd in srunG.nodes():
        if nd[-1].isalpha():
            MonorunsToMonomers[nd] = MonorunsToMonomers[nd[:-1]]

    BuildAndShorIterativeMonorunGraph(mnrunG, MonorunsToMonomers, monocen, outdir)
    BuildAndShorIterativeMonorunGraph(srunG, MonorunsToMonomers, monocen, outdir, "SplitIterMnrun")