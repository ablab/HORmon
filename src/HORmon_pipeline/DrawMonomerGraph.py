#!/usr/bin/env python3
import math

import networkx as nx
from networkx.drawing.nx_agraph import write_dot
from subprocess import check_call
import numpy as np
import os
import HORmon_pipeline.utils as utils
import HORmon_pipeline.TriplesMatrix as tm

def DrawMonomerGraph(G, outdir):
    write_dot(G, os.path.join(outdir, "graph.dot"))
    try:
        check_call(['dot', '-Tpng', os.path.join(outdir, "graph.dot"), '-o', os.path.join(outdir, "graph.png")])
    except Exception:
        return


def getEdgeThr(k2cnt):
    vcnt = {v : 0 for v, u in k2cnt.keys()}
    for vu, cnt in k2cnt.items():
        vcnt[vu[0]] += cnt
    minW = min(100, min([cnt for v, cnt in vcnt.items()])*0.9)
    return minW


def buildk_graph(kcnt1, kcnt2, k=1, nodeThr=100, edgeThr=100):
    mn_set = {tuple(list(x)[:-1]) for x, y in kcnt2.items() if y > nodeThr} | \
             {tuple(list(x)[1:]) for x, y in kcnt2.items() if y > nodeThr}

    if k == 1:
        cntmn2 = {(mn,): cnt for mn, cnt in kcnt1.items()}
        kcnt1 = cntmn2

    G = nx.DiGraph()
    for vert in mn_set:
        cnt_mn = 0
        if vert in kcnt1:
            cnt_mn = kcnt1[vert]
            lg = 0
            if cnt_mn > 0:
                lg = math.log(cnt_mn)
        clr = ["red", "#cd5700", "orange", "#ffb300", "yellow", "#ceff1d", "#a7fc00", "#00ff00", "#e0ffff", "#f5fffa",
               "#f5fffa", "#f5fffa", "#f5fffa", "#f5fffa"]
        vert = "-".join(list(vert))
        curc = clr[int(lg)]
        G.add_node(vert, style="filled", fillcolor=curc, label=f'"{vert}[{str(int(cnt_mn))}]"')

    ecnt = 0
    for vt1 in sorted(mn_set):
        for vt2 in sorted(mn_set):
            if list(vt1)[1:] != list(vt2)[:-1]:
                continue
            scr = 0
            if (*vt1, vt2[-1]) in kcnt2:
                scr = kcnt2[(*vt1, vt2[-1])]

            thr_wg = [200, 100, 50, 25, 1] #[100000000, 1000, 500, 100, 1]
            wgs = [7, 5, 3, 1, 0]
            wg = 3
            while (scr > thr_wg[wg]):
                wg -= 1

            if scr > nodeThr:
                ecnt += 1
                vrt1 = "-".join(list(vt1))
                vrt2 = "-".join(list(vt2))
                if scr < 1:
                    G.add_edge(vrt1, vrt2, label=str(int(scr)), penwidth=str(wgs[wg]), constraint="false")
                else:
                    G.add_edge(vrt1, vrt2, label=str(int(scr)), penwidth=str(wgs[wg]))

    return G


def BuildMonomerGraph(path_to_mn, sdout, nodeThr=100, edgeThr="auto"):
    k_cnt = tm.calc_mn_order_stat(sdout, maxk=2)
    if edgeThr=="auto":
        edgeThr = getEdgeThr(k_cnt[1])
    else:
        edgeThr = int(edgeThr)

    G = buildk_graph(k_cnt[0], k_cnt[1], 1, nodeThr, edgeThr)
    return G

def BuildAndDrawMonomerGraph(path_to_mn, sdout, outdir, nodeThr=100, edgeThr="auto"):
    G = BuildMonomerGraph(path_to_mn, sdout, nodeThr=nodeThr, edgeThr=edgeThr)
    DrawMonomerGraph(G, outdir)