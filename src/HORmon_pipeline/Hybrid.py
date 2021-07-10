#!/usr/bin/env python3

import HORmon_pipeline.utils as utils
import HORmon_pipeline.TriplesMatrix as TriplesMatrix

def get_hybrid_len(main_mn, mn1, mn2):
    resDiv = 500
    mn_identity = 100
    bst_res = (0, 0)
    for prfx in range(30, len(mn1.seq)):
        suffix = len(str(main_mn.seq)) - prfx
        if suffix < 30:
            break
        hbr = str(mn1.seq)[:prfx] + str(mn2.seq)[len(mn2.seq) - suffix:]
        cur_identity = utils.seq_identity(hbr, str(main_mn.seq))
        if mn_identity > cur_identity:
            mn_identity = cur_identity
            bst_res =  (prfx, suffix)

    if (mn_identity * 2 <= resDiv):
        return (bst_res[0], bst_res[1], mn_identity)
    return (0, 0, 100)


def isHybrid(main_mn, mn1, mn2, main_rd):
    pr, sf, idn = get_hybrid_len(main_mn, mn1, mn2)
    if idn < 6 and idn + 2 < main_rd:
        print("Hybrid", main_mn.id, main_rd, idn, mn1.id, mn2.id, pr, sf)
        return True
    return False

def getHybridINFO(mnpath, decpath):
    vcnt = TriplesMatrix.calc_mn_order_stat(decpath, maxk=2)[0]
    hybridSet = set()
    monCA = utils.load_fasta(mnpath)

    rd = {mn.id: 100 for mn in monCA}
    for i in range(len(monCA)):
        for j in range(len(monCA)):
            if i != j:
                rd[monCA[i].id] = min(rd[monCA[i].id], utils.seq_identity(monCA[i].seq, monCA[j].seq))

    for i in range(len(monCA)):
        for j in range(len(monCA)):
            for g in range(len(monCA)):
                if i != j and j != g and i != g:
                    if vcnt[monCA[j].id] > vcnt[monCA[i].id] and vcnt[monCA[g].id] > vcnt[monCA[i].id]:
                        if isHybrid(monCA[i], monCA[j], monCA[g], rd[monCA[i].id]):
                            hybridSet.add(monCA[i].id)
    return hybridSet