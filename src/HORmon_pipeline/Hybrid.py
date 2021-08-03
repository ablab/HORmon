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


def isHybridContext(main_mn, mn1, mn2, kcnt2, edgeThr=100):
    mnlst = {x[0] for x in kcnt2.keys()}
    for mn in mnlst:
        if (mn, mn1) in kcnt2 and kcnt2[(mn, mn1)] > edgeThr and \
                (mn, main_mn) in kcnt2 and kcnt2[(mn, main_mn)] > edgeThr:
           break
    else:
        return False

    for mn in mnlst:
        if (mn2, mn) in kcnt2 and kcnt2[(mn2, mn)] > edgeThr and \
                (main_mn, mn) in kcnt2 and kcnt2[(main_mn, mn)] > edgeThr:
           break
    else:
        return False
    return True

def isHybrid(main_mn, mn1, mn2, main_rd, kcnt2):
    pr, sf, idn = get_hybrid_len(main_mn, mn1, mn2)
    if idn < 6 and idn * 1.5 < main_rd and isHybridContext(main_mn.id, mn1.id, mn2.id, kcnt2):
        print("Hybrid", main_mn.id, main_rd, idn, mn1.id, mn2.id, pr, sf)
        return True
    return False

def getHybridINFO(mnpath, decpath):
    kcnt = TriplesMatrix.calc_mn_order_stat(decpath, maxk=2)
    vcnt = kcnt[0]
    kcnt2 = kcnt[1]
    hybridSet = set()
    hybridIdn = {}
    hybridDict = {}
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
                    if isHybrid(monCA[i], monCA[j], monCA[g], rd[monCA[i].id], kcnt2):
                        pr, sf, idn = get_hybrid_len(monCA[i], monCA[j], monCA[g])
                        if monCA[i].id in hybridSet:
                            continue
                        hybridSet.add(monCA[i].id)

    for i in range(len(monCA)):
        for j in range(len(monCA)):
            for g in range(len(monCA)):
                if not(monCA[i].id in hybridSet and (monCA[j].id not in hybridSet) and (monCA[g].id not in hybridSet)):
                    continue

                if i != j and j != g and i != g:
                    if isHybrid(monCA[i], monCA[j], monCA[g], rd[monCA[i].id], kcnt2):
                        pr, sf, idn = get_hybrid_len(monCA[i], monCA[j], monCA[g])
                        if monCA[i].id in hybridIdn and idn > hybridIdn[monCA[i].id]:
                            continue

                        hybridIdn[monCA[i].id] = idn
                        hybridDict[monCA[i].id] = (monCA[j].id, monCA[g].id)

    return hybridSet, hybridDict