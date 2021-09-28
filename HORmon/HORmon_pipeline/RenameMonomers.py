#!/usr/bin/env python3

import os
import HORmon.HORmon_pipeline.utils as utils

def getNameById(id):
    if id == 0:
        return 'A'

    fn = ""
    while id > 0:
        fn = str(chr(ord('A') + int(id % 26))) + fn
        id //= 26

    return fn


def RenameMonomers(HORs, HybridDict):
    newNames = {}

    lstid = 0
    for i in range(len(HORs)):
        for j in range(len(HORs[i])):
            if HORs[i][j] not in newNames:
                newNames[HORs[i][j]] = getNameById(lstid)
                lstid += 1

    print("HORs", HORs)
    print("New names:", newNames)
    print("Hybrid dict:", HybridDict)
    for mn, vl in HybridDict.items():
        nm = getNameById(lstid)
        lstid += 1
        nm += f'({newNames[vl[0]]}/{newNames[vl[1]]})'
        newNames[mn] = nm
    return newNames


def saveNewMn(mn_path, newNames, outdir):
    mns = utils.load_fasta(mn_path)
    for i in range(len(mns)):
        if mns[i].id in newNames:
            mns[i].id = newNames[mns[i].id]

    utils.savemn(os.path.join(outdir, "mn.fa"), mns)
    return os.path.join(outdir, "mn.fa")

def updateHORs(HORs, newNames):
    for i in range(len(HORs)):
        for j in range(len(HORs[i])):
            HORs[i][j] = newNames[HORs[i][j]]
    return HORs