#!/usr/bin/env python3

import os
import sys

import HORmon_pipeline.DrawMonomerGraph as dmg

mon_path, fdec, outdir = sys.argv[1], sys.argv[2], sys.argv[3]

if (not os.path.exists(outdir)):
	os.makedirs(outdir)
G = dmg.BuildAndDrawMonomerGraph(mon_path, fdec, outdir)