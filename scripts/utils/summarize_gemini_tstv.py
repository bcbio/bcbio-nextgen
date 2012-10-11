#!/usr/bin/env python
"""Provide table summarizing Transition/Transversion ratios for variants.
"""
import os
import sys
import glob

import sh
import pandas

def main(work_dir):
    stats = ["tstv", "tstv-coding", "tstv-noncoding", "vars-by-sample"]
    ext = "-gemini.db"
    out = []
    for fname in glob.glob(os.path.join(work_dir, "*{0}".format(ext))):
        cur = {"name": os.path.basename(fname).replace(ext, "")}
        for stat in stats:
            cur[stat] = calculate_gemini_stat(stat, fname)
        out.append(cur)
    df = pandas.DataFrame(out)
    print df

def calculate_gemini_stat(stat, fname):
    out = sh.gemini("stats", "--{0}".format(stat), fname)
    return float(out.split("\n")[1].split("\t")[-1])

if __name__ == "__main__":
    main(os.getcwd())
    
