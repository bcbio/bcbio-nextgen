"""Extract variant statistics from GEMINI database.
"""
import os
import subprocess

import toolz as tz
import yaml

from bcbio import utils
from bcbio.provenance import do
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils

def run(bam_file, data, out_dir):
    """Retrieve high level variant statistics from Gemini.
    """
    out = {}
    gemini_dbs = [d for d in
                  [tz.get_in(["population", "db"], x) for x in data.get("variants", [])] if d]
    if len(gemini_dbs) > 0:
        out_dir = utils.safe_makedir(out_dir)
        gemini_db = gemini_dbs[0]
        gemini_stat_file = os.path.join(out_dir, "%s-stats.yaml" % os.path.splitext(os.path.basename(gemini_db))[0])
        if not utils.file_uptodate(gemini_stat_file, gemini_db):
            gemini = config_utils.get_program("gemini", data["config"])
            tstv = subprocess.check_output([gemini, "stats", "--tstv", gemini_db])
            gt_counts = subprocess.check_output([gemini, "stats", "--gts-by-sample", gemini_db])
            dbsnp_count = subprocess.check_output([gemini, "query", gemini_db, "-q",
                                                   "SELECT count(*) FROM variants WHERE in_dbsnp==1"])
            out["Transition/Transversion"] = tstv.split("\n")[1].split()[-1]
            for line in gt_counts.split("\n"):
                parts = line.rstrip().split()
                if len(parts) > 0 and parts[0] != "sample":
                    name, hom_ref, het, hom_var, _, total = parts
                    out[name] = {}
                    out[name]["Variations (heterozygous)"] = int(het)
                    out[name]["Variations (homozygous)"] = int(hom_var)
                    # same total variations for all samples, keep that top level as well.
                    out["Variations (total)"] = int(total)
            out["Variations (in dbSNP)"] = int(dbsnp_count.strip())
            if out.get("Variations (total)") > 0:
                out["Variations (in dbSNP) pct"] = "%.1f%%" % (out["Variations (in dbSNP)"] /
                                                               float(out["Variations (total)"]) * 100.0)
            with open(gemini_stat_file, "w") as out_handle:
                yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
        else:
            with open(gemini_stat_file) as in_handle:
                out = yaml.safe_load(in_handle)

    res = {}
    for k, v in out.iteritems():
        if not isinstance(v, dict):
            res.update({k: v})
        if k == dd.get_sample_name(data):
            res.update(v)
    return res
