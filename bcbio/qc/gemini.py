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
        name = dd.get_sample_name(data)
        out_dir = utils.safe_makedir(out_dir)
        gemini_db = gemini_dbs[0]
        gemini_stat_file = os.path.join(out_dir, "%s-%s-stats.yaml" % (os.path.splitext(os.path.basename(gemini_db))[0], name))
        if not utils.file_uptodate(gemini_stat_file, gemini_db):
            gemini = config_utils.get_program("gemini", data["config"])
            cmd = [gemini, "query", gemini_db, "-q",
                   "SELECT count(*) FROM variants",
                   "--gt-filter",
                   "gt_types.%s != HOM_REF" % name]
            gt_counts = subprocess.check_output(cmd)
            cmd = [gemini, "query", gemini_db, "-q",
                   "SELECT count(*) FROM variants WHERE in_dbsnp==1",
                   "--gt-filter", "gt_types.%s != HOM_REF" % name]
            dbsnp_counts = subprocess.check_output(cmd)
            cmd = [gemini, "query", gemini_db, "-q",
                   "SELECT count(*) FROM variants ",
                   "--gt-filter", "gt_types.%s == HET" % name]
            het_counts = subprocess.check_output(cmd)
            cmd = [gemini, "query", gemini_db, "-q",
                   "SELECT count(*) FROM variants ",
                   "--gt-filter", "gt_types.%s == HOM_ALT" % name]
            hom_alt_counts = subprocess.check_output(cmd)
            out["Variations (heterozygous)"] = int(het_counts.strip()) if het_counts else 0
            out["Variations (homozygous)"] = int(hom_alt_counts.strip()) if hom_alt_counts else 0
            out["Variations (total)"] = int(gt_counts.strip()) if gt_counts else 0
            out["Variations (in dbSNP)"] = int(dbsnp_counts.strip()) if dbsnp_counts else 0
            if out.get("Variations (total)") > 0:
                out["Variations (in dbSNP) pct"] = "%.1f%%" % (out["Variations (in dbSNP)"] /
                                                               float(out["Variations (total)"]) * 100.0)
            with open(gemini_stat_file, "w") as out_handle:
                yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
        else:
            with open(gemini_stat_file) as in_handle:
                out = yaml.safe_load(in_handle)
    return out
