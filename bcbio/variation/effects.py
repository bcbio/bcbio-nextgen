"""Calculate potential effects of variations using external programs.

Supported:
  snpEff: http://sourceforge.net/projects/snpeff/
"""
import os
import csv
import subprocess

from bcbio.utils import file_transaction

# ## snpEff variant effects

# remap Galaxy genome names to the ones used by snpEff. Not nice code.
SNPEFF_GENOME_REMAP = {
        "GRCh37": "GRCh37.63",
        "hg19" : "GRCh37.63",
        "mm9" : "NCBIM37.63",
        "araTha_tair9": "athalianaTair10",
        "araTha_tair10": "athalianaTair10",
        }

def snpeff_effects(vcf_in, genome, config):
    """Prepare tab-delimited file for variant effects using snpEff.
    """
    interval_file = config["algorithm"].get("hybrid_target", None)
    if _vcf_has_items(vcf_in):
        se_interval = (_convert_to_snpeff_interval(interval_file, vcf_in)
                       if interval_file else None)
        try:
            genome = SNPEFF_GENOME_REMAP[genome]
            vcf_file = _run_snpeff(vcf_in, genome, se_interval, "vcf", config)
            effects_file = _run_snpeff(vcf_in, genome, se_interval, "txt", config)
        finally:
            for fname in [se_interval]:
                if fname and os.path.exists(fname):
                    os.remove(fname)
        return vcf_file, effects_file
    else:
        return None, None

def _run_snpeff(snp_in, genome, se_interval, out_format, config):
    snpeff_jar = os.path.join(config["program"]["snpEff"], "snpEff.jar")
    java_memory = config["algorithm"].get("java_memory", None)
    snpeff_config = "%s.config" % os.path.splitext(snpeff_jar)[0]
    ext = "vcf" if out_format == "vcf" else "tsv"
    out_file = "%s-effects.%s" % (os.path.splitext(snp_in)[0], ext)
    if not os.path.exists(out_file):
        cl = ["java"]
        if java_memory:
            cl += ["-Xmx%s" % java_memory]
        cl += ["-jar", snpeff_jar, "-1", "-i", "vcf", "-c", snpeff_config,
               "-o", out_format, genome, snp_in]
        if se_interval:
            cl.extend(["-filterInterval", se_interval])
        print " ".join(cl)
        with file_transaction(out_file):
            with open(out_file, "w") as out_handle:
                subprocess.check_call(cl, stdout=out_handle)
    return out_file

def _vcf_has_items(in_file):
    if os.path.exists(in_file):
        with open(in_file) as in_handle:
            for line in in_handle:
                if line.strip() and not line.startswith("#"):
                    return True
    return False

def _convert_to_snpeff_interval(in_file, base_file):
    """Handle wide variety of BED-like inputs, converting to BED-3.
    """
    out_file = "%s-snpeff-intervals.bed" % os.path.splitext(base_file)[0]
    if not os.path.exists(out_file):
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            with open(in_file) as in_handle:
                for line in (l for l in in_handle if not l.startswith(("@", "#"))):
                    parts = line.split()
                    writer.writerow(parts[:3])
    return out_file
