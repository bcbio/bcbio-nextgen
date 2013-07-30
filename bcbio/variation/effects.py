"""Calculate potential effects of variations using external programs.

Supported:
  snpEff: http://sourceforge.net/projects/snpeff/
"""
import os
import csv
import glob
import subprocess
import collections

from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils

# ## snpEff variant effects

# remap Galaxy genome names to the ones used by snpEff. Not nice code.
SnpEffGenome = collections.namedtuple("SnpEffGenome", ["base", "default_version", "is_human"])
SNPEFF_GENOME_REMAP = {
        "GRCh37": SnpEffGenome("GRCh37.", "68", True),
        "b37": SnpEffGenome("GRCh37.", "68", True),
        "hg19" : SnpEffGenome("hg19", "", True),
        "mm9" : SnpEffGenome("NCBIM37.", "68", False),
        "mm10" : SnpEffGenome("GRCm38.", "71", False),
        "araTha_tair9": SnpEffGenome("athalianaTair9", "", False),
        "araTha_tair10": SnpEffGenome("athalianaTair10", "", False),
        }

def _find_snpeff_datadir(config_file):
    with open(config_file) as in_handle:
        for line in in_handle:
            if line.startswith("data_dir"):
                data_dir = config_utils.expand_path(line.split("=")[-1].strip())
                if not data_dir.startswith("/"):
                    data_dir = os.path.join(os.path.dirname(config_file), data_dir)
                return data_dir
    raise ValueError("Did not find data directory in snpEff config file: %s" % config_file)

def _installed_snpeff_genome(config_file, base_name):
    """Find the most recent installed genome for snpEff with the given name.
    """
    data_dir = _find_snpeff_datadir(config_file)
    dbs = [d for d in sorted(glob.glob(os.path.join(data_dir, "%s*" % base_name)), reverse=True)
           if os.path.isdir(d)]
    if len(dbs) == 0:
        raise ValueError("No database found in %s for %s" % (data_dir, base_name))
    else:
        return os.path.split(dbs[0])[-1]

def _get_snpeff_genome(gname, config):
    """Generalize retrieval of the snpEff genome to use for an input name.

    This tries to find the snpEff configuration file and identify the
    installed genome corresponding to the input genome name.
    """
    try:
        ginfo = SNPEFF_GENOME_REMAP[gname]
    except KeyError:
        ginfo = SNPEFF_GENOME_REMAP[gname.split("-")[0]]
    snpeff_config_file = os.path.join(config_utils.get_program("snpEff", config, "dir"),
                                      "snpEff.config")
    if os.path.exists(snpeff_config_file):
        return _installed_snpeff_genome(snpeff_config_file, ginfo.base)
    else:
        return "%s%s" % (ginfo.base, ginfo.default_version)

def snpeff_effects(data):
    """Annotate input VCF file with effects calculated by snpEff.
    """
    vcf_in = data["vrn_file"]
    interval_file = data["config"]["algorithm"].get("hybrid_target", None)
    if vcf_has_items(vcf_in):
        se_interval = (_convert_to_snpeff_interval(interval_file, vcf_in)
                       if interval_file else None)
        try:
            vcf_file = _run_snpeff(vcf_in, _get_snpeff_genome(data["genome_build"], data["config"]),
                                   se_interval, "vcf", data)
        finally:
            for fname in [se_interval]:
                if fname and os.path.exists(fname):
                    os.remove(fname)
        return vcf_file

def _snpeff_args_from_config(data):
    """Retrieve snpEff arguments supplied through input configuration.
    """
    config = data["config"]
    args = []
    # General supplied arguments
    resources = config_utils.get_resources("snpEff", config)
    if resources.get("options"):
        args += [str(x) for x in resources.get("options", [])]
    # cancer specific calling arguments
    if data.get("metadata", {}).get("phenotype") in ["tumor", "normal"]:
        args += ["-cancer"]
    # Provide options tuned to reporting variants in clinical environments
    if config["algorithm"].get("clinical_reporting"):
        args += ["-canon", "-hgvs"]
    return args

def _run_snpeff(snp_in, genome, se_interval, out_format, data):
    config = data["config"]
    snpeff_jar = config_utils.get_jar("snpEff",
                                      config_utils.get_program("snpEff", config, "dir"))
    config_file = "%s.config" % os.path.splitext(snpeff_jar)[0]
    resources = config_utils.get_resources("snpEff", config)
    ext = "vcf" if out_format == "vcf" else "tsv"
    out_file = "%s-effects.%s" % (os.path.splitext(snp_in)[0], ext)
    if not file_exists(out_file):
        cl = ["java"]
        cl += resources.get("jvm_opts", ["-Xms750m", "-Xmx5g"])
        cl += ["-jar", snpeff_jar, "eff", "-c", config_file,
               "-noLog", "-1", "-i", "vcf", "-o", out_format, genome, snp_in]
        if se_interval:
            cl.extend(["-filterInterval", se_interval])
        cl += _snpeff_args_from_config(data)
        with file_transaction(out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                subprocess.check_call(cl, stdout=out_handle)
    return out_file

def vcf_has_items(in_file):
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
