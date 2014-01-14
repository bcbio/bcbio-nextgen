"""Calculate potential effects of variations using external programs.

Supported:
  snpEff: http://sourceforge.net/projects/snpeff/
"""
import os
import csv
import glob

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils, tools
from bcbio.provenance import do
from bcbio.variation import vcfutils

# ## snpEff variant effects

def snpeff_effects(data):
    """Annotate input VCF file with effects calculated by snpEff.
    """
    vcf_in = data["vrn_file"]
    interval_file = data["config"]["algorithm"].get("variant_regions", None)
    if vcfutils.vcf_has_variants(vcf_in):
        se_interval = (_convert_to_snpeff_interval(interval_file, vcf_in)
                       if interval_file else None)
        try:
            vcf_file = _run_snpeff(vcf_in, se_interval, "vcf", data)
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
    resources = config_utils.get_resources("snpeff", config)
    if resources.get("options"):
        args += [str(x) for x in resources.get("options", [])]
    # cancer specific calling arguments
    if data.get("metadata", {}).get("phenotype") in ["tumor", "normal"]:
        args += ["-cancer"]
    # Provide options tuned to reporting variants in clinical environments
    if config["algorithm"].get("clinical_reporting"):
        args += ["-canon", "-hgvs"]
    return args

def get_db(ref_file, resources, config=None):
    """Retrieve a snpEff database name and location relative to reference file.
    """
    snpeff_db = resources.get("aliases", {}).get("snpeff")
    if snpeff_db:
        snpeff_base_dir = utils.safe_makedir(os.path.normpath(os.path.join(
            os.path.dirname(os.path.dirname(ref_file)), "snpeff")))
        # back compatible retrieval of genome from installation directory
        if config and not os.path.exists(os.path.join(snpeff_base_dir, snpeff_db)):
            snpeff_base_dir, snpeff_db = _installed_snpeff_genome(snpeff_db, config)
    else:
        snpeff_base_dir = None
    return snpeff_db, snpeff_base_dir

def get_cmd(cmd_name, datadir, config):
    """Retrieve snpEff base command line, handling command line and jar based installs.
    """
    resources = config_utils.get_resources("snpeff", config)
    memory = " ".join(resources.get("jvm_opts", ["-Xms750m", "-Xmx5g"]))
    try:
        snpeff = config_utils.get_program("snpeff", config)
        cmd = "{snpeff} {memory} {cmd_name} -dataDir {datadir}"
    except config_utils.CmdNotFound:
        snpeff_jar = config_utils.get_jar("snpEff",
                                          config_utils.get_program("snpeff", config, "dir"))
        config_file = "%s.config" % os.path.splitext(snpeff_jar)[0]
        cmd = "java -jar {snpeff_jar} {cmd_name} -c {config_file} -dataDir {datadir}"
    return cmd.format(**locals())

def _run_snpeff(snp_in, se_interval, out_format, data):
    snpeff_db, datadir = get_db(data["sam_ref"], data["genome_resources"], data["config"])
    assert datadir is not None, \
        "Did not find snpEff resources in genome configuration: %s" % data["genome_resources"]
    assert os.path.exists(os.path.join(datadir, snpeff_db)), \
        "Did not find %s snpEff genome data in %s" % (snpeff_db, datadir)
    snpeff_cmd = get_cmd("eff", datadir, data["config"])
    ext = utils.splitext_plus(snp_in)[1] if out_format == "vcf" else ".tsv"
    out_file = "%s-effects%s" % (utils.splitext_plus(snp_in)[0], ext)
    if not utils.file_exists(out_file):
        interval = "-filterinterval %s" % (se_interval) if se_interval else ""
        config_args = " ".join(_snpeff_args_from_config(data))
        if ext.endswith(".gz"):
            bgzip_cmd = "| %s -c" % tools.get_bgzip_cmd(data["config"])
        else:
            bgzip_cmd = ""
        with file_transaction(out_file) as tx_out_file:
            cmd = ("{snpeff_cmd} {interval} {config_args} -noLog -1 -i vcf -o {out_format} "
                   "{snpeff_db} {snp_in} {bgzip_cmd} > {tx_out_file}")
            do.run(cmd.format(**locals()), "snpEff effects", data)
    if ext.endswith(".gz"):
        out_file = vcfutils.bgzip_and_index(out_file, data["config"])
    return out_file

def _convert_to_snpeff_interval(in_file, base_file):
    """Handle wide variety of BED-like inputs, converting to BED-3.
    """
    out_file = "%s-snpeff-intervals.bed" % utils.splitext_plus(base_file)[0]
    if not os.path.exists(out_file):
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            with open(in_file) as in_handle:
                for line in (l for l in in_handle if not l.startswith(("@", "#"))):
                    parts = line.split()
                    writer.writerow(parts[:3])
    return out_file

# ## back-compatibility

def _find_snpeff_datadir(config_file):
    with open(config_file) as in_handle:
        for line in in_handle:
            if line.startswith("data_dir"):
                data_dir = config_utils.expand_path(line.split("=")[-1].strip())
                if not data_dir.startswith("/"):
                    data_dir = os.path.join(os.path.dirname(config_file), data_dir)
                return data_dir
    raise ValueError("Did not find data directory in snpEff config file: %s" % config_file)

def _installed_snpeff_genome(base_name, config):
    """Find the most recent installed genome for snpEff with the given name.
    """
    snpeff_config_file = os.path.join(config_utils.get_program("snpEff", config, "dir"),
                                      "snpEff.config")
    data_dir = _find_snpeff_datadir(snpeff_config_file)
    dbs = [d for d in sorted(glob.glob(os.path.join(data_dir, "%s*" % base_name)), reverse=True)
           if os.path.isdir(d)]
    if len(dbs) == 0:
        raise ValueError("No database found in %s for %s" % (data_dir, base_name))
    else:
        return data_dir, os.path.split(dbs[0])[-1]
