"""Calculate potential effects of variations using external programs.

Supported:
  snpEff: http://sourceforge.net/projects/snpeff/
  VEP: http://www.ensembl.org/info/docs/tools/vep/index.html
"""
import os
import glob
import shutil
import subprocess

import toolz as tz
import yaml

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils, tools
from bcbio.provenance import do
from bcbio.variation import vcfutils

# ## Ensembl VEP

def vep_version(config):
    try:
        vep = config_utils.get_program("variant_effect_predictor.pl", config)
        help_str = subprocess.check_output([vep, "--help"])
        for line in help_str.split("\n"):
            if line.startswith("version"):
                return line.split()[-1].strip()
        return None
    except config_utils.CmdNotFound:
        return None
        return False

def _special_dbkey_maps(dbkey, ref_file):
    """Avoid duplicate VEP information for databases with chromosome differences like hg19/GRCh37.
    """
    remaps = {"hg19": "GRCh37"}
    if dbkey in remaps:
        base_dir = os.path.normpath(os.path.join(os.path.dirname(ref_file), os.pardir))
        vep_dir = os.path.normpath(os.path.join(base_dir, "vep"))
        other_dir = os.path.relpath(os.path.normpath(os.path.join(base_dir, os.pardir, remaps[dbkey], "vep")),
                                    base_dir)
        if not os.path.lexists(vep_dir):
            os.symlink(other_dir, vep_dir)
        return vep_dir
    else:
        return None

def prep_vep_cache(dbkey, ref_file, config=None):
    """Ensure correct installation of VEP cache file.
    """
    if config is None: config = {}
    resource_file = os.path.join(os.path.dirname(ref_file), "%s-resources.yaml" % dbkey)
    vepv = vep_version(config)
    if os.path.exists(resource_file) and vepv:
        with open(resource_file) as in_handle:
            resources = yaml.load(in_handle)
        ensembl_name = tz.get_in(["aliases", "ensembl"], resources)
        symlink_dir = _special_dbkey_maps(dbkey, ref_file)
        if symlink_dir:
            return symlink_dir, ensembl_name
        elif ensembl_name:
            vep_dir = utils.safe_makedir(os.path.normpath(os.path.join(
                os.path.dirname(os.path.dirname(ref_file)), "vep")))
            out_dir = os.path.join(vep_dir, ensembl_name, vepv)
            if not os.path.exists(out_dir):
                cmd = ["vep_install.pl", "-a", "c", "-s", ensembl_name,
                       "-c", vep_dir]
                do.run(cmd, "Prepare VEP directory for %s" % ensembl_name)
                cmd = ["vep_convert_cache.pl", "-species", ensembl_name, "-version", vepv,
                       "-d", vep_dir]
                do.run(cmd, "Convert VEP cache to tabix %s" % ensembl_name)
            tmp_dir = os.path.join(vep_dir, "tmp")
            if os.path.exists(tmp_dir):
                shutil.rmtree(tmp_dir)
            return vep_dir, ensembl_name
    return None, None

def run_vep(data):
    """Annotate input VCF file with Ensembl variant effect predictor.
    """
    out_file = utils.append_stem(data["vrn_file"], "-vepeffects")
    assert data["vrn_file"].endswith(".gz") and out_file.endswith(".gz")
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            vep_dir, ensembl_name = prep_vep_cache(data["genome_build"],
                                                   tz.get_in(["reference", "fasta", "base"], data))
            if vep_dir:
                vep = config_utils.get_program("variant_effect_predictor.pl", data["config"])
                cmd = [vep, "--vcf", "-o", "stdout", "--fork", "4",
                       "--species", ensembl_name,
                       "--no_stats",
                       "--cache", "--offline", "--dir", vep_dir,
                       "--sift", "b", "--polyphen", "b", "--symbol", "--numbers", "--biotype", "--total_length",
                       "--fields",
                       "Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE"]
                cmd = "gunzip -c %s | %s | bgzip -c > %s" % (data["vrn_file"], " ".join(cmd), tx_out_file)
                do.run(cmd, "Ensembl variant effect predictor", data)
    if utils.file_exists(out_file):
        vcfutils.bgzip_and_index(out_file, data["config"])
        return out_file

# ## snpEff variant effects

def snpeff_effects(data):
    """Annotate input VCF file with effects calculated by snpEff.
    """
    vcf_in = data["vrn_file"]
    if vcfutils.vcf_has_variants(vcf_in):
        vcf_file = _run_snpeff(vcf_in, "vcf", data)
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
    if vcfutils.get_paired_phenotype(data):
        args += ["-cancer"]
    # Provide options tuned to reporting variants in clinical environments
    if config["algorithm"].get("clinical_reporting"):
        args += ["-canon", "-hgvs"]
    return args

def get_db(data):
    """Retrieve a snpEff database name and location relative to reference file.
    """
    snpeff_db = utils.get_in(data, ("genome_resources", "aliases", "snpeff"))
    snpeff_base_dir = None
    if snpeff_db:
        snpeff_base_dir = utils.get_in(data, ("reference", "snpeff", snpeff_db, "base"))
        if not snpeff_base_dir:
            ref_file = utils.get_in(data, ("reference", "fasta", "base"))
            snpeff_base_dir = utils.safe_makedir(os.path.normpath(os.path.join(
                os.path.dirname(os.path.dirname(ref_file)), "snpeff")))
            # back compatible retrieval of genome from installation directory
            if "config" in data and not os.path.exists(os.path.join(snpeff_base_dir, snpeff_db)):
                snpeff_base_dir, snpeff_db = _installed_snpeff_genome(snpeff_db, data["config"])
    return snpeff_db, snpeff_base_dir

def get_snpeff_files(data):
    try:
        snpeff_db, datadir = get_db(data)
    except ValueError:
        snpeff_db = None
    if snpeff_db:
        return {snpeff_db: {"base": datadir,
                            "indexes": glob.glob(os.path.join(datadir, snpeff_db, "*"))}}
    else:
        return {}

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
        cmd = "java {memory} -jar {snpeff_jar} {cmd_name} -c {config_file} -dataDir {datadir}"
    return cmd.format(**locals())

def _run_snpeff(snp_in, out_format, data):
    snpeff_db, datadir = get_db(data)
    assert datadir is not None, \
        "Did not find snpEff resources in genome configuration: %s" % data["genome_resources"]
    assert os.path.exists(os.path.join(datadir, snpeff_db)), \
        "Did not find %s snpEff genome data in %s" % (snpeff_db, datadir)
    snpeff_cmd = get_cmd("eff", datadir, data["config"])
    ext = utils.splitext_plus(snp_in)[1] if out_format == "vcf" else ".tsv"
    out_file = "%s-effects%s" % (utils.splitext_plus(snp_in)[0], ext)
    if not utils.file_exists(out_file):
        config_args = " ".join(_snpeff_args_from_config(data))
        if ext.endswith(".gz"):
            bgzip_cmd = "| %s -c" % tools.get_bgzip_cmd(data["config"])
        else:
            bgzip_cmd = ""
        with file_transaction(out_file) as tx_out_file:
            cmd = ("{snpeff_cmd} {config_args} -noLog -1 -i vcf -o {out_format} "
                   "{snpeff_db} {snp_in} {bgzip_cmd} > {tx_out_file}")
            do.run(cmd.format(**locals()), "snpEff effects", data)
    if ext.endswith(".gz"):
        out_file = vcfutils.bgzip_and_index(out_file, data["config"])
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
    snpeff_config_file = os.path.join(config_utils.get_program("snpeff", config, "dir"),
                                      "snpEff.config")
    data_dir = _find_snpeff_datadir(snpeff_config_file)
    dbs = [d for d in sorted(glob.glob(os.path.join(data_dir, "%s*" % base_name)), reverse=True)
           if os.path.isdir(d)]
    if len(dbs) == 0:
        raise ValueError("No database found in %s for %s" % (data_dir, base_name))
    else:
        return data_dir, os.path.split(dbs[0])[-1]
