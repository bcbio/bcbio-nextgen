"""Calculate potential effects of variations using external programs.

Supported:
  snpEff: http://sourceforge.net/projects/snpeff/
  VEP: http://www.ensembl.org/info/docs/tools/vep/index.html
"""
from distutils.version import LooseVersion
import os
import glob
import shutil
import subprocess
import string

import toolz as tz
import yaml

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils, tools
from bcbio.provenance import do, programs
from bcbio.variation import vcfutils

# ## High level

def add_to_vcf(in_file, data, effect_todo=None):
    if effect_todo is None:
        effect_todo = get_type(data)
    if effect_todo:
        stats = None
        if effect_todo == "snpeff":
            ann_vrn_file, stats_file = snpeff_effects(in_file, data)
            if utils.file_exists(stats_file):
                stats = {"effects-stats": stats_file}
        elif effect_todo == "vep":
            ann_vrn_file = run_vep(in_file, data)
        else:
            raise ValueError("Unexpected variant effects configuration: %s" % effect_todo)
        if ann_vrn_file:
            return ann_vrn_file, stats
    return None, None

def get_type(data):
    """Retrieve the type of effects calculation to do.
    """
    if data["analysis"].lower().startswith("var"):
        return tz.get_in(("config", "algorithm", "effects"), data, "snpeff")

# ## Ensembl VEP

def _special_dbkey_maps(dbkey, ref_file):
    """Avoid duplicate VEP information for databases with chromosome differences like hg19/GRCh37.
    """
    remaps = {"hg19": "GRCh37",
              "hg38-noalt": "hg38"}
    if dbkey in remaps:
        base_dir = os.path.normpath(os.path.join(os.path.dirname(ref_file), os.pardir))
        vep_dir = os.path.normpath(os.path.join(base_dir, "vep"))
        other_dir = os.path.relpath(os.path.normpath(os.path.join(base_dir, os.pardir, remaps[dbkey], "vep")),
                                    base_dir)
        if os.path.exists(os.path.join(base_dir, other_dir)):
            if not os.path.lexists(vep_dir):
                os.symlink(other_dir, vep_dir)
            return vep_dir
        else:
            return None
    else:
        return None


def prep_vep_cache(dbkey, ref_file, tooldir=None, config=None):
    """Ensure correct installation of VEP cache file.
    """
    if config is None: config = {}
    resource_file = os.path.join(os.path.dirname(ref_file), "%s-resources.yaml" % dbkey)
    if os.path.exists(resource_file):
        with open(resource_file) as in_handle:
            resources = yaml.load(in_handle)
        ensembl_name = tz.get_in(["aliases", "ensembl"], resources)
        symlink_dir = _special_dbkey_maps(dbkey, ref_file)
        if ensembl_name and ensembl_name.find("_vep_") == -1:
            raise ValueError("%s has ensembl an incorrect value."
                             "It should have _vep_ in the name."
                             "Remove line or fix the name to avoid error.")
        if symlink_dir and ensembl_name:
            species, vepv = ensembl_name.split("_vep_")
            return symlink_dir, species
        elif ensembl_name:
            species, vepv = ensembl_name.split("_vep_")
            vep_dir = utils.safe_makedir(os.path.normpath(os.path.join(
                os.path.dirname(os.path.dirname(ref_file)), "vep")))
            out_dir = os.path.join(vep_dir, species, vepv)
            if not os.path.exists(out_dir):
                tmp_dir = utils.safe_makedir(os.path.join(vep_dir, species, "txtmp"))
                eversion = vepv.split("_")[0]
                url = "ftp://ftp.ensembl.org/pub/release-%s/variation/VEP/%s.tar.gz" % (eversion, ensembl_name)
                with utils.chdir(tmp_dir):
                    subprocess.check_call(["wget", "--no-check-certificate", "-c", url])
                vep_path = "%s/bin/" % tooldir if tooldir else ""
                perl_exports = utils.get_perl_exports()
                cmd = ["%svep_install.pl" % vep_path, "-a", "c", "-s", ensembl_name,
                       "-c", vep_dir, "-u", tmp_dir]
                do.run("%s && %s" % (perl_exports, " ".join(cmd)), "Prepare VEP directory for %s" % ensembl_name)
                cmd = ["%svep_convert_cache.pl" % vep_path, "-species", species, "-version", vepv,
                       "-d", vep_dir]
                do.run("%s && %s" % (perl_exports, " ".join(cmd)), "Convert VEP cache to tabix %s" % ensembl_name)
                for tmp_fname in os.listdir(tmp_dir):
                    os.remove(os.path.join(tmp_dir, tmp_fname))
                os.rmdir(tmp_dir)
            tmp_dir = os.path.join(vep_dir, "tmp")
            if os.path.exists(tmp_dir):
                shutil.rmtree(tmp_dir)
            return vep_dir, species
    return None, None

def run_vep(in_file, data):
    """Annotate input VCF file with Ensembl variant effect predictor.
    """
    if not vcfutils.vcf_has_variants(in_file):
        return None
    out_file = utils.append_stem(in_file, "-vepeffects")
    assert in_file.endswith(".gz") and out_file.endswith(".gz")
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            vep_dir, ensembl_name = prep_vep_cache(data["genome_build"],
                                                   tz.get_in(["reference", "fasta", "base"], data))
            if vep_dir:
                cores = tz.get_in(("config", "algorithm", "num_cores"), data, 1)
                fork_args = ["--fork", str(cores)] if cores > 1 else []
                vep = config_utils.get_program("variant_effect_predictor.pl", data["config"])
                is_human = tz.get_in(["genome_resources", "aliases", "human"], data, False)
                if is_human:
                    dbnsfp_args, dbnsfp_fields = _get_dbnsfp(data)
                    loftee_args, loftee_fields = _get_loftee(data)
                    prediction_args = ["--sift", "b", "--polyphen", "b"]
                    prediction_fields = ["PolyPhen", "SIFT"]
                else:
                    dbnsfp_args, dbnsfp_fields = [], []
                    loftee_args, loftee_fields = [], []
                    prediction_args, prediction_fields = [], []
                std_fields = ["Consequence", "Codons", "Amino_acids", "Gene", "SYMBOL", "Feature",
                              "EXON"] + prediction_fields + ["Protein_position", "BIOTYPE", "CANONICAL", "CCDS"]
                resources = config_utils.get_resources("vep", data["config"])
                extra_args = [str(x) for x in resources.get("options", [])]
                cmd = [vep, "--vcf", "-o", "stdout", "-i", in_file] + fork_args + extra_args + \
                      ["--species", ensembl_name,
                       "--no_stats",
                       "--cache", "--offline", "--dir", vep_dir,
                       "--symbol", "--numbers", "--biotype", "--total_length", "--canonical", "--ccds",
                       "--fields", ",".join(std_fields + dbnsfp_fields + loftee_fields)] + \
                       prediction_args + dbnsfp_args + loftee_args

                if tz.get_in(("config", "algorithm", "clinical_reporting"), data, False):

                    # In case of clinical reporting, we need one and only one
                    # variant per gene
                    # From the VEP docs:
                    # "Pick once line of consequence data per variant,
                    # including transcript-specific columns. Consequences are
                    # chosen by the canonical, biotype status and length of the
                    # transcript, along with the ranking of the consequence
                    # type according to this table. This is the best method to
                    # use if you are interested only in one consequence per
                    #  variant.

                    cmd += ["--pick"]

                    # TODO investigate hgvs reporting but requires indexing the reference file
                    # cmd += ["--hgvs", "--shift-hgvs", "--fasta", dd.get_ref_file(data)]
                perl_exports = utils.get_perl_exports()
                # Remove empty fields (';;') which can cause parsing errors downstream
                cmd = "%s && %s | sed '/^#/! s/;;/;/g' | bgzip -c > %s" % (perl_exports, " ".join(cmd), tx_out_file)
                do.run(cmd, "Ensembl variant effect predictor", data)
    if utils.file_exists(out_file):
        vcfutils.bgzip_and_index(out_file, data["config"])
        return out_file

def _get_dbnsfp(data):
    """Retrieve dbNSFP file options for VEP if downloaded and available.

    Uses high level combined annotations from this GEMINI discussion as a
    starting point:
    https://groups.google.com/d/msg/gemini-variation/WeZ6C2YvfUA/mII9uum_pGoJ
    """
    dbnsfp_file = tz.get_in(("genome_resources", "variation", "dbnsfp"), data)
    if dbnsfp_file and os.path.exists(dbnsfp_file):
        annotations = ["RadialSVM_score", "RadialSVM_pred", "LR_score", "LR_pred",
                       "CADD_raw", "CADD_phred", "Reliability_index"]
        return ["--plugin", "dbNSFP,%s,%s" % (dbnsfp_file, ",".join(annotations))], annotations
    else:
        return [], []

def _get_loftee(data):
    """Retrieve loss of function plugin parameters for LOFTEE.
    https://github.com/konradjk/loftee
    """
    ancestral_file = tz.get_in(("genome_resources", "variation", "ancestral"), data)
    if not ancestral_file or not os.path.exists(ancestral_file):
        ancestral_file = "false"
    annotations = ["LoF", "LoF_filter", "LoF_flags"]
    args = ["--plugin", "LoF,human_ancestor_fa:%s" % ancestral_file]
    return args, annotations

# ## snpEff variant effects

def snpeff_version(args=None, data=None):
    raw_version = programs.get_version_manifest("snpeff", data=data)
    if not raw_version:
        raw_version = ""
    snpeff_version = "".join([x for x in str(raw_version)
                              if x in set(string.digits + ".")])
    return snpeff_version

def snpeff_effects(vcf_in, data):
    """Annotate input VCF file with effects calculated by snpEff.
    """
    if vcfutils.vcf_has_variants(vcf_in):
        return _run_snpeff(vcf_in, "vcf", data)
    else:
        return None, None

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

    # Skip HGVS if running structural variant calling due to errors
    # https://github.com/chapmanb/bcbio-nextgen/issues/1205
    # https://github.com/pcingola/SnpEff/issues/128
    svcaller = tz.get_in(["config", "algorithm", "svcaller_active"], data)
    if svcaller:
        args += ["-noHgvs"]
    # Provide options tuned to reporting variants in clinical environments
    elif config["algorithm"].get("clinical_reporting"):
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

def get_cmd(cmd_name, datadir, config, out_file):
    """Retrieve snpEff base command line.
    """
    resources = config_utils.get_resources("snpeff", config)
    memory = " ".join(resources.get("jvm_opts", ["-Xms750m", "-Xmx5g"]))
    snpeff = config_utils.get_program("snpEff", config)
    java_args = "-Djava.io.tmpdir=%s" % utils.safe_makedir(os.path.join(os.path.dirname(out_file), "tmp"))
    cmd = "{snpeff} {memory} {java_args} {cmd_name} -dataDir {datadir}"
    return cmd.format(**locals())

def _run_snpeff(snp_in, out_format, data):
    """Run effects prediction with snpEff, skipping if snpEff database not present.
    """
    snpeff_db, datadir = get_db(data)
    if not snpeff_db:
        return None, None

    assert os.path.exists(os.path.join(datadir, snpeff_db)), \
        "Did not find %s snpEff genome data in %s" % (snpeff_db, datadir)
    ext = utils.splitext_plus(snp_in)[1] if out_format == "vcf" else ".tsv"
    out_file = "%s-effects%s" % (utils.splitext_plus(snp_in)[0], ext)
    stats_file = "%s-stats.html" % utils.splitext_plus(out_file)[0]
    if not utils.file_exists(out_file):
        config_args = " ".join(_snpeff_args_from_config(data))
        if ext.endswith(".gz"):
            bgzip_cmd = "| %s -c" % tools.get_bgzip_cmd(data["config"])
        else:
            bgzip_cmd = ""
        with file_transaction(data, out_file) as tx_out_file:
            snpeff_cmd = get_cmd("eff", datadir, data["config"], tx_out_file)
            cmd = ("{snpeff_cmd} {config_args} -noLog -i vcf -o {out_format} "
                   "-s {stats_file} {snpeff_db} {snp_in} {bgzip_cmd} > {tx_out_file}")
            do.run(cmd.format(**locals()), "snpEff effects", data)
    if ext.endswith(".gz"):
        out_file = vcfutils.bgzip_and_index(out_file, data["config"])
    return out_file, stats_file

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
