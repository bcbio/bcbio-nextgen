"""Calculate potential effects of variations using external programs.

Supported:
  snpEff: http://sourceforge.net/projects/snpeff/
  VEP: http://www.ensembl.org/info/docs/tools/vep/index.html
"""
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
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do, programs
from bcbio.variation import vcfutils

# ## High level

def add_to_vcf(in_file, data, effect_todo=None):
    if effect_todo is None:
        effect_todo = get_type(data)
    if effect_todo:
        stats = None
        if effect_todo == "snpeff":
            ann_vrn_file, stats_files = snpeff_effects(in_file, data)
            if stats_files:
                stats = {}
                for key, val in zip(["effects-stats", "effects-stats-csv"], stats_files):
                    if utils.file_exists(val):
                        stats[key] = val
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
                cmd = ["%svep_install" % vep_path, "-a", "c", "-s", ensembl_name,
                       "-c", vep_dir, "-u", tmp_dir]
                do.run("%s && %s" % (perl_exports, " ".join(cmd)), "Prepare VEP directory for %s" % ensembl_name)
                cmd = ["%svep_convert_cache" % vep_path, "-species", species, "-version", vepv,
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
                vep = config_utils.get_program("vep", data["config"])
                is_human = tz.get_in(["genome_resources", "aliases", "human"], data, False)
                config_args = []
                if is_human:
                    plugin_fns = { "loftee": _get_loftee, "maxentscan": _get_maxentscan, "genesplicer": _get_genesplicer}
                    plugins = ["loftee"]
                    if "vep_splicesite_annotations" in dd.get_tools_on(data):
                        plugins += ["maxentscan","genesplicer"]
                    for plugin in plugins:
                        plugin_args = plugin_fns[plugin](data)
                        config_args += plugin_args
                    config_args += ["--sift", "b", "--polyphen", "b"]
                    # Use HGVS by default, requires indexing the reference genome
                    config_args += ["--hgvs", "--shift_hgvs", "1", "--fasta", dd.get_ref_file(data)]
                if (dd.get_effects_transcripts(data).startswith("canonical")
                      or tz.get_in(("config", "algorithm", "clinical_reporting"), data)):
                    config_args += ["--pick"]
                resources = config_utils.get_resources("vep", data["config"])
                extra_args = [str(x) for x in resources.get("options", [])]
                cmd = [vep, "--vcf", "-o", "stdout", "-i", in_file] + fork_args + extra_args + \
                      ["--species", ensembl_name,
                       "--no_stats", "--cache", "--merged"
                        "--offline", "--dir", vep_dir,
                       "--symbol", "--numbers", "--biotype", "--total_length", "--canonical",
                       "--gene_phenotype", "--ccds", "--uniprot", "--domains", "--regulatory",
                       "--protein", "--tsl", "--appris", "--af", "--max_af" "--af_1kg", "--af_esp", "--af_exac",
                       "--pubmed", "--variant_class"] + config_args
                perl_exports = utils.get_perl_exports()
                # Remove empty fields (';;') which can cause parsing errors downstream
                cmd = "%s && %s | sed '/^#/! s/;;/;/g' | bgzip -c > %s" % (perl_exports, " ".join(cmd), tx_out_file)
                do.run(cmd, "Ensembl variant effect predictor", data)
    if utils.file_exists(out_file):
        vcfutils.bgzip_and_index(out_file, data["config"])
        return out_file

def _get_loftee(data):
    """Retrieve loss of function plugin parameters for LOFTEE.
    https://github.com/konradjk/loftee
    """
    ancestral_file = tz.get_in(("genome_resources", "variation", "ancestral"), data)
    if not ancestral_file or not os.path.exists(ancestral_file):
        ancestral_file = "false"
    args = ["--plugin", "LoF,human_ancestor_fa:%s" % ancestral_file]
    return args

def _get_maxentscan(data):
    """
    The plugin executes the logic from one of the scripts depending on which
    splice region the variant overlaps:
        score5.pl : last 3 bases of exon    --> first 6 bases of intron
        score3.pl : last 20 bases of intron --> first 3 bases of exon
    The plugin reports the reference, alternate and difference (REF - ALT) maximumentropy scores.
    https://github.com/Ensembl/VEP_plugins/blob/master/MaxEntScan.pm
    """

    maxentscan_dir = os.path.dirname(os.path.realpath(config_utils.get_program("maxentscan_score3.pl", data["config"])))
    if maxentscan_dir and os.path.exists(maxentscan_dir):
        return ["--plugin", "MaxEntScan,%s" % (maxentscan_dir)]
    else:
        return []

def _get_genesplicer(data):
    """
    This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
    runs GeneSplicer (https://ccb.jhu.edu/software/genesplicer/) to get splice site predictions.
    https://github.com/Ensembl/VEP_plugins/blob/master/GeneSplicer.pm
    """

    genesplicer_exec = os.path.realpath(config_utils.get_program("genesplicer", data["config"]))
    genesplicer_training = tz.get_in(("genome_resources", "variation", "genesplicer"), data)
    if genesplicer_exec and os.path.exists(genesplicer_exec) and genesplicer_training and os.path.exists(genesplicer_training) :
        return ["--plugin", "GeneSplicer,%s,%s" % (genesplicer_exec,genesplicer_training)]
    else:
        return []

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
    args = ["-hgvs"]
    # General supplied arguments
    resources = config_utils.get_resources("snpeff", config)
    if resources.get("options"):
        args += [str(x) for x in resources.get("options", [])]
    # cancer specific calling arguments
    if vcfutils.get_paired_phenotype(data):
        args += ["-cancer"]

    effects_transcripts = dd.get_effects_transcripts(data)
    if effects_transcripts in set(["canonical_cancer"]):
        _, snpeff_base_dir = get_db(data)
        canon_list_file = os.path.join(snpeff_base_dir, "transcripts", "%s.txt" % effects_transcripts)
        if not utils.file_exists(canon_list_file):
            raise ValueError("Cannot find expected file for effects_transcripts: %s" % canon_list_file)
        args += ["-canonList", canon_list_file]
    elif effects_transcripts == "canonical" or tz.get_in(("config", "algorithm", "clinical_reporting"), data):
        args += ["-canon"]
    return args

def get_db(data):
    """Retrieve a snpEff database name and location relative to reference file.
    """
    snpeff_db = utils.get_in(data, ("genome_resources", "aliases", "snpeff"))
    snpeff_base_dir = None
    if snpeff_db:
        snpeff_base_dir = utils.get_in(data, ("reference", "snpeff", snpeff_db, "base"))
        if not snpeff_base_dir:
            # We need to mask '.' characters for CWL/WDL processing, check for them here
            snpeff_base_dir = utils.get_in(data, ("reference", "snpeff", snpeff_db.replace("_", "."), "base"))
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
        # Clean problem characters for CWL/WDL representation
        clean_snpeff_db = snpeff_db.replace(".", "_")
        return {clean_snpeff_db: {"base": datadir,
                                  "indexes": glob.glob(os.path.join(datadir, snpeff_db, "*"))}}
    else:
        return {}

def _get_snpeff_cmd(cmd_name, datadir, data, out_file):
    """Retrieve snpEff base command line.
    """
    resources = config_utils.get_resources("snpeff", data["config"])
    jvm_opts = resources.get("jvm_opts", ["-Xms750m", "-Xmx3g"])
    # scale by cores, defaulting to 2x base usage to ensure we have enough memory
    # for single core runs to use with human genomes
    jvm_opts = config_utils.adjust_opts(jvm_opts, {"algorithm": {"memory_adjust":
                                                                 {"direction": "increase",
                                                                  "magnitude": max(2, dd.get_cores(data))}}})
    memory = " ".join(jvm_opts)
    snpeff = config_utils.get_program("snpEff", data["config"])
    java_args = "-Djava.io.tmpdir=%s" % utils.safe_makedir(os.path.join(os.path.dirname(out_file), "tmp"))
    export = utils.local_path_export()
    cmd = "{export} {snpeff} {memory} {java_args} {cmd_name} -dataDir {datadir}"
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
    csv_file = "%s-stats.csv" % utils.splitext_plus(out_file)[0]
    if not utils.file_exists(out_file):
        config_args = " ".join(_snpeff_args_from_config(data))
        if ext.endswith(".gz"):
            bgzip_cmd = "| %s -c" % tools.get_bgzip_cmd(data["config"])
        else:
            bgzip_cmd = ""
        with file_transaction(data, out_file) as tx_out_file:
            snpeff_cmd = _get_snpeff_cmd("eff", datadir, data, tx_out_file)
            cmd = ("{snpeff_cmd} {config_args} -noLog -i vcf -o {out_format} "
                   "-csvStats {csv_file} -s {stats_file} {snpeff_db} {snp_in} {bgzip_cmd} > {tx_out_file}")
            do.run(cmd.format(**locals()), "snpEff effects", data)
    if ext.endswith(".gz"):
        out_file = vcfutils.bgzip_and_index(out_file, data["config"])
    return out_file, [stats_file, csv_file]

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
