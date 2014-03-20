"""Provide infrastructure to allow exploration of variations within populations.

Uses the gemini framework (https://github.com/arq5x/gemini) to build SQLite
database of variations for query and evaluation.
"""
import collections
from distutils.version import LooseVersion
import glob
import os
import subprocess

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.provenance import do, programs
from bcbio.variation import vcfutils

def prep_gemini_db(fnames, call_info, samples, data):
    """Prepare a gemini database from VCF inputs prepared with snpEff.
    """
    out_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "gemini"))
    name, caller, is_batch = call_info
    gemini_db = os.path.join(out_dir, "%s-%s.db" % (name, caller))
    gemini_vcf = get_multisample_vcf(fnames, name, caller, data)
    use_gemini_quick = (do_db_build(samples, check_gemini=False) and
                        any(vcfutils.vcf_has_variants(f) for f in fnames))
    if not utils.file_exists(gemini_db) and use_gemini_quick:
        use_gemini = do_db_build(samples) and any(vcfutils.vcf_has_variants(f) for f in fnames)
        if use_gemini:
            with file_transaction(gemini_db) as tx_gemini_db:
                gemini = config_utils.get_program("gemini", data["config"])
                if "program_versions" in data["config"].get("resources", {}):
                    gemini_ver = programs.get_version("gemini", config=data["config"])
                else:
                    gemini_ver = None
                # Recent versions of gemini allow loading only passing variants
                load_opts = ""
                if not gemini_ver or LooseVersion(gemini_ver) > LooseVersion("0.6.2.1"):
                    load_opts += " --passonly"
                # For small test files, skip gene table loading which takes a long time
                if gemini_ver and LooseVersion(gemini_ver) > LooseVersion("0.6.4"):
                    if _is_small_vcf(gemini_vcf):
                        load_opts += " --skip-gene-tables"
                    if "/test_automated_output/" in gemini_vcf:
                        load_opts += " --test-mode"
                num_cores = data["config"]["algorithm"].get("num_cores", 1)
                cmd = "{gemini} load {load_opts} -v {gemini_vcf} -t snpEff --cores {num_cores} {tx_gemini_db}"
                cmd = cmd.format(**locals())
                do.run(cmd, "Create gemini database for %s %s" % (name, caller), data)
    return [[(name, caller), {"db": gemini_db if utils.file_exists(gemini_db) else None,
                              "vcf": gemini_vcf if is_batch else None}]]

def _is_small_vcf(vcf_file):
    """Check for small VCFs which we want to analyze quicker.
    """
    count = 0
    small_thresh = 250
    with utils.open_gzipsafe(vcf_file) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                count += 1
            if count > small_thresh:
                return False
    return True

def get_multisample_vcf(fnames, name, caller, data):
    """Retrieve a multiple sample VCF file in a standard location.

    Handles inputs with multiple repeated input files from batches.
    """
    unique_fnames = []
    for f in fnames:
        if f not in unique_fnames:
            unique_fnames.append(f)
    out_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "gemini"))
    if len(unique_fnames) > 1:
        gemini_vcf = os.path.join(out_dir, "%s-%s.vcf.gz" % (name, caller))
        return vcfutils.merge_variant_files(unique_fnames, gemini_vcf, data["sam_ref"],
                                            data["config"])
    else:
        gemini_vcf = os.path.join(out_dir, "%s-%s%s" % (name, caller, utils.splitext_plus(unique_fnames[0])[1]))
        utils.symlink_plus(unique_fnames[0], gemini_vcf)
        return gemini_vcf

def _has_gemini(config):
    try:
        gemini = config_utils.get_program("gemini", config)
    except config_utils.CmdNotFound:
        return False
    try:
        p = subprocess.Popen([gemini, "-h"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        p.wait()
        p.stdout.close()
        if p.returncode not in [0, 1]:
            return False
    except OSError:
        return False
    return True

def do_db_build(samples, check_gemini=True, need_bam=True):
    """Confirm we should build a gemini database: need gemini + human samples.
    """
    genomes = set()
    for data in samples:
        if not need_bam or data.get("work_bam"):
            genomes.add(data["genome_build"])
    if len(genomes) == 1:
        return (samples[0]["genome_resources"].get("aliases", {}).get("human", False)
                and (not check_gemini or _has_gemini(samples[0]["config"])))
    else:
        return False

def get_gemini_files(data):
    """Enumerate available gemini data files in a standard installation.
    """
    try:
        from gemini import annotations, config
    except ImportError:
        return {}
    return {"base": config.read_gemini_config()["annotation_dir"],
            "files": annotations.get_anno_files().values()}

def _group_by_batches(samples, check_fn):
    """Group data items into batches, providing details to retrieve results.
    """
    batch_groups = collections.defaultdict(list)
    singles = []
    out_retrieve = []
    extras = []
    for data in [x[0] for x in samples]:
        if check_fn(data):
            batch = data.get("metadata", {}).get("batch")
            name = str(data["name"][-1])
            if batch:
                out_retrieve.append((str(batch), data))
            else:
                out_retrieve.append((name, data))
            for vrn in data["variants"]:
                if batch:
                    batch_groups[(str(batch), vrn["variantcaller"])].append((vrn["vrn_file"], data))
                else:
                    singles.append((name, vrn["variantcaller"], data, vrn["vrn_file"]))
        else:
            extras.append(data)
    return batch_groups, singles, out_retrieve, extras

def _has_variant_calls(data):
    return data["work_bam"] and data.get("vrn_file") and vcfutils.vcf_has_variants(data["vrn_file"])

def prep_db_parallel(samples, parallel_fn):
    """Prepares gemini databases in parallel, handling jointly called populations.
    """
    batch_groups, singles, out_retrieve, extras = _group_by_batches(samples, _has_variant_calls)
    to_process = []
    has_batches = False
    for (name, caller), info in batch_groups.iteritems():
        fnames = [x[0] for x in info]
        to_process.append([fnames, (str(name), caller, True), [x[1] for x in info], info[0][1]])
        has_batches = True
    for name, caller, data, fname in singles:
        to_process.append([[fname], (str(name), caller, False), [data], data])
    if len(samples) > 0 and not do_db_build([x[0] for x in samples], check_gemini=False) and not has_batches:
        return samples
    output = parallel_fn("prep_gemini_db", to_process)
    out_fetch = {}
    for batch_id, out_file in output:
        out_fetch[tuple(batch_id)] = out_file
    out = []
    for batch_name, data in out_retrieve:
        out_variants = []
        for vrn in data["variants"]:
            vrn["population"] = out_fetch[(batch_name, vrn["variantcaller"])]
            out_variants.append(vrn)
        data["variants"] = out_variants
        out.append([data])
    for x in extras:
        out.append([x])
    return out
