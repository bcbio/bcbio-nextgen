"""Provide infrastructure to allow exploration of variations within populations.

Uses the gemini framework (https://github.com/arq5x/gemini) to build SQLite
database of variations for query and evaluation.
"""
import collections
import csv
from distutils.version import LooseVersion
import os
import subprocess

import toolz as tz

from bcbio import install, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do, programs
from bcbio.variation import multiallelic, vcfutils

def prep_gemini_db(fnames, call_info, samples, extras):
    """Prepare a gemini database from VCF inputs prepared with snpEff.
    """
    data = samples[0]
    out_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "gemini"))
    name, caller, is_batch = call_info
    gemini_db = os.path.join(out_dir, "%s-%s.db" % (name, caller))
    multisample_vcf = get_multisample_vcf(fnames, name, caller, data)
    gemini_vcf = multiallelic.to_single(multisample_vcf, data)
    use_gemini_quick = (do_db_build(samples) and
                        any(vcfutils.vcf_has_variants(f) for f in fnames))
    if not utils.file_exists(gemini_db) and use_gemini_quick:
        use_gemini = do_db_build(samples) and any(vcfutils.vcf_has_variants(f) for f in fnames)
        if use_gemini:
            ped_file = create_ped_file(samples + extras, gemini_vcf)
            gemini_db = create_gemini_db(gemini_vcf, data, gemini_db, ped_file)
    return [[(name, caller), {"db": gemini_db if utils.file_exists(gemini_db) else None,
                              "vcf": multisample_vcf if is_batch else None}]]

def create_gemini_db(gemini_vcf, data, gemini_db=None, ped_file=None):
    if not gemini_db:
        gemini_db = "%s.db" % utils.splitext_plus(gemini_vcf)[0]
    if not utils.file_exists(gemini_db):
        if not vcfutils.vcf_has_variants(gemini_vcf):
            return None
        with file_transaction(data, gemini_db) as tx_gemini_db:
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
            # Skip CADD or gerp-bp if neither are loaded
            if gemini_ver and LooseVersion(gemini_ver) >= LooseVersion("0.7.0"):
                gemini_dir = install.get_gemini_dir(data)
                for skip_cmd, check_file in [("--skip-cadd", "whole_genome_SNVs.tsv.compressed.gz")]:
                    if not os.path.exists(os.path.join(gemini_dir, check_file)):
                        load_opts += " %s" % skip_cmd
            # skip gerp-bp which slows down loading
            load_opts += " --skip-gerp-bp "
            num_cores = data["config"]["algorithm"].get("num_cores", 1)
            tmpdir = os.path.dirname(tx_gemini_db)
            eanns = _get_effects_flag(data)
            # Apply custom resource specifications, allowing use of alternative annotation_dir
            resources = config_utils.get_resources("gemini", data["config"])
            gemini_opts = " ".join([str(x) for x in resources["options"]]) if resources.get("options") else ""
            cmd = ("{gemini} {gemini_opts} load {load_opts} -v {gemini_vcf} {eanns} --cores {num_cores} "
                   "--tempdir {tmpdir} {tx_gemini_db}")
            cmd = cmd.format(**locals())
            do.run(cmd, "Create gemini database for %s" % gemini_vcf, data)
            if ped_file:
                cmd = [gemini, "amend", "--sample", ped_file, tx_gemini_db]
                do.run(cmd, "Add PED file to gemini database", data)
    return gemini_db

def _get_effects_flag(data):
    effects_config = tz.get_in(("config", "algorithm", "effects"), data, "snpeff")
    if effects_config == "snpeff":
        return "-t snpEff"
    elif effects_config == "vep":
        return "-t VEP"
    else:
        return ""

def get_affected_status(data):
    """Retrieve the affected/unaffected status of sample.

    Uses unaffected (1), affected (2), unknown (0) coding from PED files:

    http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
    """
    affected = set(["tumor", "affected"])
    unaffected = set(["normal", "unaffected"])
    phenotype = str(tz.get_in(["metadata", "phenotype"], data, "")).lower()
    if phenotype in affected:
        return 2
    elif phenotype in unaffected:
        return 1
    else:
        return 0

def get_gender(data):
    """Retrieve gender from metadata, codified as male/female/unknown.
    """
    g = dd.get_gender(data)
    if g and str(g).lower() in ["male", "m"]:
        return "male"
    elif g and str(g).lower() in ["female", "f"]:
        return "female"
    else:
        return "unknown"

def create_ped_file(samples, base_vcf):
    """Create a GEMINI-compatible PED file, including gender, family and phenotype information.

    Checks for a specified `ped` file in metadata, and will use sample information from this file
    before reconstituting from metadata information.
    """
    out_file = "%s.ped" % utils.splitext_plus(base_vcf)[0]
    sample_ped_lines = {}
    header = ["#Family_ID", "Individual_ID", "Paternal_ID", "Maternal_ID", "Sex", "Phenotype", "Ethnicity"]
    for md_ped in list(set([x for x in [tz.get_in(["metadata", "ped"], data)
                                        for data in samples] if x is not None])):
        with open(md_ped) as in_handle:
            reader = csv.reader(in_handle, dialect="excel-tab")
            for parts in reader:
                if parts[0].startswith("#") and len(parts) > len(header):
                    header = header + parts[len(header):]
                else:
                    sample_ped_lines[parts[1]] = parts
    if not utils.file_exists(out_file):
        with file_transaction(samples[0], out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                writer = csv.writer(out_handle, dialect="excel-tab")
                writer.writerow(header)
                batch = _find_shared_batch(samples)
                for data in samples:
                    gender = {"male": 1, "female": 2, "unknown": 0}.get(get_gender(data))
                    sname = dd.get_sample_name(data)
                    if sname in sample_ped_lines:
                        writer.writerow(sample_ped_lines[sname])
                    else:
                        writer.writerow([batch, sname, "-9", "-9",
                                         gender, get_affected_status(data), "-9"])
    return out_file

def _find_shared_batch(samples):
    for data in samples:
        batch = tz.get_in(["metadata", "batch"], data, dd.get_sample_name(data))
        if not isinstance(batch, (list, tuple)):
            return batch

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
        vrn_file_batch = None
        for variant in data["variants"]:
            if variant["variantcaller"] == caller and variant.get("vrn_file_batch"):
                vrn_file_batch = variant["vrn_file_batch"]
        if vrn_file_batch:
            utils.symlink_plus(vrn_file_batch, gemini_vcf)
            return gemini_vcf
        else:
            return vcfutils.merge_variant_files(unique_fnames, gemini_vcf, data["sam_ref"],
                                                data["config"])
    else:
        gemini_vcf = os.path.join(out_dir, "%s-%s%s" % (name, caller, utils.splitext_plus(unique_fnames[0])[1]))
        utils.symlink_plus(unique_fnames[0], gemini_vcf)
        return gemini_vcf

def _has_gemini(data):
    from bcbio import install
    gemini_dir = install.get_gemini_dir(data)
    return ((os.path.exists(gemini_dir) and len(os.listdir(gemini_dir)) > 0)
            and os.path.exists(os.path.join(os.path.dirname(gemini_dir), "gemini-config.yaml")))

def do_db_build(samples, need_bam=True, gresources=None):
    """Confirm we should build a gemini database: need gemini + human samples + not in tool_skip.
    """
    genomes = set()
    for data in samples:
        if not need_bam or data.get("align_bam"):
            genomes.add(data["genome_build"])
        if "gemini" in utils.get_in(data, ("config", "algorithm", "tools_off"), []):
            return False
    if len(genomes) == 1:
        if not gresources:
            gresources = samples[0]["genome_resources"]
        return (tz.get_in(["aliases", "human"], gresources, False)
                and _has_gemini(samples[0]))
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
            batch = tz.get_in(["metadata", "batch"], data)
            name = str(data["name"][-1])
            if batch:
                out_retrieve.append((str(batch), data))
            else:
                out_retrieve.append((name, data))
            for vrn in data["variants"]:
                if vrn.get("population", True):
                    if batch:
                        batch_groups[(str(batch), vrn["variantcaller"])].append((vrn["vrn_file"], data))
                    else:
                        singles.append((name, vrn["variantcaller"], data, vrn["vrn_file"]))
        else:
            extras.append(data)
    return batch_groups, singles, out_retrieve, extras

def _has_variant_calls(data):
    if data.get("align_bam"):
        for vrn in data["variants"]:
            if vrn.get("vrn_file") and vcfutils.vcf_has_variants(vrn["vrn_file"]):
                return True
    return False

def prep_db_parallel(samples, parallel_fn):
    """Prepares gemini databases in parallel, handling jointly called populations.
    """
    batch_groups, singles, out_retrieve, extras = _group_by_batches(samples, _has_variant_calls)
    to_process = []
    has_batches = False
    for (name, caller), info in batch_groups.iteritems():
        fnames = [x[0] for x in info]
        to_process.append([fnames, (str(name), caller, True), [x[1] for x in info], extras])
        has_batches = True
    for name, caller, data, fname in singles:
        to_process.append([[fname], (str(name), caller, False), [data], extras])
    if len(samples) > 0 and not do_db_build([x[0] for x in samples]) and not has_batches:
        return samples
    output = parallel_fn("prep_gemini_db", to_process)
    out_fetch = {}
    for batch_id, out_file in output:
        out_fetch[tuple(batch_id)] = out_file
    out = []
    for batch_name, data in out_retrieve:
        out_variants = []
        for vrn in data["variants"]:
            use_population = vrn.pop("population", True)
            if use_population:
                vrn["population"] = out_fetch[(batch_name, vrn["variantcaller"])]
            out_variants.append(vrn)
        data["variants"] = out_variants
        out.append([data])
    for x in extras:
        out.append([x])
    return out
