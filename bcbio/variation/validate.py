"""Perform validation of final calls against known reference materials.

Automates the process of checking pipeline results against known valid calls
to identify discordant variants. This provides a baseline for ensuring the
validity of pipeline updates and algorithm changes.
"""
import collections
import contextlib
import csv
import os
import shutil
import subprocess
import time

from pysam import VariantFile
import toolz as tz
import yaml

from bcbio import broad, utils
from bcbio.bam import callable
from bcbio.distributed.transaction import file_transaction
from bcbio.heterogeneity import bubbletree
from bcbio.pipeline import config_utils, shared
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import bedutils, validateplot, vcfutils, multi, naming

# ## Individual sample comparisons

def _get_validate(data):
    """Retrieve items to validate, from single samples or from combined joint calls.
    """
    if data.get("vrn_file") and tz.get_in(["config", "algorithm", "validate"], data):
        return data
    elif "group_orig" in data:
        for sub in multi.get_orig_items(data):
            if "validate" in sub["config"]["algorithm"]:
                sub_val = utils.deepish_copy(sub)
                sub_val["vrn_file"] = data["vrn_file"]
                return sub_val
    return None

def normalize_input_path(x, data):
    """Normalize path for input files, handling relative paths.
    Looks for non-absolute paths in local and fastq directories
    """
    if x is None:
        return None
    elif os.path.isabs(x):
        return os.path.normpath(x)
    else:
        for d in [data["dirs"].get("fastq"), data["dirs"].get("work")]:
            if d:
                cur_x = os.path.normpath(os.path.join(d, x))
                if os.path.exists(cur_x):
                    return cur_x
        raise IOError("Could not find validation file %s" % x)

def _gunzip(f, data):
    if f is None:
        return None
    elif f.endswith(".gz"):
        out_file = f.replace(".gz", "")
        if not utils.file_exists(out_file):
            with file_transaction(data, out_file) as tx_out_file:
                cmd = "gunzip -c {f} > {tx_out_file}"
                do.run(cmd.format(**locals()), "gunzip input file")
        return out_file
    else:
        return f

def _get_caller(data):
    callers = [tz.get_in(["config", "algorithm", "jointcaller"], data),
               tz.get_in(["config", "algorithm", "variantcaller"], data),
               "precalled"]
    return [c for c in callers if c][0]

def _get_caller_supplement(caller, data):
    """Some callers like MuTect incorporate a second caller for indels.
    """
    if caller == "mutect":
        icaller = tz.get_in(["config", "algorithm", "indelcaller"], data)
        if icaller:
            caller = "%s/%s" % (caller, icaller)
    return caller

def _normalize_cwl_inputs(items):
    """Extract variation and validation data from CWL input list of batched samples.
    """
    with_validate = []
    vrn_files = []
    for data in items:
        if tz.get_in(["config", "algorithm", "validate"], data):
            with_validate.append(data)
        if data.get("vrn_file"):
            vrn_files.append(data["vrn_file"])
    if len(with_validate) == 0:
        return items[0]
    else:
        assert len(set([tz.get_in(["config", "algorithm", "validate"], data) for data in with_validate])) == 1
        assert len(set(vrn_files)) == 1
        data = with_validate[0]
        data["vrn_file"] = vrn_files[0]
        return data

def compare_to_rm(data):
    """Compare final variant calls against reference materials of known calls.
    """
    if isinstance(data, (list, tuple)):
        data = _normalize_cwl_inputs(data)
    toval_data = _get_validate(data)
    if toval_data:
        caller = _get_caller(toval_data)
        sample = dd.get_sample_name(toval_data)
        base_dir = utils.safe_makedir(os.path.join(toval_data["dirs"]["work"], "validate", sample, caller))

        if isinstance(toval_data["vrn_file"], (list, tuple)):
            raise NotImplementedError("Multiple input files for validation: %s" % toval_data["vrn_file"])
        else:
            vrn_file = os.path.abspath(toval_data["vrn_file"])
        rm_file = normalize_input_path(toval_data["config"]["algorithm"]["validate"], toval_data)
        rm_interval_file = _gunzip(normalize_input_path(toval_data["config"]["algorithm"].get("validate_regions"),
                                                        toval_data),
                                   toval_data)
        rm_interval_file = bedutils.clean_file(rm_interval_file, toval_data,
                                               bedprep_dir=utils.safe_makedir(os.path.join(base_dir, "bedprep")))
        rm_file = naming.handle_synonyms(rm_file, dd.get_ref_file(data), data["genome_build"], base_dir, data)
        rm_interval_file = (naming.handle_synonyms(rm_interval_file, dd.get_ref_file(data),
                                                   data["genome_build"], base_dir, data)
                            if rm_interval_file else None)
        vmethod = tz.get_in(["config", "algorithm", "validate_method"], data, "rtg")
        if not vcfutils.vcf_has_variants(vrn_file):
            # RTG can fail on totally empty files. Skip these since we have nothing.
            pass
        # empty validation file, every call is a false positive
        elif not vcfutils.vcf_has_variants(rm_file):
            eval_files = _setup_call_fps(vrn_file, rm_interval_file, base_dir, toval_data)
            data["validate"] = _rtg_add_summary_file(eval_files, base_dir, toval_data)
        elif vmethod == "rtg":
            eval_files = _run_rtg_eval(vrn_file, rm_file, rm_interval_file, base_dir, toval_data)
            data["validate"] = _rtg_add_summary_file(eval_files, base_dir, toval_data)
        elif vmethod == "hap.py":
            data["validate"] = _run_happy_eval(vrn_file, rm_file, rm_interval_file, base_dir, toval_data)
        elif vmethod == "bcbio.variation":
            data["validate"] = _run_bcbio_variation(vrn_file, rm_file, rm_interval_file, base_dir,
                                                    sample, caller, toval_data)
    return [[data]]

# ## Empty truth sets

def _setup_call_fps(vrn_file, rm_bed, base_dir, data):
    """Create set of false positives for inputs with empty truth sets.
    """
    out_file = os.path.join(base_dir, "fp.vcf.gz")
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = ("bcftools view -R {rm_bed} -f 'PASS,.' {vrn_file} -O z -o {tx_out_file}")
            do.run(cmd.format(**locals()), "Prepare false positives with empty reference", data)
    return {"fp": out_file}

# ## Real Time Genomics vcfeval

def _get_sample_and_caller(data):
    return [tz.get_in(["metadata", "validate_sample"], data) or dd.get_sample_name(data),
            _get_caller_supplement(_get_caller(data), data)]

def _rtg_add_summary_file(eval_files, base_dir, data):
    """Parse output TP FP and FN files to generate metrics for plotting.
    """
    out_file = os.path.join(base_dir, "validate-summary.csv")
    if not utils.file_uptodate(out_file, eval_files.get("tp", eval_files["fp"])):
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                writer = csv.writer(out_handle)
                writer.writerow(["sample", "caller", "vtype", "metric", "value"])
                base = _get_sample_and_caller(data)
                for metric in ["tp", "fp", "fn"]:
                    for vtype, bcftools_types in [("SNPs", "--types snps"),
                                                  ("Indels", "--exclude-types snps")]:
                        in_file = eval_files.get(metric)
                        if in_file and os.path.exists(in_file):
                            cmd = ("bcftools view {bcftools_types} {in_file} | grep -v ^# | wc -l")
                            count = int(subprocess.check_output(cmd.format(**locals()), shell=True))
                        else:
                            count = 0
                        writer.writerow(base + [vtype, metric, count])
    eval_files["summary"] = out_file
    return eval_files

def _prepare_inputs(vrn_file, rm_file, rm_interval_file, base_dir, data):
    """Prepare input VCF and BED files for validation.
    """
    if not rm_file.endswith(".vcf.gz") or not os.path.exists(rm_file + ".tbi"):
        rm_file = vcfutils.bgzip_and_index(rm_file, data["config"], out_dir=base_dir)
    if len(vcfutils.get_samples(vrn_file)) > 1:
        base, ext = utils.splitext_plus(os.path.basename(vrn_file))
        sample_file = os.path.join(base_dir, "%s-%s%s" % (base, dd.get_sample_name(data), ext))
        vrn_file = vcfutils.select_sample(vrn_file, dd.get_sample_name(data), sample_file, data["config"])
    if not vrn_file.endswith(".vcf.gz") or not os.path.exists(vrn_file + ".tbi"):
        vrn_file = vcfutils.bgzip_and_index(vrn_file, data["config"], out_dir=base_dir)

    interval_bed = _get_merged_intervals(rm_interval_file, base_dir, data)
    return vrn_file, rm_file, interval_bed

def _run_rtg_eval(vrn_file, rm_file, rm_interval_file, base_dir, data):
    """Run evaluation of a caller against the truth set using rtg vcfeval.
    """
    out_dir = os.path.join(base_dir, "rtg")
    if not utils.file_exists(os.path.join(out_dir, "done")):
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        vrn_file, rm_file, interval_bed = _prepare_inputs(vrn_file, rm_file, rm_interval_file, base_dir, data)

        rtg_ref = tz.get_in(["reference", "rtg"], data)
        assert rtg_ref and os.path.exists(rtg_ref), ("Did not find rtg indexed reference file for validation:\n%s\n"
                                                     "Run bcbio_nextgen.py upgrade --data --aligners rtg" % rtg_ref)
        # handle CWL where we have a reference to a single file in the RTG directory
        if os.path.isfile(rtg_ref):
            rtg_ref = os.path.dirname(rtg_ref)

        # get core and memory usage from standard configuration
        threads = min(dd.get_num_cores(data), 6)
        resources = config_utils.get_resources("rtg", data["config"])
        memory = config_utils.adjust_opts(resources.get("jvm_opts", ["-Xms500m", "-Xmx1500m"]),
                                          {"algorithm": {"memory_adjust": {"magnitude": threads,
                                                                           "direction": "increase"}}})
        jvm_stack = [x for x in memory if x.startswith("-Xms")]
        jvm_mem = [x for x in memory if x.startswith("-Xmx")]
        jvm_stack = jvm_stack[0] if len(jvm_stack) > 0 else "-Xms500m"
        jvm_mem = jvm_mem[0].replace("-Xmx", "") if len(jvm_mem) > 0 else "3g"
        cmd = ["rtg", "vcfeval", "--threads", str(threads),
               "-b", rm_file, "--bed-regions", interval_bed,
               "-c", vrn_file, "-t", rtg_ref, "-o", out_dir]
        cmd += ["--vcf-score-field='%s'" % (_pick_best_quality_score(vrn_file))]
        mem_export = "export RTG_JAVA_OPTS='%s' && export RTG_MEM=%s" % (jvm_stack, jvm_mem)
        cmd = mem_export + " && " + " ".join(cmd)
        do.run(cmd, "Validate calls using rtg vcfeval", data)
    out = {"fp": os.path.join(out_dir, "fp.vcf.gz"),
           "fn": os.path.join(out_dir, "fn.vcf.gz")}
    tp_calls = os.path.join(out_dir, "tp.vcf.gz")
    tp_baseline = os.path.join(out_dir, "tp-baseline.vcf.gz")
    if os.path.exists(tp_baseline):
        out["tp"] = tp_baseline
        out["tp-calls"] = tp_calls
    else:
        out["tp"] = tp_calls
    return out

def _pick_best_quality_score(vrn_file):
    """Flexible quality score selection, picking the best available.

    Implementation based on discussion:

    https://github.com/chapmanb/bcbio-nextgen/commit/a538cecd86c0000d17d3f9d4f8ac9d2da04f9884#commitcomment-14539249

    (RTG=AVR/GATK=VQSLOD/MuTect=t_lod_fstar, otherwise GQ, otherwise QUAL, otherwise DP.)

    For MuTect, it's not clear how to get t_lod_fstar, the right quality score, into VCF cleanly.
    MuTect2 has TLOD in the INFO field.
    """
    # pysam fails on checking reference contigs if input is empty
    if not vcfutils.vcf_has_variants(vrn_file):
        return "DP"
    to_check = 25
    scores = collections.defaultdict(int)
    try:
        in_handle = VariantFile(vrn_file)
    except ValueError:
        raise ValueError("Failed to parse input file in preparation for validation: %s" % vrn_file)
    with contextlib.closing(in_handle) as val_in:
        for i, rec in enumerate(val_in):
            if i > to_check:
                break
            if rec.info.get("VQSLOD") is not None:
                scores["INFO=VQSLOD"] += 1
            if rec.info.get("TLOD") is not None:
                scores["INFO=TLOD"] += 1
            for skey in ["AVR", "GQ", "DP"]:
                if len(rec.samples) > 0 and rec.samples[0].get(skey) is not None:
                    scores[skey] += 1
            if rec.qual:
                scores["QUAL"] += 1
    for key in ["AVR", "INFO=VQSLOD", "INFO=TLOD", "GQ", "QUAL", "DP"]:
        if scores[key] > 0:
            return key
    raise ValueError("Did not find quality score for validation from %s" % vrn_file)

def _get_merged_intervals(rm_interval_file, base_dir, data):
    """Retrieve intervals to run validation on, merging reference and callable BED files.
    """
    a_intervals = get_analysis_intervals(data)
    if a_intervals:
        final_intervals = shared.remove_lcr_regions(a_intervals, [data])
        if rm_interval_file:
            caller = _get_caller(data)
            sample = dd.get_sample_name(data)
            combo_intervals = os.path.join(base_dir, "%s-%s-%s-wrm.bed" %
                                           (utils.splitext_plus(os.path.basename(final_intervals))[0],
                                            sample, caller))
            if not utils.file_uptodate(combo_intervals, final_intervals):
                with file_transaction(data, combo_intervals) as tx_out_file:
                    with utils.chdir(os.path.dirname(tx_out_file)):
                        # Copy files locally to avoid issues on shared filesystems
                        # where BEDtools has trouble accessing the same base
                        # files from multiple locations
                        a = os.path.basename(final_intervals)
                        b = os.path.basename(rm_interval_file)
                        try:
                            shutil.copyfile(final_intervals, a)
                        except IOError:
                            time.sleep(60)
                            shutil.copyfile(final_intervals, a)
                        try:
                            shutil.copyfile(rm_interval_file, b)
                        except IOError:
                            time.sleep(60)
                            shutil.copyfile(rm_interval_file, b)
                        cmd = ("bedtools intersect -nonamecheck -a {a} -b {b} > {tx_out_file}")
                        do.run(cmd.format(**locals()), "Intersect callable intervals for rtg vcfeval")
            final_intervals = combo_intervals
    else:
        assert rm_interval_file, "No intervals to subset analysis with"
        final_intervals = shared.remove_lcr_regions(rm_interval_file, [data])
    return final_intervals

def get_analysis_intervals(data):
    """Retrieve analysis regions for the current variant calling pipeline.
    """
    if data.get("ensemble_bed"):
        return data["ensemble_bed"]
    elif dd.get_callable_regions(data):
        return dd.get_callable_regions(data)
    elif data.get("align_bam"):
        return callable.sample_callable_bed(data["align_bam"], dd.get_ref_file(data), data)
    elif data.get("work_bam"):
        return callable.sample_callable_bed(data["work_bam"], dd.get_ref_file(data), data)
    elif data.get("work_bam_callable"):
        return callable.sample_callable_bed(data["work_bam_callable"], dd.get_ref_file(data), data)
    elif tz.get_in(["config", "algorithm", "callable_regions"], data):
        return tz.get_in(["config", "algorithm", "callable_regions"], data)
    elif tz.get_in(["config", "algorithm", "variant_regions"], data):
        return tz.get_in(["config", "algorithm", "variant_regions"], data)

# ## hap.py

def _run_happy_eval(vrn_file, rm_file, rm_interval_file, base_dir, data):
    """Validation with hap.py: https://github.com/Illumina/hap.py

    XXX Does not yet parse out metrics for plotting.
    """
    out_dir = utils.safe_makedir(os.path.join(base_dir, "happy"))
    out_prefix = os.path.join(out_dir, "val")
    if not utils.file_exists(out_prefix + ".summary.csv"):
        vrn_file, rm_file, interval_bed = _prepare_inputs(vrn_file, rm_file, rm_interval_file, base_dir, data)
        cmd = ["hap.py", "-V", "-f", interval_bed, "-r", dd.get_ref_file(data),
               "-l", ",".join(_get_location_list(interval_bed)),
               "-o", out_prefix, rm_file, vrn_file]
        do.run(cmd, "Validate calls using hap.py", data)
    return {"vcf": out_prefix + ".vcf.gz"}

def _get_location_list(interval_bed):
    """Retrieve list of locations to analyze from input BED file.
    """
    import pybedtools
    regions = collections.OrderedDict()
    for region in pybedtools.BedTool(interval_bed):
        regions[str(region.chrom)] = None
    return regions.keys()

# ## bcbio.variation comparison -- deprecated approach

def _run_bcbio_variation(vrn_file, rm_file, rm_interval_file, base_dir, sample, caller, data):
    """Run validation of a caller against the truth set using bcbio.variation.
    """
    val_config_file = _create_validate_config_file(vrn_file, rm_file, rm_interval_file,
                                                   base_dir, data)
    work_dir = os.path.join(base_dir, "work")
    out = {"summary": os.path.join(work_dir, "validate-summary.csv"),
           "grading": os.path.join(work_dir, "validate-grading.yaml"),
           "discordant": os.path.join(work_dir, "%s-eval-ref-discordance-annotate.vcf" % sample)}
    if not utils.file_exists(out["discordant"]) or not utils.file_exists(out["grading"]):
        bcbio_variation_comparison(val_config_file, base_dir, data)
    out["concordant"] = filter(os.path.exists,
                                [os.path.join(work_dir, "%s-%s-concordance.vcf" % (sample, x))
                                 for x in ["eval-ref", "ref-eval"]])[0]
    return out

def bcbio_variation_comparison(config_file, base_dir, data):
    """Run a variant comparison using the bcbio.variation toolkit, given an input configuration.
    """
    tmp_dir = utils.safe_makedir(os.path.join(base_dir, "tmp"))
    resources = config_utils.get_resources("bcbio_variation", data["config"])
    jvm_opts = resources.get("jvm_opts", ["-Xms750m", "-Xmx2g"])
    cmd = ["bcbio-variation"] + jvm_opts + broad.get_default_jvm_opts(tmp_dir) + \
          ["variant-compare", config_file]
    do.run(cmd, "Comparing variant calls using bcbio.variation", data)

def _create_validate_config_file(vrn_file, rm_file, rm_interval_file,
                                 base_dir, data):
    config_dir = utils.safe_makedir(os.path.join(base_dir, "config"))
    config_file = os.path.join(config_dir, "validate.yaml")
    if not utils.file_uptodate(config_file, vrn_file):
        with file_transaction(data, config_file) as tx_config_file:
            with open(tx_config_file, "w") as out_handle:
                out = _create_validate_config(vrn_file, rm_file, rm_interval_file,
                                              base_dir, data)
                yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return config_file

def _create_validate_config(vrn_file, rm_file, rm_interval_file, base_dir, data):
    """Create a bcbio.variation configuration input for validation.
    """
    ref_call = {"file": str(rm_file), "name": "ref", "type": "grading-ref",
                "fix-sample-header": True, "remove-refcalls": True}
    a_intervals = get_analysis_intervals(data)
    if a_intervals:
        a_intervals = shared.remove_lcr_regions(a_intervals, [data])
    if rm_interval_file:
        ref_call["intervals"] = rm_interval_file
    eval_call = {"file": vrn_file, "name": "eval", "remove-refcalls": True}
    exp = {"sample": data["name"][-1],
           "ref": dd.get_ref_file(data),
           "approach": "grade",
           "calls": [ref_call, eval_call]}
    if a_intervals:
        exp["intervals"] = os.path.abspath(a_intervals)
    if data.get("align_bam"):
        exp["align"] = data["align_bam"]
    elif data.get("work_bam"):
        exp["align"] = data["work_bam"]
    return {"dir": {"base": base_dir, "out": "work", "prep": "work/prep"},
            "experiments": [exp]}

# ## Summarize comparisons

def _flatten_grading(stats):
    vtypes = ["snp", "indel"]
    cat = "concordant"
    for vtype in vtypes:
        yield vtype, cat, stats[cat][cat].get(vtype, 0)
    for vtype in vtypes:
        for vclass, vitems in sorted(stats["discordant"].get(vtype, {}).iteritems()):
            for vreason, val in sorted(vitems.iteritems()):
                yield vtype, "discordant-%s-%s" % (vclass, vreason), val
            yield vtype, "discordant-%s-total" % vclass, sum(vitems.itervalues())

def _has_grading_info(samples):
    for data in (x[0] for x in samples):
        for variant in data.get("variants", []):
            if variant.get("validate"):
                return True
    return False

def _group_validate_samples(samples):
    extras = []
    validated = collections.defaultdict(list)
    for data in (x[0] for x in samples):
        is_v = False
        for variant in data.get("variants", []):
            if variant.get("validate"):
                is_v = True
        if is_v:
            for batch_key in (["metadata", "validate_batch"], ["metadata", "batch"],
                              ["description"]):
                vname = tz.get_in(batch_key, data)
                if vname:
                    break
            if isinstance(vname, (list, tuple)):
                vname = vname[0]
            validated[vname].append(data)
        else:
            extras.append([data])
    return validated, extras

def summarize_grading(samples):
    """Provide summaries of grading results across all samples.
    """
    if not _has_grading_info(samples):
        return samples
    validate_dir = utils.safe_makedir(os.path.join(samples[0][0]["dirs"]["work"], "validate"))
    header = ["sample", "caller", "variant.type", "category", "value"]
    validated, out = _group_validate_samples(samples)
    for vname, vitems in validated.iteritems():
        out_csv = os.path.join(validate_dir, "grading-summary-%s.csv" % vname)
        with open(out_csv, "w") as out_handle:
            writer = csv.writer(out_handle)
            writer.writerow(header)
            plot_data = []
            plot_files = []
            for data in sorted(vitems, key=lambda x: x.get("lane", dd.get_sample_name(x))):
                for variant in data.get("variants", []):
                    if variant.get("validate"):
                        variant["validate"]["grading_summary"] = out_csv
                        if tz.get_in(["validate", "grading"], variant):
                            for row in _get_validate_plotdata_yaml(variant, data):
                                writer.writerow(row)
                                plot_data.append(row)
                        elif tz.get_in(["validate", "summary"], variant):
                            plot_files.append(variant["validate"]["summary"])
        if plot_files:
            plots = validateplot.classifyplot_from_plotfiles(plot_files, out_csv)
        elif plot_data:
            plots = validateplot.create(plot_data, header, 0, data["config"],
                                        os.path.splitext(out_csv)[0])
        else:
            plots = None
        for data in vitems:
            for variant in data.get("variants", []):
                if variant.get("validate"):
                    variant["validate"]["grading_plots"] = plots
            out.append([data])
    return out

def _get_validate_plotdata_yaml(variant, data):
    """Retrieve validation plot data from grading YAML file (old style).
    """
    with open(variant["validate"]["grading"]) as in_handle:
        grade_stats = yaml.load(in_handle)
    for sample_stats in grade_stats:
        sample = sample_stats["sample"]
        for vtype, cat, val in _flatten_grading(sample_stats):
            yield [sample, variant.get("variantcaller", ""),
                   vtype, cat, val]

# ## Summarize by frequency

def freq_summary(val_file, call_file, truth_file, target_name):
    """Summarize true and false positive calls by variant type and frequency.

    Resolve differences in true/false calls based on output from hap.py:
    https://github.com/sequencing/hap.py
    """
    out_file = "%s-freqs.csv" % utils.splitext_plus(val_file)[0]
    truth_freqs = _read_truth_freqs(truth_file)
    call_freqs = _read_call_freqs(call_file, target_name)
    with VariantFile(val_file) as val_in:
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle)
            writer.writerow(["vtype", "valclass", "freq"])
            for rec in val_in:
                call_type = _classify_rec(rec)
                val_type = _get_validation_status(rec)
                key = _get_key(rec)
                freq = truth_freqs.get(key, call_freqs.get(key, 0.0))
                writer.writerow([call_type, val_type, freq])
    return out_file

def _get_key(rec):
    return (rec.contig, rec.pos, rec.ref, rec.alts[0])

def _classify_rec(rec):
    """Determine class of variant in the record.
    """
    if max([len(x) for x in rec.alleles]) == 1:
        return "snp"
    else:
        return "indel"

def _get_validation_status(rec):
    """Retrieve the status of the validation, supporting hap.py output
    """
    return rec.info["type"]

def _read_call_freqs(in_file, sample_name):
    """Identify frequencies for calls in the input file.
    """
    out = {}
    with VariantFile(in_file) as call_in:
        for rec in call_in:
            if rec.filter.keys() == ["PASS"]:
                for name, sample in rec.samples.items():
                    if name == sample_name:
                        alt, depth = bubbletree.sample_alt_and_depth(sample)
                        if depth > 0:
                            out[_get_key(rec)] = float(alt) / float(depth)
    return out

def _read_truth_freqs(in_file):
    """Read frequency of calls from truth VCF.

    Currently handles DREAM data, needs generalization for other datasets.
    """
    out = {}
    with VariantFile(in_file) as bcf_in:
        for rec in bcf_in:
            freq = float(rec.info.get("VAF", 1.0))
            out[_get_key(rec)] = freq
    return out
