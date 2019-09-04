"""High level summaries of samples and programs with MultiQC.

https://github.com/ewels/MultiQC
"""
import collections
import glob
import io
import json
import mimetypes
import os
import pandas as pd
import shutil
import numpy as np
from collections import OrderedDict

import pybedtools
import six
import toolz as tz
import yaml

from bcbio import utils
from bcbio.cwl import cwlutils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.log import logger
from bcbio.provenance import do, programs
from bcbio.provenance import data as provenancedata
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.bam import ref
from bcbio.qc.qsignature import get_qsig_multiqc_files
from bcbio.structural import annotate
from bcbio.utils import walk_json
from bcbio.variation import bedutils
from bcbio.qc.variant import get_active_vcinfo
from bcbio.upload import get_all_upload_paths_from_sample
from bcbio.variation import coverage

def summary(*samples):
    """Summarize all quality metrics together"""
    samples = list(utils.flatten(samples))
    work_dir = dd.get_work_dir(samples[0])
    multiqc = config_utils.get_program("multiqc", samples[0]["config"])
    if not multiqc:
        logger.debug("multiqc not found. Update bcbio_nextgen.py tools to fix this issue.")
    out_dir = utils.safe_makedir(os.path.join(work_dir, "qc", "multiqc"))
    out_data = os.path.join(out_dir, "multiqc_data")
    out_file = os.path.join(out_dir, "multiqc_report.html")
    file_list = os.path.join(out_dir, "list_files.txt")
    work_samples = cwlutils.unpack_tarballs([utils.deepish_copy(x) for x in samples], samples[0])
    work_samples = _summarize_inputs(work_samples, out_dir)
    if not utils.file_exists(out_file):
        with tx_tmpdir(samples[0], work_dir) as tx_out:
            in_files = _get_input_files(work_samples, out_dir, tx_out)
            in_files += _merge_metrics(work_samples, out_dir)
            if _one_exists(in_files):
                with utils.chdir(out_dir):
                    _create_config_file(out_dir, work_samples)
                    input_list_file = _create_list_file(in_files, file_list)
                    if dd.get_tmp_dir(samples[0]):
                        export_tmp = "export TMPDIR=%s && " % dd.get_tmp_dir(samples[0])
                    else:
                        export_tmp = ""
                    locale_export = utils.locale_export()
                    path_export = utils.local_path_export()
                    other_opts = config_utils.get_resources("multiqc", samples[0]["config"]).get("options", [])
                    other_opts = " ".join([str(x) for x in other_opts])
                    cmd = ("{path_export}{export_tmp}{locale_export} "
                           "{multiqc} -f -l {input_list_file} {other_opts} -o {tx_out}")
                    do.run(cmd.format(**locals()), "Run multiqc")
                    if utils.file_exists(os.path.join(tx_out, "multiqc_report.html")):
                        shutil.move(os.path.join(tx_out, "multiqc_report.html"), out_file)
                        shutil.move(os.path.join(tx_out, "multiqc_data"), out_data)
    samples = _group_by_sample_and_batch(samples)
    if utils.file_exists(out_file) and samples:
        data_files = set()
        for i, data in enumerate(samples):
            data_files.add(os.path.join(out_dir, "report", "metrics", dd.get_sample_name(data) + "_bcbio.txt"))
        data_files.add(os.path.join(out_dir, "report", "metrics", "target_info.yaml"))
        data_files.add(os.path.join(out_dir, "multiqc_config.yaml"))
        [data_files.add(f) for f in glob.glob(os.path.join(out_dir, "multiqc_data", "*"))]
        data_files = [f for f in data_files if f and utils.file_exists(f)]
        if "summary" not in samples[0]:
            samples[0]["summary"] = {}
        samples[0]["summary"]["multiqc"] = {"base": out_file, "secondary": data_files}

        data_json = os.path.join(out_dir, "multiqc_data", "multiqc_data.json")
        data_json_final = _save_uploaded_data_json(samples, data_json, os.path.join(out_dir, "multiqc_data"))
        if data_json_final:
            samples[0]["summary"]["multiqc"]["secondary"].append(data_json_final)

        # Prepare final file list and inputs for downstream usage
        file_list_final = _save_uploaded_file_list(samples, file_list, out_dir)
        if file_list_final:
            samples[0]["summary"]["multiqc"]["secondary"].append(file_list_final)
            if any([cwlutils.is_cwl_run(d) for d in samples]):
                for indir in ["inputs", "report"]:
                    tarball = os.path.join(out_dir, "multiqc-%s.tar.gz" % (indir))
                    if not utils.file_exists(tarball):
                        with utils.chdir(out_dir):
                            cmd = ["tar", "-czvpf", tarball, indir]
                            do.run(cmd, "Compress multiqc inputs: %s" % indir)
                    samples[0]["summary"]["multiqc"]["secondary"].append(tarball)

    if any([cwlutils.is_cwl_run(d) for d in samples]):
        samples = _add_versions(samples)

    return [[data] for data in samples]

def _add_versions(samples):
    """Add tool and data versions to the summary.
    """
    samples[0]["versions"] = {"tools": programs.write_versions(samples[0]["dirs"], samples[0]["config"]),
                              "data": provenancedata.write_versions(samples[0]["dirs"], samples)}
    return samples

def _summarize_inputs(samples, out_dir):
    """Summarize inputs for MultiQC reporting in display.
    """
    logger.info("summarize target information")
    if samples[0].get("analysis", "").lower() in ["variant", "variant2"]:
        metrics_dir = utils.safe_makedir(os.path.join(out_dir, "report", "metrics"))
        samples = _merge_target_information(samples, metrics_dir)

    logger.info("summarize fastqc")
    out_dir = utils.safe_makedir(os.path.join(out_dir, "report", "fastqc"))
    with utils.chdir(out_dir):
        _merge_fastqc(samples)

    preseq_samples = [s for s in samples if tz.get_in(["config", "algorithm", "preseq"], s)]
    if preseq_samples:
        logger.info("summarize preseq")
        out_dir = utils.safe_makedir(os.path.join(out_dir, "report", "preseq"))
        with utils.chdir(out_dir):
            _merge_preseq(preseq_samples)
    return samples

def _save_uploaded_data_json(samples, data_json_work, out_dir):
    """ Fixes all absolute work-rooted paths to relative final-rooted paths
    """
    if not utils.file_exists(data_json_work):
        return None

    upload_path_mapping = dict()
    for sample in samples:
        upload_path_mapping.update(get_all_upload_paths_from_sample(sample))
    if not upload_path_mapping:
        return data_json_work

    with io.open(data_json_work, encoding="utf-8") as f:
        data = json.load(f, object_pairs_hook=OrderedDict)
    upload_base = samples[0]["upload"]["dir"]
    data = walk_json(data, lambda s: _work_path_to_rel_final_path(s, upload_path_mapping, upload_base))

    data_json_final = os.path.join(out_dir, "multiqc_data_final.json")
    with io.open(data_json_final, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4)
    return data_json_final

def _save_uploaded_file_list(samples, file_list_work, out_dir):
    """ Fixes all absolute work-rooted paths to relative final-rooted paths

    For CWL, prepare paths relative to output directory.
    """
    if not utils.file_exists(file_list_work):
        return None

    if any([cwlutils.is_cwl_run(d) for d in samples]):
        upload_paths = []
        with open(file_list_work) as f:
            for p in (l.strip() for l in f.readlines() if os.path.exists(l.strip())):
                if p.startswith(out_dir):
                    upload_paths.append(p.replace(out_dir + "/", ""))
    else:
        upload_path_mapping = dict()
        for sample in samples:
            upload_path_mapping.update(get_all_upload_paths_from_sample(sample))
        if not upload_path_mapping:
            return None

        with open(file_list_work) as f:
            paths = [l.strip() for l in f.readlines() if os.path.exists(l.strip())]
        upload_paths = [p for p in [
            _work_path_to_rel_final_path(path, upload_path_mapping, samples[0]["upload"]["dir"])
            for path in paths
        ] if p]
        if not upload_paths:
            return None

    file_list_final = os.path.join(out_dir, "list_files_final.txt")
    with open(file_list_final, "w") as f:
        for path in upload_paths:
            f.write(path + '\n')
    return file_list_final

def _work_path_to_rel_final_path(path, upload_path_mapping, upload_base_dir):
    """ Check if `path` is a work-rooted path, and convert to a relative final-rooted path
    """
    if not path or not isinstance(path, str):
        return path
    upload_path = None

    # First, check in the mapping: if it's there is a direct reference and
    # it's a file, we immediately return it (saves lots of iterations)
    if upload_path_mapping.get(path) is not None and os.path.isfile(path):
        upload_path = upload_path_mapping[path]
    else:
        # Not a file: check for elements in the mapping that contain
        # it
        paths_to_check = [key for key in upload_path_mapping
                          if path.startswith(key)]

        if paths_to_check:
            for work_path in paths_to_check:
                if os.path.isdir(work_path):
                    final_path = upload_path_mapping[work_path]
                    upload_path = path.replace(work_path, final_path)
                    break

    if upload_path is not None:
        return os.path.relpath(upload_path, upload_base_dir)
    else:
        return None

def _one_exists(input_files):
    """
    at least one file must exist for multiqc to run properly
    """
    for f in input_files:
        if os.path.exists(f):
            return True
    return False

def _get_input_files(samples, base_dir, tx_out_dir):
    """Retrieve input files, keyed by sample and QC method name.

    Stages files into the work directory to ensure correct names for
    MultiQC sample assessment when running with CWL.
    """
    in_files = collections.defaultdict(list)
    for data in samples:
        sum_qc = tz.get_in(["summary", "qc"], data, {})
        if sum_qc in [None, "None"]:
            sum_qc = {}
        elif isinstance(sum_qc, six.string_types):
            sum_qc = {dd.get_algorithm_qc(data)[0]: sum_qc}
        elif not isinstance(sum_qc, dict):
            raise ValueError("Unexpected summary qc: %s" % sum_qc)
        for program, pfiles in sum_qc.items():
            if isinstance(pfiles, dict):
                pfiles = [pfiles["base"]] + pfiles.get("secondary", [])
            # CWL: presents output files as single file plus associated secondary files
            elif isinstance(pfiles, six.string_types):
                if os.path.exists(pfiles):
                    pfiles = [os.path.join(basedir, f) for basedir, subdir, filenames in os.walk(os.path.dirname(pfiles)) for f in filenames]
                else:
                    pfiles = []
            in_files[(dd.get_sample_name(data), program)].extend(pfiles)
    staged_files = []
    for (sample, program), files in in_files.items():
        cur_dir = utils.safe_makedir(os.path.join(base_dir, "inputs", sample, program))
        for f in files:
            if _check_multiqc_input(f) and _is_good_file_for_multiqc(f):
                if _in_temp_directory(f) or any([cwlutils.is_cwl_run(d) for d in samples]):
                    staged_f = os.path.join(cur_dir, os.path.basename(f))
                    shutil.copy(f, staged_f)
                    staged_files.append(staged_f)
                else:
                    staged_files.append(f)
    staged_files.extend(get_qsig_multiqc_files(samples))
    # Back compatible -- to migrate to explicit specifications in input YAML
    if not any([cwlutils.is_cwl_run(d) for d in samples]):
        staged_files += ["trimmed", "htseq-count/*summary"]
        # Add in created target_info file
        if os.path.isfile(os.path.join(base_dir, "report", "metrics", "target_info.yaml")):
            staged_files += [os.path.join(base_dir, "report", "metrics", "target_info.yaml")]
    return sorted(list(set(staged_files)))

def _in_temp_directory(f):
    return any(x.startswith("tmp") for x in f.split("/"))

def _get_batches(data):
    batches = dd.get_batch(data) or dd.get_sample_name(data)
    if not isinstance(batches, (list, tuple)):
        batches = [batches]
    return batches

def _group_by_sample_and_batch(samples):
    """Group samples split by QC method back one per sample-batch.
    """
    out = collections.defaultdict(list)
    for data in samples:
        out[(dd.get_sample_name(data), dd.get_align_bam(data), tuple(_get_batches(data)))].append(data)
    return [xs[0] for xs in out.values()]

def _create_list_file(paths, out_file):
    with open(out_file, "w") as f:
        for path in paths:
            f.write(path + '\n')
    return out_file

def _create_config_file(out_dir, samples):
    """Provide configuration file hiding duplicate columns.

    Future entry point for providing top level configuration of output reports.
    """
    out_file = os.path.join(out_dir, "multiqc_config.yaml")
    out = {"table_columns_visible": dict()}

    # Avoid duplicated bcbio columns with qualimap
    if any(("qualimap" in dd.get_tools_on(d) or "qualimap_full" in dd.get_tools_on(d)) for d in samples):
        # Hiding metrics duplicated by Qualimap
        out["table_columns_visible"]["bcbio"] = {"Average_insert_size": False}
        out["table_columns_visible"]["FastQC"] = {"percent_gc": False}

        # Setting up thresholds for Qualimap depth cutoff calculations, based on sample avg depths
        avg_depths = [tz.get_in(["summary", "metrics", "Avg_coverage"], s) for s in samples]
        avg_depths = [x for x in avg_depths if x]
        # Picking all thresholds up to the highest sample average depth
        thresholds = [t for t in coverage.DEPTH_THRESHOLDS if not avg_depths or t <= max(avg_depths)]
        # ...plus one more
        if len(thresholds) < len(coverage.DEPTH_THRESHOLDS):
            thresholds.append(coverage.DEPTH_THRESHOLDS[len(thresholds)])

        # Showing only thresholds surrounding any of average depths
        thresholds_hidden = []
        for i, t in enumerate(thresholds):
            if t > 20:  # Not hiding anything below 20x
                if any(thresholds[i-1] <= c < thresholds[i] for c in avg_depths if c and i-1 >= 0) or \
                   any(thresholds[i] <= c < thresholds[i+1] for c in avg_depths if c and i+1 < len(thresholds)):
                    pass
                else:
                    thresholds_hidden.append(t)

        # Hide coverage unless running full qualimap, downsampled inputs are confusing
        if not any(("qualimap_full" in dd.get_tools_on(d)) for d in samples):
            thresholds_hidden = thresholds + thresholds_hidden
            thresholds_hidden.sort()
            thresholds = []
        out['qualimap_config'] = {
            'general_stats_coverage': [str(t) for t in thresholds],
            'general_stats_coverage_hidden': [str(t) for t in thresholds_hidden]}

    # Avoid confusing peddy outputs, sticking to ancestry and sex prediction
    out["table_columns_visible"]["Peddy"] = {"family_id": False, "sex_het_ratio": False,
                                             "error_sex_check": False}

    # Setting the module order
    module_order = []
    module_order.extend([
        "bcbio",
        "samtools",
        "goleft_indexcov",
        "peddy"
    ])
    out['bcftools'] = {'write_separate_table': True}
    # if germline calling was performed:
    if any("germline" in (get_active_vcinfo(s) or {}) or  # tumor-only somatic with germline extraction
           dd.get_phenotype(s) == "germline" or           # or paired somatic with germline calling for normal
           _has_bcftools_germline_stats(s)                # CWL organized statistics
           for s in samples):
        # Split somatic and germline variant stats into separate multiqc submodules,
        # with somatic going into General Stats, and germline going into a separate table:
        module_order.extend([{
            'bcftools': {
                'name': 'Bcftools (somatic)',
                'info': 'Bcftools stats for somatic variant calls only.',
                'path_filters': ['*_bcftools_stats.txt'],
                'write_general_stats': True,
            }},
            {'bcftools': {
                'name': 'Bcftools (germline)',
                'info': 'Bcftools stats for germline variant calls only.',
                'path_filters': ['*_bcftools_stats_germline.txt'],
                'write_general_stats': False
            }},
        ])
    else:
        module_order.append("bcftools")
    module_order.extend([
        "salmon",
        "picard",
        "qualimap",
        "snpeff",
        "fastqc",
        "preseq",
    ])
    out["module_order"] = module_order

    preseq_samples = [s for s in samples if tz.get_in(["config", "algorithm", "preseq"], s)]
    if preseq_samples:
        out["preseq"] = _make_preseq_multiqc_config(preseq_samples)

    with open(out_file, "w") as out_handle:
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file

def _has_bcftools_germline_stats(data):
    """Check for the presence of a germline stats file, CWL compatible.
    """
    stats_file = tz.get_in(["summary", "qc"], data)
    if isinstance(stats_file, dict):
        stats_file = tz.get_in(["variants", "base"], stats_file)
    if not stats_file:
        stats_file = ""
    return stats_file.find("bcftools_stats_germline") > 0

def _check_multiqc_input(path):
    """Check if file exists, and return empty if it doesn't"""
    if utils.file_exists(path):
        return path

# ## report and coverage

def _is_good_file_for_multiqc(fpath):
    """Returns False if the file is binary or image."""
    # Use mimetypes to exclude binary files where possible
    (ftype, encoding) = mimetypes.guess_type(fpath)
    if encoding is not None:
        return False
    if ftype is not None and ftype.startswith('image'):
        return False
    return True

def _parse_disambiguate(disambiguatestatsfilename):
    """Parse disambiguation stats from given file.
    """
    disambig_stats = [0, 0, 0]
    with open(disambiguatestatsfilename, "r") as in_handle:
        for i, line in enumerate(in_handle):
            fields = line.strip().split("\t")
            if i == 0:
                assert fields == ['sample', 'unique species A pairs', 'unique species B pairs', 'ambiguous pairs']
            else:
                disambig_stats = [x + int(y) for x, y in zip(disambig_stats, fields[1:])]
    return disambig_stats

def _add_disambiguate(sample):
    # check if disambiguation was run
    if "disambiguate" in sample:
        if utils.file_exists(sample["disambiguate"]["summary"]):
            disambigStats = _parse_disambiguate(sample["disambiguate"]["summary"])
            sample["summary"]["metrics"]["Disambiguated %s reads" % str(sample["genome_build"])] = disambigStats[0]
            disambigGenome = (sample["config"]["algorithm"]["disambiguate"][0]
                              if isinstance(sample["config"]["algorithm"]["disambiguate"], (list, tuple))
                              else sample["config"]["algorithm"]["disambiguate"])
            sample["summary"]["metrics"]["Disambiguated %s reads" % disambigGenome] = disambigStats[1]
            sample["summary"]["metrics"]["Disambiguated ambiguous reads"] = disambigStats[2]
    return sample

def _fix_duplicated_rate(dt):
    """Get RNA duplicated rate if exists and replace by samtools metric"""
    if "Duplication_Rate_of_Mapped" in dt:
        dt["Duplicates_pct"] = 100.0 * dt["Duplication_Rate_of_Mapped"]
    return dt

def _merge_metrics(samples, out_dir):
    """Merge metrics from multiple QC steps
    """
    logger.info("summarize metrics")
    out_dir = utils.safe_makedir(os.path.join(out_dir, "report", "metrics"))
    sample_metrics = collections.defaultdict(dict)
    for s in samples:
        s = _add_disambiguate(s)
        m = tz.get_in(['summary', 'metrics'], s)
        if isinstance(m, six.string_types):
            m = json.loads(m)
        if m:
            for me in list(m.keys()):
                if isinstance(m[me], list) or isinstance(m[me], dict) or isinstance(m[me], tuple):
                    m.pop(me, None)
            sample_metrics[dd.get_sample_name(s)].update(m)
    out = []
    for sample_name, m in sample_metrics.items():
        sample_file = os.path.join(out_dir, "%s_bcbio.txt" % sample_name)
        with file_transaction(samples[0], sample_file) as tx_out_file:
            dt = pd.DataFrame(m, index=['1'])
            dt.columns = [k.replace(" ", "_").replace("(", "").replace(")", "") for k in dt.columns]
            dt['sample'] = sample_name
            dt['rRNA_rate'] = m.get('rRNA_rate', "NA")
            dt['RiP_pct'] = "%.3f" % (int(m.get("RiP", 0)) / float(m.get("Total_reads", 1)) * 100)
            dt = _fix_duplicated_rate(dt)
            dt.transpose().to_csv(tx_out_file, sep="\t", header=False)
        out.append(sample_file)
    return out

def _merge_fastqc(samples):
    """
    merge all fastqc samples into one by module
    """
    fastqc_list = collections.defaultdict(list)
    seen = set()
    for data in samples:
        name = dd.get_sample_name(data)
        if name in seen:
            continue
        seen.add(name)
        fns = glob.glob(os.path.join(dd.get_work_dir(data), "qc", dd.get_sample_name(data), "fastqc") + "/*")
        for fn in fns:
            if fn.endswith("tsv"):
                metric = os.path.basename(fn)
                fastqc_list[metric].append([name, fn])

    for metric in fastqc_list:
        dt_by_sample = []
        for fn in fastqc_list[metric]:
            dt = pd.read_csv(fn[1], sep="\t")
            dt['sample'] = fn[0]
            dt_by_sample.append(dt)
        dt = utils.rbind(dt_by_sample)
        dt.to_csv(metric, sep="\t", index=False, mode ='w')
    return samples

def _merge_preseq(samples):
    metrics = [utils.get_in(s, ("summary", "metrics")) for s in samples]
    real_counts_file = os.path.abspath(os.path.join("preseq_real_counts.txt"))
    with file_transaction(samples[0], real_counts_file) as tx_out_file:
        with open(tx_out_file, "w") as f:
            for s, m in zip(samples, metrics):
                line = dd.get_sample_name(s) + "\t" + str(m["Preseq_read_count"])
                if m.get("Preseq_unique_count") is not None:
                    line += "\t" + str(m["Preseq_unique_count"])
                line += "\n"
                f.write(line)
    samples[0]["summary"]["qc"]["preseq"]["secondary"] = [real_counts_file]

def _make_preseq_multiqc_config(samples):
    metrics = [utils.get_in(s, ("summary", "metrics")) for s in samples]
    out = {"read_length": float(np.median([m["Preseq_read_length"] for m in metrics]))}

    genome_sizes = list(set(m["Preseq_genome_size"] for m in metrics))
    if len(genome_sizes) == 1:
        out["genome_size"] = genome_sizes[0]

    return out

def _merge_target_information(samples, metrics_dir):
    out_file = os.path.abspath(os.path.join(metrics_dir, "target_info.yaml"))
    if utils.file_exists(out_file):
        return samples

    genomes = set(dd.get_genome_build(data) for data in samples)
    coverage_beds = set(dd.get_coverage(data) for data in samples)
    original_variant_regions = set(dd.get_variant_regions_orig(data) for data in samples)

    data = samples[0]
    info = {}

    # Reporting in MultiQC only if the genome is the same across all samples
    if len(genomes) == 1:
        info["genome_info"] = {
            "name": dd.get_genome_build(data),
            "size": sum([c.size for c in ref.file_contigs(dd.get_ref_file(data), data["config"])]),
        }

    # Reporting in MultiQC only if the target is the same across all samples
    vcr_orig = None
    if len(original_variant_regions) == 1 and list(original_variant_regions)[0] is not None:
        vcr_orig = list(original_variant_regions)[0]
        vcr_clean = bedutils.clean_file(vcr_orig, data)
        info["variants_regions_info"] = {
            "bed": vcr_orig,
            "size": sum(len(x) for x in pybedtools.BedTool(dd.get_variant_regions_merged(data))),
            "regions": pybedtools.BedTool(vcr_clean).count(),
        }
        gene_num = annotate.count_genes(vcr_clean, data)
        if gene_num is not None:
            info["variants_regions_info"]["genes"] = gene_num
    else:
        info["variants_regions_info"] = {
            "bed": "callable regions",
        }
    # Reporting in MultiQC only if the target is the same across samples
    if len(coverage_beds) == 1:
        cov_bed = list(coverage_beds)[0]
        if cov_bed not in [None, "None"]:
            if vcr_orig and vcr_orig == cov_bed:
                info["coverage_bed_info"] = info["variants_regions_info"]
            else:
                clean_bed = bedutils.clean_file(cov_bed, data, prefix="cov-", simple=True)
                info["coverage_bed_info"] = {
                    "bed": cov_bed,
                    "size": pybedtools.BedTool(cov_bed).total_coverage(),
                    "regions": pybedtools.BedTool(clean_bed).count(),
                }
                gene_num = annotate.count_genes(clean_bed, data)
                if gene_num is not None:
                    info["coverage_bed_info"]["genes"] = gene_num
        else:
            info["coverage_bed_info"] = info["variants_regions_info"]

    coverage_intervals = set(data["config"]["algorithm"]["coverage_interval"] for data in samples)
    if len(coverage_intervals) == 1:
        info["coverage_interval"] = list(coverage_intervals)[0]

    if info:
        with open(out_file, "w") as out_handle:
            yaml.safe_dump(info, out_handle)

    return samples
