"""High level summaries of samples and programs with MultiQC.

https://github.com/ewels/MultiQC
"""
import collections
import glob
import mimetypes
import os
import pandas as pd
import shutil

import pybedtools
import toolz as tz
import yaml

from bcbio import utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.log import logger
from bcbio.provenance import do
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.bam import ref
from bcbio.structural import annotate, regions

def summary(*samples):
    """Summarize all quality metrics together"""
    samples = utils.unpack_worlds(samples)
    work_dir = dd.get_work_dir(samples[0])
    multiqc = config_utils.get_program("multiqc", samples[0]["config"])
    if not multiqc:
        logger.debug("multiqc not found. Update bcbio_nextgen.py tools to fix this issue.")
    file_fapths = []
    opts = ""
    out_dir = os.path.join(work_dir, "multiqc")
    out_data = os.path.join(work_dir, "multiqc", "multiqc_data")
    out_file = os.path.join(out_dir, "multiqc_report.html")
    samples = _report_summary(samples, os.path.join(out_dir, "report"))
    for data in samples:
        for program, pfiles in tz.get_in(["summary", "qc"], data, {}).iteritems():
            if isinstance(pfiles, dict):
                pfiles = [pfiles["base"]] + pfiles["secondary"]
            elif isinstance(pfiles, basestring):
                pfiles = [pfiles]
            file_fapths.extend(pfiles)
    file_fapths.append(os.path.join(out_dir, "report", "metrics", "target_info.yaml"))
    # XXX temporary workaround until we can handle larger inputs through MultiQC
    file_fapths = list(set(file_fapths))
    # Back compatible -- to migrate to explicit specifications in input YAML
    file_fapths += ["trimmed", "htseq-count/*summary"]
    if not utils.file_exists(out_file):
        with utils.chdir(work_dir):
            file_fapths = [fpath for fpath in file_fapths if _check_multiqc_input(fpath) and _is_good_file_for_multiqc(fpath)]
            input_list_file = _create_list_file(file_fapths)
            export_tmp = ""
            if dd.get_tmp_dir(samples[0]):
                export_tmp = "export TMPDIR=%s &&" % dd.get_tmp_dir(samples[0])
            if input_list_file:
                cmd = "{export_tmp} {multiqc} -f -l {input_list_file} -o {tx_out} {opts}"
                with tx_tmpdir(data, work_dir) as tx_out:
                    do.run(cmd.format(**locals()), "Run multiqc")
                    if utils.file_exists(os.path.join(tx_out, "multiqc_report.html")):
                        shutil.move(os.path.join(tx_out, "multiqc_report.html"), out_file)
                        shutil.move(os.path.join(tx_out, "multiqc_data"), out_data)
    out = []
    for i, data in enumerate(samples):
        if i == 0:
            if utils.file_exists(out_file):
                data_files = glob.glob(os.path.join(out_dir, "multiqc_data", "*.txt"))
                data_files += glob.glob(os.path.join(out_dir, "report", "*", "*.bed"))
                data_files += glob.glob(os.path.join(out_dir, "report", "*", "*.txt"))
                data_files += glob.glob(os.path.join(out_dir, "report", "*", "*.tsv"))
                data_files += glob.glob(os.path.join(out_dir, "report", "*.R*"))
                if "summary" not in data:
                    data["summary"] = {}
                data["summary"]["multiqc"] = {"base": out_file, "secondary": data_files}
        out.append(data)
    return [[fpath] for fpath in out]

def _create_list_file(paths):
    out_file = os.path.join(os.getcwd(), "list_files.txt")
    is_any = False
    with open(out_file, "w") as outh:
        for f in paths:
            if f:
                is_any = True
                print >>outh, f
    if not is_any:
        return None
    return out_file

def _check_multiqc_input(path):
    """Check if file exists, and return empty if it doesn't"""
    if utils.file_exists(path):
        if path.find("bcfstats") == -1:
            return path

# ## report and coverage

def _is_good_file_for_multiqc(fapth):
    """Returns False if the file is binary or image."""
    # Use mimetypes to exclude binary files where possible
    (ftype, encoding) = mimetypes.guess_type(fapth)
    if encoding is not None:
        return False
    if ftype is not None and ftype.startswith('image'):
        return False
    return True

def _report_summary(samples, out_dir):
    """
    Run coverage report with bcbiocov package
    """
    try:
        import bcbreport.prepare as bcbreport
    except ImportError:
        logger.info("skipping report. No bcbreport installed.")
        return samples
    # samples = utils.unpack_worlds(samples)
    work_dir = dd.get_work_dir(samples[0])
    parent_dir = utils.safe_makedir(out_dir)
    with utils.chdir(parent_dir):
        logger.info("copy qsignature")
        qsignature_fn = os.path.join(work_dir, "qc", "qsignature", "qsignature.ma")
        if qsignature_fn:  # this need to be inside summary/qc dict
            if utils.file_exists(qsignature_fn) and not utils.file_exists("qsignature.ma"):
                shutil.copy(qsignature_fn, "bcbio_qsignature.ma")

        out_dir = utils.safe_makedir("fastqc")
        logger.info("summarize fastqc")
        with utils.chdir(out_dir):
            _merge_fastqc(samples)

        logger.info("summarize metrics")
        samples = _merge_metrics(samples)

        logger.info("summarize target information")
        samples = _merge_target_information(samples)

        out_dir = utils.safe_makedir("coverage")
        logger.info("summarize coverage")
        for data in samples:
            pfiles = tz.get_in(["summary", "qc", "coverage"], data, [])
            if isinstance(pfiles, dict):
                pfiles = [pfiles["base"]] + pfiles["secondary"]
            elif pfiles:
                pfiles = [pfiles]
            for fn in pfiles:
                if os.path.basename(fn).find("coverage_fixed") > -1:
                    utils.copy_plus(fn, os.path.join(out_dir, os.path.basename(fn)))

        out_dir = utils.safe_makedir("variants")
        logger.info("summarize variants")
        for data in samples:
            pfiles = tz.get_in(["summary", "qc", "variants"], data, [])
            if isinstance(pfiles, dict):
                pfiles = [pfiles["base"]] + pfiles["secondary"]
            elif pfiles:
                pfiles = [pfiles]
            for fn in pfiles:
                if os.path.basename(fn).find("gc-depth-parse.tsv") > -1:
                    utils.copy_plus(fn, os.path.join(out_dir, os.path.basename(fn)))
        bcbreport.report(parent_dir)
        out_report = os.path.join(parent_dir, "qc-coverage-report.html")
        if not utils.file_exists(out_report):
            rmd_file = os.path.join(parent_dir, "report-ready.Rmd")
            run_file = "%s-run.R" % (os.path.splitext(out_report)[0])
            with open(run_file, "w") as out_handle:
                out_handle.write("""library(rmarkdown)\nrender("%s")\n""" % rmd_file)
            cmd = "%s %s" % (utils.Rscript_cmd(), run_file)
            # Skip automated generation of coverage report to avoid error
            # messages. We need to generalize coverage reporting and re-include.
            # try:
            #     do.run(cmd, "Prepare coverage summary", log_error=False)
# except subprocess.CalledProcessError as msg:
            #     logger.info("Skipping generation of coverage report: %s" % (str(msg)))
            if utils.file_exists("report-ready.html"):
                shutil.move("report-ready.html", out_report)
    return samples

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
        dt["Duplicates_pct"] = dt["Duplication_Rate_of_Mapped"]
    return dt

def _merge_metrics(samples):
    """
    parse project.yaml file to get metrics for each bam
    """
    out_file = os.path.join("metrics", "metrics.tsv")
    dt_together = []
    cov = {}
    with file_transaction(out_file) as out_tx:
        for s in samples:
            sample_name = dd.get_sample_name(s)
            s = _add_disambiguate(s)
            if sample_name in cov:
                continue
            m = tz.get_in(['summary', 'metrics'], s)
            sample_file = os.path.abspath(os.path.join("metrics", "%s_bcbio.txt" % sample_name))
            if not tz.get_in(['summary', 'qc'], s):
                s['summary'] = {"qc": {}}
            if m:
                for me in m.keys():
                    if isinstance(m[me], list) or isinstance(m[me], dict) or isinstance(m[me], tuple):
                        m.pop(me, None)
                dt = pd.DataFrame(m, index=['1'])
                dt.columns = [k.replace(" ", "_").replace("(", "").replace(")", "") for k in dt.columns]
                dt['sample'] = sample_name
                dt['rRNA_rate'] = m.get('rRNA_rate', "NA")
                dt = _fix_duplicated_rate(dt)
                dt.transpose().to_csv(sample_file, sep="\t", header=False)
                dt_together.append(dt)
                s['summary']['qc'].update({'bcbio':{'base': sample_file, 'secondary': []}})
        if len(dt_together) > 0:
            dt_together = utils.rbind(dt_together)
            dt_together.to_csv(out_tx, index=False, sep="\t")

    out = []
    for s in samples:
        out.append(s)
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

def _merge_target_information(samples):
    out_file = os.path.abspath(os.path.join("metrics", "target_info.yaml"))
    if utils.file_exists(out_file):
        return samples

    genomes = set(dd.get_genome_build(data) for data in samples)
    coverage_beds = set(dd.get_coverage(data) for data in samples)
    original_variant_regions = set(dd.get_variant_regions_orig(data) for data in samples)

    data = samples[0]
    info = {}

    # Reporting in MultiQC only if the genome is the sample across samples
    if len(genomes) == 1:
        info["genome_info"] = {
            "name": dd.get_genome_build(data),
            "size": sum([c.size for c in ref.file_contigs(dd.get_ref_file(data), data["config"])]),
        }

    # Reporting in MultiQC only if the target is the sample across samples
    vcr_orig = None
    if len(original_variant_regions) == 1 and list(original_variant_regions)[0] is not None:
        vcr_orig = list(original_variant_regions)[0]
        vcr_ann = annotate.add_genes(vcr_orig, data)
        info["variants_regions_info"] = {
            "bed": vcr_orig,
            "size": sum(len(x) for x in pybedtools.BedTool(dd.get_variant_regions_merged(data))),
            "regions": pybedtools.BedTool(vcr_orig).count(),
            "genes": len(list(set(r.name for r in pybedtools.BedTool(vcr_ann) if r.name and r.name != "."))),
        }
    else:
        info["variants_regions_info"] = {
            "bed": "callable regions",
        }
    # Reporting in MultiQC only if the target is the sample across samples
    if len(coverage_beds) == 1:
        cov_bed = list(coverage_beds)[0]
        if cov_bed is not None:
            if vcr_orig and vcr_orig == cov_bed:
                info["coverage_bed_info"] = info["variants_regions_info"]
            else:
                ann_bed = annotate.add_genes(cov_bed, data)
                info["coverage_bed_info"] = {
                    "bed": cov_bed,
                    "size": pybedtools.BedTool(cov_bed).total_coverage(),
                    "regions": pybedtools.BedTool(cov_bed).count(),
                    "genes": len(list(set(r.name for r in pybedtools.BedTool(ann_bed) if r.name and r.name != "."))),
                }
        else:
            info["coverage_bed_info"] = info["variants_regions_info"]

    coverage_intervals = set(data["config"]["algorithm"]["coverage_interval"] for data in samples)
    if len(coverage_intervals) == 1:
        info["coverage_interval"] = list(coverage_intervals)[0]

    if info:
        with open(out_file, "w") as out_handle:
            yaml.safe_dump(info, out_handle)

    return samples
