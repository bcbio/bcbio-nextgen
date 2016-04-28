"""Handle extraction of final files from processing pipelines into storage.
"""
import datetime
import os

import toolz as tz

from bcbio import log, utils
from bcbio.upload import shared, filesystem, galaxy, s3
from bcbio.pipeline import run_info
import bcbio.pipeline.datadict as dd

_approaches = {"filesystem": filesystem,
               "galaxy": galaxy,
               "s3": s3}

def project_from_sample(sample):
    upload_config = sample.get("upload")
    if upload_config:
        approach = _approaches[upload_config.get("method", "filesystem")]
        for finfo in _get_files_project(sample, upload_config):
            approach.update_file(finfo, None, upload_config)
    return [[sample]]

def from_sample(sample):
    """Upload results of processing from an analysis pipeline sample.
    """
    upload_config = sample.get("upload")
    if upload_config:
        approach = _approaches[upload_config.get("method", "filesystem")]
        for finfo in _get_files(sample):
            approach.update_file(finfo, sample, upload_config)
    return [[sample]]

# ## File information from sample

def _get_files(sample):
    """Retrieve files for the sample, dispatching by analysis type.

    Each file is a dictionary containing the path plus associated
    metadata about the file and pipeline versions.
    """
    analysis = sample.get("analysis")
    if analysis.lower() in ["variant", "snp calling", "variant2", "standard"]:
        return _get_files_variantcall(sample)
    elif analysis.lower() in ["rna-seq", "fastrna-seq"]:
        return _get_files_rnaseq(sample)
    elif analysis.lower() in ["smallrna-seq"]:
        return _get_files_srnaseq(sample)
    elif analysis.lower() in ["chip-seq"]:
        return _get_files_chipseq(sample)
    else:
        return []

def _get_files_rnaseq(sample):
    out = []
    algorithm = sample["config"]["algorithm"]
    out = _maybe_add_summary(algorithm, sample, out)
    out = _maybe_add_alignment(algorithm, sample, out)
    out = _maybe_add_transcriptome_alignment(sample, out)
    out = _maybe_add_disambiguate(algorithm, sample, out)
    out = _maybe_add_counts(algorithm, sample, out)
    out = _maybe_add_cufflinks(algorithm, sample, out)
    out = _maybe_add_oncofuse(algorithm, sample, out)
    out = _maybe_add_rnaseq_variant_file(algorithm, sample, out)
    out = _maybe_add_sailfish_files(algorithm, sample, out)
    out = _maybe_add_salmon_files(algorithm, sample, out)
    return _add_meta(out, sample)

def _get_files_srnaseq(sample):
    out = []
    algorithm = sample["config"]["algorithm"]
    out = _maybe_add_summary(algorithm, sample, out)
    out = _maybe_add_trimming(algorithm, sample, out)
    out = _maybe_add_seqbuster(algorithm, sample, out)
    out = _maybe_add_trna(algorithm, sample, out)
    return _add_meta(out, sample)

def _get_files_chipseq(sample):
    out = []
    algorithm = sample["config"]["algorithm"]
    out = _maybe_add_summary(algorithm, sample, out)
    out = _maybe_add_alignment(algorithm, sample, out)
    out = _maybe_add_peaks(algorithm, sample, out)
    return _add_meta(out, sample)

def _add_meta(xs, sample=None, config=None):
    out = []
    for x in xs:
        x["mtime"] = shared.get_file_timestamp(x["path"])
        if sample and "sample" not in x:
            if isinstance(sample["name"], (tuple, list)):
                name = sample["name"][-1]
            else:
                name = "%s-%s" % (sample["name"],
                                  run_info.clean_name(sample["description"]))
            x["sample"] = name
        if config:
            if "fc_name" in config and "fc_date" in config:
                x["run"] = "%s_%s" % (config["fc_date"], config["fc_name"])
            else:
                x["run"] = "project_%s" % datetime.datetime.now().strftime("%Y-%m-%d")
        out.append(x)
    return out

def _get_files_variantcall(sample):
    """Return output files for the variant calling pipeline.
    """
    out = []
    algorithm = sample["config"]["algorithm"]
    out = _maybe_add_summary(algorithm, sample, out)
    out = _maybe_add_alignment(algorithm, sample, out)
    out = _maybe_add_disambiguate(algorithm, sample, out)
    out = _maybe_add_variant_file(algorithm, sample, out)
    out = _maybe_add_sv(algorithm, sample, out)
    out = _maybe_add_hla(algorithm, sample, out)
    out = _maybe_add_heterogeneity(algorithm, sample, out)
    out = _maybe_add_validate(algorithm, sample, out)
    return _add_meta(out, sample)

def _maybe_add_validate(algorith, sample, out):
    for i, plot in enumerate(tz.get_in(("validate", "grading_plots"), sample, [])):
        ptype = os.path.splitext(plot)[-1].replace(".", "")
        out.append({"path": plot,
                    "type": ptype,
                    "ext": "validate%s" % ("" if i == 0 else "-%s" % (i + 1))})
    return out

def _maybe_add_rnaseq_variant_file(algorithm, sample, out):
    if sample.get("vrn_file"):
        out.append({"path": sample.get("vrn_file"),
                    "type": "vcf",
                    "ext": "vcf"})
    return out

def _maybe_add_variant_file(algorithm, sample, out):
    if sample.get("align_bam") is not None and sample.get("vrn_file"):
        for x in sample["variants"]:
            if not _sample_variant_file_in_population(x):
                out.extend(_get_variant_file(x, ("vrn_file",)))
            if x.get("bed_file"):
                out.append({"path": x["bed_file"],
                            "type": "bed",
                            "ext": "%s-callregions" % x["variantcaller"],
                            "variantcaller": x["variantcaller"]})
            if x.get("vrn_stats"):
                for extra, fname in x["vrn_stats"].items():
                    ext = utils.splitext_plus(fname)[-1].replace(".", "")
                    out.append({"path": fname,
                                "type": ext,
                                "ext": "%s-%s" % (x["variantcaller"], extra),
                                "variantcaller": x["variantcaller"]})
            if x.get("germline") and os.path.exists(x["germline"]):
                out.extend(_get_variant_file(x, ("germline",), "-germline"))
    return out

def _maybe_add_hla(algorithm, sample, out):
    if sample.get("align_bam") is not None and sample.get("hla") and "call_file" in sample["hla"]:
        out.append({"path": sample["hla"]["call_file"],
                    "type": "csv",
                    "ext": "hla-%s" % (sample["hla"]["hlacaller"])})
    return out

def _maybe_add_heterogeneity(algorithm, sample, out):
    for hetinfo in sample.get("heterogeneity", []):
        report = hetinfo.get("report")
        if report and os.path.exists(report):
            out.append({"path": report,
                        "type": utils.splitext_plus(report)[-1].replace(".", "").replace("-", ""),
                        "ext": "%s-report" % (hetinfo["caller"])})
        for plot_type, plot_file in hetinfo.get("plots", {}).items():
            if plot_file and os.path.exists(plot_file):
                out.append({"path": plot_file,
                            "type": utils.splitext_plus(plot_file)[-1].replace(".", ""),
                            "ext": "%s-%s-plot" % (hetinfo["caller"], plot_type)})
    return out

def _maybe_add_sv(algorithm, sample, out):
    if sample.get("align_bam") is not None and sample.get("sv"):
        for svcall in sample["sv"]:
            for key in ["vrn_file", "cnr", "cns", "seg", "gainloss",
                        "segmetrics", "vrn_bed", "vrn_bedpe"]:
                out.extend(_get_variant_file(svcall, (key,)))
            if "plot" in svcall:
                for plot_name, fname in svcall["plot"].items():
                    ext = os.path.splitext(fname)[-1].replace(".", "")
                    out.append({"path": fname,
                                "type": ext,
                                "ext": "%s-%s" % (svcall["variantcaller"], plot_name),
                                "variantcaller": svcall["variantcaller"]})
            for extra in ["subclones", "contamination"]:
                svfile = svcall.get(extra)
                if svfile and os.path.exists(svfile):
                    ext = os.path.splitext(svfile)[-1].replace(".", "")
                    out.append({"path": svfile,
                                "type": ext,
                                "ext": "%s-%s" % (svcall["variantcaller"], extra),
                                "variantcaller": svcall["variantcaller"]})
        if "sv-validate" in sample:
            for vkey in ["csv", "plot", "df"]:
                vfile = tz.get_in(["sv-validate", vkey], sample)
                if vfile:
                    to_u = []
                    if isinstance(vfile, dict):
                        for svtype, fname in vfile.items():
                            to_u.append((fname, "-%s" % svtype))
                    else:
                        to_u.append((vfile, "-%s" % vkey if vkey in ["df"] else ""))
                    for vfile, ext in to_u:
                        vext = os.path.splitext(vfile)[-1].replace(".", "")
                        out.append({"path": vfile,
                                    "type": vext,
                                    "ext": "sv-validate%s" % ext})
    return out

def _sample_variant_file_in_population(x):
    """Check if a sample file is the same as the population file.

    This is true for batches where we don't extract into samples.
    '"""
    if "population" in x:
        a = _get_variant_file(x, ("population", "vcf"))
        b = _get_variant_file(x, ("vrn_file",))
        if a and b and len(a) > 0 and len(b) > 0 and os.path.getsize(a[0]["path"]) == os.path.getsize(b[0]["path"]):
            return True
    return False

def _get_variant_file(x, key, suffix=""):
    """Retrieve VCF file with the given key if it exists, handling bgzipped.
    """
    out = []
    fname = utils.get_in(x, key)
    upload_key = list(key)
    upload_key[-1] = "do_upload"
    do_upload = tz.get_in(tuple(upload_key), x, True)
    if fname and do_upload:
        if fname.endswith(".vcf.gz"):
            out.append({"path": fname,
                        "type": "vcf.gz",
                        "ext": "%s%s" % (x["variantcaller"], suffix),
                        "variantcaller": x["variantcaller"]})
            if utils.file_exists(fname + ".tbi"):
                out.append({"path": fname + ".tbi",
                            "type": "vcf.gz.tbi",
                            "index": True,
                            "ext": "%s%s" % (x["variantcaller"], suffix),
                            "variantcaller": x["variantcaller"]})
        elif fname.endswith((".vcf", ".bed", ".bedpe", ".bedgraph", ".cnr", ".cns", ".cnn", ".txt", ".tsv")):
            ftype = utils.splitext_plus(fname)[-1][1:]
            if ftype == "txt":
                ftype = fname.split("-")[-1]
            out.append({"path": fname,
                        "type": ftype,
                        "ext": "%s%s" % (x["variantcaller"], suffix),
                        "variantcaller": x["variantcaller"]})
    return out

def _maybe_add_sailfish_files(algorithm, sample, out):
    if dd.get_sailfish_dir(sample):
        out.append({"path": dd.get_sailfish_dir(sample),
                    "type": "directory",
                    "ext": "sailfish"})
    return out

def _maybe_add_salmon_files(algorithm, sample, out):
    salmon_dir = os.path.join(dd.get_work_dir(sample), "salmon", dd.get_sample_name(sample), "quant")
    if os.path.exists(salmon_dir):
        out.append({"path": salmon_dir,
                    "type": "directory",
                    "ext": "salmon"})
    return out

def _flatten_file_with_secondary(input, out_dir):
    """Flatten file representation with secondary indices (CWL-like)
    """
    out = []
    orig_dir = os.path.dirname(input["base"])
    for finfo in [input["base"]] + input["secondary"]:
        cur_dir = os.path.dirname(finfo)
        if cur_dir != orig_dir and cur_dir.startswith(orig_dir):
            cur_out_dir = os.path.join(out_dir, cur_dir.replace(orig_dir + "/", ""))
        else:
            cur_out_dir = out_dir
        out.append({"path": finfo, "dir": cur_out_dir})
    return out

def _maybe_add_summary(algorithm, sample, out):
    out = []
    if "summary" in sample:
        if sample["summary"].get("pdf"):
            out.append({"path": sample["summary"]["pdf"],
                        "type": "pdf",
                        "ext": "summary"})
        if sample["summary"].get("qc"):
            for program, finfo in sample["summary"]["qc"].iteritems():
                out.extend(_flatten_file_with_secondary(finfo, os.path.join("qc", program)))
        if utils.get_in(sample, ("summary", "researcher")):
            out.append({"path": sample["summary"]["researcher"],
                        "type": "tsv",
                        "sample": run_info.clean_name(utils.get_in(sample, ("upload", "researcher"))),
                        "ext": "summary"})
    return out

def _maybe_add_alignment(algorithm, sample, out):
    if _has_alignment_file(algorithm, sample):
        for (fname, ext, isplus) in [(sample.get("work_bam"), "ready", False),
                                     (utils.get_in(sample, ("work_bam-plus", "disc")), "disc", True),
                                     (utils.get_in(sample, ("work_bam-plus", "sr")), "sr", True)]:
            if fname and os.path.exists(fname):
                if fname.endswith("bam"):
                    ftype, fext = "bam", ".bai"
                elif fname.endswith("cram"):
                    ftype, fext = "cram", ".crai"
                else:
                    raise ValueError("Unexpected alignment file type %s" % fname)
                out.append({"path": fname,
                            "type": ftype,
                            "plus": isplus,
                            "ext": ext})
                if utils.file_exists(fname + fext):
                    out.append({"path": fname + fext,
                                "type": ftype + fext,
                                "plus": isplus,
                                "index": True,
                                "ext": ext})
    return out

def _maybe_add_disambiguate(algorithm, sample, out):
    if "disambiguate" in sample:
        for extra_name, fname in sample["disambiguate"].items():
            ftype = os.path.splitext(fname)[-1].replace(".", "")
            fext = ".bai" if ftype == "bam" else ""
            if fname and os.path.exists(fname):
                out.append({"path": fname,
                            "type": ftype,
                            "plus": True,
                            "ext": "disambiguate-%s" % extra_name})
                if fext and utils.file_exists(fname + fext):
                    out.append({"path": fname + fext,
                                "type": ftype + fext,
                                "plus": True,
                                "index": True,
                                "ext": "disambiguate-%s" % extra_name})
    return out

def _maybe_add_transcriptome_alignment(sample, out):
    transcriptome_bam = dd.get_transcriptome_bam(sample)
    if transcriptome_bam and utils.file_exists(transcriptome_bam):
        out.append({"path": transcriptome_bam,
                    "type": "bam",
                    "ext": "transcriptome"})
    return out

def _maybe_add_counts(algorithm, sample, out):
    if not dd.get_count_file(sample):
        return out
    out.append({"path": sample["count_file"],
             "type": "counts",
             "ext": "ready"})
    stats_file = os.path.splitext(sample["count_file"])[0] + ".stats"
    if utils.file_exists(stats_file):
        out.append({"path": stats_file,
                    "type": "count_stats",
                    "ext": "ready"})
    return out

def _maybe_add_oncofuse(algorithm, sample, out):
    if sample.get("oncofuse_file", None) is not None:
        out.append({"path": sample["oncofuse_file"],
                    "type": "oncofuse_outfile",
                    "ext": "ready"})
    return out

def _maybe_add_cufflinks(algorithm, sample, out):
    if "cufflinks_dir" in sample:
        out.append({"path": sample["cufflinks_dir"],
                    "type": "directory",
                    "ext": "cufflinks"})
    return out

def _maybe_add_trimming(algorithm, sample, out):
    fn = sample["collapse"] + "_size_stats"
    if utils.file_exists(fn):
        out.append({"path": fn,
                    "type": "trimming_stats",
                    "ext": "ready"})
        return out

def _maybe_add_seqbuster(algorithm, sample, out):
    if "seqbuster" not in sample:
        return out
    fn = sample["seqbuster"]
    if utils.file_exists(fn):
        out.append({"path": fn,
                    "type": "counts",
                    "ext": "mirbase-ready"})
    fn = sample.get("seqbuster_novel")
    if fn and utils.file_exists(fn):
        out.append({"path": fn,
                    "type": "counts",
                    "ext": "novel-ready"})
    return out

def _maybe_add_trna(algorithm, sample, out):
    if "trna" not in sample:
        return out
    fn = sample["trna"]
    if utils.file_exists(fn):
        out.append({"path": fn,
                    "type": "directory",
                    "ext": "tdrmapper"})
    return out

def _maybe_add_peaks(algorithm, sample, out):
    fns = sample.get("peaks_file", [])
    for fn in fns:
        if utils.file_exists(fn):
            name, ext = utils.splitext_plus(fn)
            caller = _get_peak_file(sample, name)
            out.append({"path": fn,
                        "type": ext,
                        "ext": caller})
    return out

def _get_peak_file(x, fn_name):
    """Get peak caller for this file name."""
    for caller in dd.get_peakcaller(x):
        if fn_name.find(caller) > -1:
            return caller
    return os.path.basename(fn_name)

def _has_alignment_file(algorithm, sample):
    return (((algorithm.get("aligner") or algorithm.get("realign")
              or algorithm.get("recalibrate") or algorithm.get("bam_clean")
              or algorithm.get("mark_duplicates"))) and
              sample.get("work_bam") is not None)

# ## File information from full project

def _get_files_project(sample, upload_config):
    """Retrieve output files associated with an entire analysis project.
    """
    out = [{"path": sample["provenance"]["programs"]}]
    if os.path.exists(tz.get_in(["provenance", "data"], sample) or ""):
        out.append({"path": sample["provenance"]["data"]})
    for fname in ["bcbio-nextgen.log", "bcbio-nextgen-commands.log"]:
        if os.path.exists(os.path.join(log.get_log_dir(sample["config"]), fname)):
            out.append({"path": os.path.join(log.get_log_dir(sample["config"]), fname),
                        "type": "external_command_log",
                        "ext": ""})

    if "summary" in sample and sample["summary"].get("project"):
        out.append({"path": sample["summary"]["project"]})
    mixup_check = tz.get_in(["summary", "mixup_check"], sample)
    if mixup_check:
        out.append({"path": sample["summary"]["mixup_check"],
                    "type": "directory", "ext": "mixup_check"})

    report = os.path.join(dd.get_work_dir(sample), "report")
    if utils.file_exists(report):
        out.append({"path": report,
                    "type": "directory", "ext": "report"})

    multiqc = tz.get_in(["summary", "multiqc"], sample)
    if multiqc:
        out.extend(_flatten_file_with_secondary(multiqc, "multiqc"))

    if sample.get("seqcluster", None):
        out.append({"path": sample["seqcluster"],
                    "type": "directory", "ext": "seqcluster"})

    if sample.get("report", None):
        out.append({"path": os.path.dirname(sample["report"]),
                    "type": "directory", "ext": "seqclusterViz"})

    for x in sample.get("variants", []):
        if "pop_db" in x:
            out.append({"path": x["pop_db"],
                        "type": "sqlite",
                        "variantcaller": x["variantcaller"]})
    for x in sample.get("variants", []):
        if "population" in x:
            pop_db = tz.get_in(["population", "db"], x)
            if pop_db:
                out.append({"path": pop_db,
                            "type": "sqlite",
                            "variantcaller": x["variantcaller"]})
            out.extend(_get_variant_file(x, ("population", "vcf")))
    for x in sample.get("variants", []):
        if x.get("validate") and x["validate"].get("grading_summary"):
            out.append({"path": x["validate"]["grading_summary"]})
            break
    if "coverage" in sample:
        cov_db = tz.get_in(["coverage", "summary"], sample)
        if cov_db:
            out.append({"path": cov_db, "type": "sqlite", "ext": "coverage"})
        all_coverage = tz.get_in(["coverage", "all"], sample)
        if all_coverage:
            out.append({"path": all_coverage, "type": "bed", "ext": "coverage"})

    if dd.get_mirna_counts(sample):
        out.append({"path": dd.get_mirna_counts(sample)})
    if dd.get_isomir_counts(sample):
        out.append({"path": dd.get_isomir_counts(sample)})
    if dd.get_novel_mirna_counts(sample):
        out.append({"path": dd.get_novel_mirna_counts(sample)})
    if dd.get_novel_isomir_counts(sample):
        out.append({"path": dd.get_novel_isomir_counts(sample)})
    if dd.get_combined_counts(sample):
        out.append({"path": dd.get_combined_counts(sample)})
    if dd.get_annotated_combined_counts(sample):
        out.append({"path": dd.get_annotated_combined_counts(sample)})
    if dd.get_combined_fpkm(sample):
        out.append({"path": dd.get_combined_fpkm(sample)})
    if dd.get_combined_fpkm_isoform(sample):
        out.append({"path": dd.get_combined_fpkm_isoform(sample)})
    if dd.get_transcript_assembler(sample):
        out.append({"path": dd.get_merged_gtf(sample)})
    if dd.get_dexseq_counts(sample):
        out.append({"path": dd.get_dexseq_counts(sample)})
    if dd.get_express_counts(sample):
        out.append({"path": dd.get_express_counts(sample)})
    if dd.get_express_fpkm(sample):
        out.append({"path": dd.get_express_fpkm(sample)})
    if dd.get_express_tpm(sample):
        out.append({"path": dd.get_express_tpm(sample)})
    if dd.get_isoform_to_gene(sample):
        out.append({"path": dd.get_isoform_to_gene(sample)})
    if dd.get_square_vcf(sample):
        out.append({"path": dd.get_square_vcf(sample)})
    if dd.get_sailfish_tidy(sample):
        out.append({"path": dd.get_sailfish_tidy(sample)})
    if dd.get_sailfish_transcript_tpm(sample):
        out.append({"path": dd.get_sailfish_transcript_tpm(sample)})
    if dd.get_sailfish_gene_tpm(sample):
        out.append({"path": dd.get_sailfish_gene_tpm(sample)})
    if dd.get_tx2gene(sample):
        out.append({"path": dd.get_tx2gene(sample)})
    return _add_meta(out, config=upload_config)
