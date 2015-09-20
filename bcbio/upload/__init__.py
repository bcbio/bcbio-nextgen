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
    elif analysis in ["RNA-seq"]:
        return _get_files_rnaseq(sample)
    elif analysis.lower() in ["smallrna-seq"]:
        return _get_files_srnaseq(sample)
    elif analysis.lower() in ["chip-seq"]:
        return _get_files_chipseq(sample)
    elif analysis.lower() in ["sailfish"]:
        return _get_files_sailfish(sample)
    else:
        return []

def _get_files_sailfish(sample):
    out = []
    out.append({"path": sample["sailfish_dir"],
                "type": "directory",
                "ext": "sailfish"})
    return _add_meta(out, sample)

def _get_files_rnaseq(sample):
    out = []
    algorithm = sample["config"]["algorithm"]
    out = _maybe_add_summary(algorithm, sample, out)
    out = _maybe_add_alignment(algorithm, sample, out)
    out = _maybe_add_disambiguate(algorithm, sample, out)
    out = _maybe_add_counts(algorithm, sample, out)
    out = _maybe_add_cufflinks(algorithm, sample, out)
    out = _maybe_add_oncofuse(algorithm, sample, out)
    out = _maybe_add_rnaseq_variant_file(algorithm, sample, out)
    return _add_meta(out, sample)

def _get_files_srnaseq(sample):
    out = []
    algorithm = sample["config"]["algorithm"]
    out = _maybe_add_summary(algorithm, sample, out)
    out = _maybe_add_trimming(algorithm, sample, out)
    out = _maybe_add_seqbuster(algorithm, sample, out)
    return _add_meta(out, sample)

def _get_files_chipseq(sample):
    out = []
    algorithm = sample["config"]["algorithm"]
    out = _maybe_add_summary(algorithm, sample, out)
    out = _maybe_add_alignment(algorithm, sample, out)
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
    return out

def _maybe_add_sv(algorithm, sample, out):
    if sample.get("align_bam") is not None and sample.get("sv"):
        for svcall in sample["sv"]:
            for key in ["vrn_file", "cnr", "cns", "cnr_bed", "cnr_bedgraph", "seg",
                        "gainloss", "segmetrics", "vrn_bed", "vrn_bedpe"]:
                out.extend(_get_variant_file(svcall, (key,)))
            if "plot" in svcall:
                for plot_name, fname in svcall["plot"].items():
                    ext = os.path.splitext(fname)[-1].replace(".", "")
                    out.append({"path": fname,
                                "type": ext,
                                "ext": "%s-%s" % (svcall["variantcaller"], plot_name),
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

def _get_variant_file(x, key):
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
                        "ext": x["variantcaller"],
                        "variantcaller": x["variantcaller"]})
            if utils.file_exists(fname + ".tbi"):
                out.append({"path": fname + ".tbi",
                            "type": "vcf.gz.tbi",
                            "index": True,
                            "ext": x["variantcaller"],
                            "variantcaller": x["variantcaller"]})
        elif fname.endswith((".vcf", ".bed", ".bedpe", ".bedgraph", ".cnr", ".cns", ".cnn", ".txt")):
            ftype = utils.splitext_plus(fname)[-1][1:]
            if ftype == "txt":
                ftype = fname.split("-")[-1]
            out.append({"path": fname,
                        "type": ftype,
                        "ext": x["variantcaller"],
                        "variantcaller": x["variantcaller"]})
    return out

def _maybe_add_summary(algorithm, sample, out):
    out = []
    if "summary" in sample:
        if sample["summary"].get("pdf"):
            out.append({"path": sample["summary"]["pdf"],
                       "type": "pdf",
                       "ext": "summary"})
        if sample["summary"].get("qc"):
            out.append({"path": sample["summary"]["qc"],
                        "type": "directory",
                        "ext": "qc"})
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

def _maybe_add_counts(algorithm, sample, out):
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
    fn = sample["seqbuster"]
    if utils.file_exists(fn):
        out.append({"path": fn,
                    "type": "counts",
                    "ext": "ready"})
        return out

def _has_alignment_file(algorithm, sample):
    return (((algorithm.get("aligner") or algorithm.get("realign")
              or algorithm.get("recalibrate") or algorithm.get("bam_clean")
              or algorithm.get("mark_duplicates")) and
              algorithm.get("merge_bamprep", True)) and
              sample.get("work_bam") is not None)

# ## File information from full project

def _get_files_project(sample, upload_config):
    """Retrieve output files associated with an entire analysis project.
    """
    out = [{"path": sample["provenance"]["programs"]}]
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

    if sample.get("seqcluster", None):
        out.append({"path": sample["seqcluster"],
                    "type": "directory", "ext": "seqcluster"})

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

    if dd.get_combined_counts(sample):
        out.append({"path": dd.get_combined_counts(sample)})
    if dd.get_annotated_combined_counts(sample):
        out.append({"path": dd.get_annotated_combined_counts(sample)})
    if dd.get_combined_fpkm(sample):
        out.append({"path": dd.get_combined_fpkm(sample)})
    if dd.get_combined_fpkm_isoform(sample):
        out.append({"path": dd.get_combined_fpkm_isoform(sample)})
    if dd.get_assembled_gtf(sample):
        out.append({"path": dd.get_assembled_gtf(sample)})
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

    return _add_meta(out, config=upload_config)
