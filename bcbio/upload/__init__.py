"""Handle extraction of final files from processing pipelines into storage.
"""
import datetime
import os

from bcbio.upload import shared, filesystem, galaxy, s3

_approaches = {"filesystem": filesystem,
               "galaxy": galaxy,
               "s3": s3}

def from_sample(sample):
    """Upload results of processing from an analysis pipeline sample.
    """
    upload_config = sample["info"].get("upload")
    if upload_config:
        approach = _approaches[upload_config.get("method", "filesystem")]
        for finfo in _get_files(sample):
            approach.update_file(finfo, sample["info"], upload_config)
        for finfo in _get_files_project(sample, upload_config):
            approach.update_file(finfo, None, upload_config)

# ## File information from sample

def _get_files(sample):
    """Retrieve files for the sample, dispatching by analysis type.

    Each file is a dictionary containing the path plus associated
    metadata about the file and pipeline versions.
    """
    analysis = sample["info"].get("analysis")
    if analysis in ["variant", "SNP calling", "variant2"]:
        return _get_files_variantcall(sample)
    else:
        return []

def _add_meta(xs, sample=None, config=None):
    out = []
    for x in xs:
        x["mtime"] = shared.get_file_timestamp(x["path"])
        if sample:
            x["sample"] = sample["name"][-1]
        if config:
            x["run"] = "%s_%s" % (config["fc_date"], config["fc_name"])
        out.append(x)
    return out

def _get_files_variantcall(sample):
    """Return output files for the variant calling pipeline.
    """
    out = []
    algorithm = sample["config"]["algorithm"]
    if algorithm.get("write_summary", True) and "summary" in sample:
        if sample["summary"].get("pdf"):
            out = [{"path": sample["summary"]["pdf"],
                    "type": "pdf",
                    "ext": "summary"}]
    if ((algorithm.get("aligner") or algorithm.get("realign") or algorithm.get("recalibrate"))
          and algorithm.get("merge_bamprep", True)) and sample["work_bam"] is not None:
        out.append({"path": sample["work_bam"],
                    "type": "bam",
                    "ext": "ready"})
        if os.path.exists(sample["work_bam"] + ".bai"):
            out.append({"path": sample["work_bam"] + ".bai",
                        "type": "bai",
                        "ext": "ready"})
    if sample["work_bam"] is not None:
        for x in sample["variants"]:
            out.append({"path": x["vrn_file"],
                        "type": "vcf",
                        "ext": x["variantcaller"],
                        "variantcaller": x["variantcaller"]})
            if x.get("bed_file"):
                out.append({"path": x["bed_file"],
                            "type": "bed",
                            "ext": "%s-callregions" % x["variantcaller"],
                            "variantcaller": x["variantcaller"]})
    return _add_meta(out, sample)

# ## File information from full project

def _get_files_project(sample, upload_config):
    """Retrieve output files associated with an entire analysis project.
    """
    out = [{"path": sample["info"]["provenance"]["programs"]}]
    if sample["summary"].get("project"):
        out.append({"path": sample["summary"]["project"]})
    for x in sample["variants"]:
        if "pop_db" in x:
            out.append({"path": x["pop_db"],
                        "type": "sqlite",
                        "variantcaller": x["variantcaller"]})
    for x in sample["variants"]:
        if x.get("validate", {}).get("grading_summary"):
            out.append({"path": x["validate"]["grading_summary"]})
            break
    return _add_meta(out, config=upload_config)
