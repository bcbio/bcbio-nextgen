"""Provide analysis of input files by chromosomal regions.

Handle splitting and analysis of files from chromosomal subsets separated by
no-read regions.
"""
import os

from bcbio.distributed.split import parallel_split_combine

# ## data preparation

def add_region_info(samples, regions):
    """Add reference to BED file of callable regions to each sample.
    """
    out = []
    for data in samples:
        data["config"]["algorithm"]["callable_regions"] = regions["analysis_bed"]
        out.append(data)
    return out

# ## BAM preparation

def to_safestr(region):
    if region[0] in ["nochrom", "noanalysis"]:
        return region[0]
    else:
        return "_".join([str(x) for x in region])

def _split_by_regions(regions, dirname, out_ext, in_key):
    """Split a BAM file data analysis into chromosomal regions.
    """
    def _do_work(data):
        bam_file = data[in_key]
        if bam_file is None:
            return None, []
        part_info = []
        base_out = os.path.splitext(os.path.basename(bam_file))[0]
        nowork = [["nochrom"], ["noanalysis", regions["noanalysis"]]]
        for region in regions["analysis"] + nowork:
            out_dir = os.path.join(data["dirs"]["work"], dirname, data["name"][-1], region[0])
            region_outfile = os.path.join(out_dir, "%s-%s%s" %
                                          (base_out, to_safestr(region), out_ext))
            part_info.append((region, region_outfile))
        out_file = os.path.join(data["dirs"]["work"], dirname, data["name"][-1],
                                "%s%s" % (base_out, out_ext))
        return out_file, part_info
    return _do_work

def parallel_prep_region(samples, regions, run_parallel):
    """Perform full pre-variant calling BAM prep work on regions.
    """
    file_key = "work_bam"
    split_fn = _split_by_regions(regions, "bamprep", "-prep.bam", file_key)
    # identify samples that do not need preparation -- no prep or
    # variant calling
    extras = []
    torun = []
    for data in [x[0] for x in samples]:
        data["align_bam"] = data["work_bam"]
        a = data["config"]["algorithm"]
        if (not a.get("mark_duplicates") and not a.get("recalibrate") and
              not a.get("realign", "gatk") and not a.get("variantcaller", "gatk")):
            extras.append([data])
        elif not data.get(file_key):
            extras.append([data])
        else:
            torun.append([data])
    return extras + parallel_split_combine(torun, split_fn, run_parallel,
                                           "piped_bamprep", None, file_key, ["config"])

def delayed_bamprep_merge(samples, run_parallel):
    """Perform a delayed merge on regional prepared BAM files.
    """
    needs_merge = False
    for data in samples:
        if (data[0]["config"]["algorithm"].get("merge_bamprep", True) and
              "combine" in data[0]):
            needs_merge = True
            break
    if needs_merge:
        return run_parallel("delayed_bam_merge", samples)
    else:
        return samples

# ## Utilities

def clean_sample_data(samples):
    """Clean unnecessary information from sample data, reducing size for message passing.
    """
    out = []
    for data in samples:
        data["dirs"] = {"work": data["dirs"]["work"], "galaxy": data["dirs"]["galaxy"],
                        "fastq": data["dirs"].get("fastq")}
        data["config"] = {"algorithm": data["config"]["algorithm"],
                          "resources": data["config"]["resources"]}
        for remove_attr in ["config_file", "regions", "algorithm"]:
            data.pop(remove_attr, None)
        out.append([data])
    return out
