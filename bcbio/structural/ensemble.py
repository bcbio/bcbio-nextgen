"""Combine multiple structural variation callers into single output file.

Takes a simple union approach for reporting the final set of calls, reporting
the evidence from each input.
"""
import collections
import fileinput
import os
import shutil

import toolz as tz
import vcf

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import shared
from bcbio.structural import validate
from bcbio.variation import bedutils

# ## Conversions to simplified BED files

MAX_SVSIZE = 1e6  # 1Mb maximum size from callers to avoid huge calls collapsing all structural variants
N_FILTER_CALLERS = 2  # Minimum number of callers for doing filtering of ensemble calls

def _vcf_to_bed(in_file, caller, out_file):
    if in_file and in_file.endswith((".vcf", "vcf.gz")):
        with utils.open_gzipsafe(in_file) as in_handle:
            with open(out_file, "w") as out_handle:
                for rec in vcf.Reader(in_handle, in_file):
                    if not rec.FILTER:
                        if (rec.samples[0].gt_type and
                              not (hasattr(rec.samples[0].data, "FT") and rec.samples[0].data.FT)):
                            start = rec.start - 1
                            end = int(rec.INFO.get("END", rec.start))
                            if end - start < MAX_SVSIZE:
                                out_handle.write("\t".join([rec.CHROM, str(start), str(end),
                                                            "%s_%s" % (_get_svtype(rec), caller)])
                                                 + "\n")

def _get_svtype(rec):
    try:
        return rec.INFO["SVTYPE"]
    except KeyError:
        return "-".join(str(x).replace("<", "").replace(">", "") for x in rec.ALT)

def _cnvbed_to_bed(in_file, caller, out_file):
    """Convert cn_mops CNV based bed files into flattened BED
    """
    import pybedtools
    with open(out_file, "w") as out_handle:
        for feat in pybedtools.BedTool(in_file):
            out_handle.write("\t".join([feat.chrom, str(feat.start), str(feat.end),
                                        "cnv%s_%s" % (feat.score, caller)])
                             + "\n")

def _copy_file(in_file, caller, out_file):
    shutil.copy(in_file, out_file)

CALLER_TO_BED = {"lumpy": _vcf_to_bed,
                 "delly": _vcf_to_bed,
                 "cn_mops": _cnvbed_to_bed,
                 "wham": _copy_file}

def _create_bed(call, sample, work_dir, data):
    """Create a simplified BED file from caller specific input.
    """
    out_file = os.path.join(work_dir, "%s-ensemble-%s.bed" % (sample, call["variantcaller"]))
    if call.get("vrn_file") and not utils.file_uptodate(out_file, call["vrn_file"]):
        with file_transaction(data, out_file) as tx_out_file:
            convert_fn = CALLER_TO_BED.get(call["variantcaller"])
            if convert_fn:
                convert_fn(call["vrn_file"], call["variantcaller"], tx_out_file)
    if utils.file_exists(out_file):
        return out_file

# ## Top level

def _combine_bed_by_size(input_beds, sample, work_dir, data):
    """Combine a set of BED files, breaking into individual size chunks.
    """
    import pybedtools
    out_file = os.path.join(work_dir, "%s-ensemble.bed" % sample)
    if len(input_beds) > 0:
        size_beds = []
        for e_start, e_end in validate.EVENT_SIZES:
            base, ext = os.path.splitext(out_file)
            size_out_file = "%s-%s_%s%s" % (base, e_start, e_end, ext)
            if not utils.file_exists(size_out_file):
                with file_transaction(data, size_out_file) as tx_out_file:
                    with shared.bedtools_tmpdir(data):
                        all_file = "%s-all.bed" % utils.splitext_plus(tx_out_file)[0]
                        has_regions = False
                        with open(all_file, "w") as out_handle:
                            for line in fileinput.input(input_beds):
                                chrom, start, end = line.split()[:3]
                                size = int(end) - int(start)
                                if size >= e_start and size < e_end:
                                    out_handle.write(line)
                                    has_regions = True
                        if has_regions:
                            pybedtools.BedTool(all_file).sort(stream=True)\
                              .merge(c=4, o="distinct", delim=",").saveas(tx_out_file)
            if utils.file_exists(size_out_file):
                size_beds.append(size_out_file)
        if len(size_beds) > 0:
            out_file = bedutils.combine(size_beds, out_file, data)
    return out_file

def _filter_ensemble(in_bed, data):
    """Filter ensemble set of calls, requiring calls supported by 2 callers.

    We filter only smaller size events, which seem to benefit the most since
    they have lower precision. We also check to be sure that the required
    number of callers actually called in each event, since some callers don't handle
    all event types.
    """
    max_size = max([xs[1] for xs in validate.EVENT_SIZES[:2]])
    out_file = "%s-filter%s" % utils.splitext_plus(in_bed)
    total_callers = collections.defaultdict(set)
    with open(in_bed) as in_handle:
        for line in in_handle:
            caller_strs = line.strip().split()[-1]
            for event, caller in [x.split("_") for x in caller_strs.split(",")]:
                total_callers[validate.cnv_to_event(event)].add(caller)

    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                with open(in_bed) as in_handle:
                    for line in in_handle:
                        chrom, start, end, caller_strs = line.strip().split()
                        size = int(end) - int(start)
                        callers = set([x.split("_")[-1] for x in caller_strs.split(",")])
                        events = set([validate.cnv_to_event(x.split("_")[0]) for x in caller_strs.split(",")])
                        pass_event_counts = [len(total_callers[e]) > N_FILTER_CALLERS for e in list(events)]
                        if len(callers) > 1 or size > max_size or not any(pass_event_counts):
                            out_handle.write(line)
    return out_file

def summarize(calls, data):
    """Summarize results from multiple callers into a single flattened BED file.
    """
    sample = tz.get_in(["rgnames", "sample"], data)
    work_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                               sample, "ensemble"))
    with shared.bedtools_tmpdir(data):
        input_beds = filter(lambda x: x is not None,
                            [_create_bed(c, sample, work_dir, data) for c in calls])
    if len(input_beds) > 0:
        out_file = _combine_bed_by_size(input_beds, sample, work_dir, data)
        if utils.file_exists(out_file):
            if len(input_beds) > N_FILTER_CALLERS:
                filter_file = _filter_ensemble(out_file, data)
            else:
                filter_file = out_file
            bedprep_dir = utils.safe_makedir(os.path.join(os.path.dirname(filter_file), "bedprep"))
            calls.append({"variantcaller": "ensemble",
                          "vrn_file": bedutils.clean_file(filter_file, data, bedprep_dir=bedprep_dir)})
    return calls
