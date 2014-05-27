"""Combine multiple structural variation callers into single output file.

Takes a simple union approach for reporting the final set of calls, reporting
the evidence from each input.
"""
import fileinput
import os

try:
    import pybedtools
except ImportError:
    pybedtools = None
import toolz as tz
import vcf

from bcbio import utils
from bcbio.distributed.transaction import file_transaction

# ## Conversions to simplified BED files

def _vcf_to_bed(in_file, caller, out_file):
    if in_file and in_file.endswith((".vcf", "vcf.gz")):
        with utils.open_gzipsafe(in_file) as in_handle:
            with open(out_file, "w") as out_handle:
                for rec in vcf.Reader(in_handle, in_file):
                    if not rec.FILTER:
                        if (rec.samples[0].gt_type and
                              not (hasattr(rec.samples[0].data, "FT") and rec.samples[0].data.FT)):
                            out_handle.write("\t".join([rec.CHROM, str(rec.start - 1),
                                                        str(rec.INFO.get("END", rec.start)),
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

    with open(out_file, "w") as out_handle:
        for feat in pybedtools.BedTool(in_file):
            out_handle.write("\t".join([feat.chrom, str(feat.start), str(feat.end),
                                        "cnv%s_%s" % (feat.score, caller)])
                             + "\n")

CALLER_TO_BED = {"lumpy": _vcf_to_bed,
                 "delly": _vcf_to_bed,
                 "cn_mops": _cnvbed_to_bed}

def _create_bed(call, base_file):
    """Create a simplified BED file from caller specific input.
    """
    out_file = "%s-%s.bed" % (utils.splitext_plus(base_file)[0], call["variantcaller"])
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            convert_fn = CALLER_TO_BED.get(call["variantcaller"])
            if convert_fn:
                convert_fn(call["vrn_file"], call["variantcaller"], tx_out_file)

    if utils.file_exists(out_file):
        return out_file

# ## Top level

def summarize(calls, data):
    """Summarize results from multiple callers into a single flattened BED file.
    """
    sample = tz.get_in(["rgnames", "sample"], data)
    work_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                               sample, "ensemble"))
    out_file = os.path.join(work_dir, "%s-ensemble.bed" % sample)
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            with utils.curdir_tmpdir(data) as tmpdir:
                pybedtools.set_tempdir(tmpdir)
                input_beds = filter(lambda x: x is not None,
                                    [_create_bed(c, out_file) for c in calls])
                if len(input_beds) > 0:
                    all_file = "%s-all.bed" % utils.splitext_plus(tx_out_file)[0]
                    with open(all_file, "w") as out_handle:
                        for line in fileinput.input(input_beds):
                            out_handle.write(line)
                    pybedtools.BedTool(all_file).sort(stream=True).merge(nms=True).saveas(tx_out_file)
    if utils.file_exists(out_file):
        calls.append({"variantcaller": "ensemble",
                      "vrn_file": out_file})
    return calls
