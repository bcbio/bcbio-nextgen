"""Handle conversions between structural variant formats.
"""
import os

import pybedtools
import vcf

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.variation import vcfutils

# ## Conversions to simplified BED files

MAX_SVSIZE = 1e6  # 1Mb maximum size from callers to avoid huge calls collapsing all structural variants

def _vcf_to_bed(in_file, caller, out_file):
    if in_file and in_file.endswith((".vcf", "vcf.gz")):
        with utils.open_gzipsafe(in_file) as in_handle:
            with open(out_file, "w") as out_handle:
                for rec in vcf.Reader(in_handle, in_file):
                    if not rec.FILTER:
                        if (rec.samples[0].gt_type != 0 and
                              not (hasattr(rec.samples[0].data, "FT") and rec.samples[0].data.FT
                                   and rec.samples[0].data.FT not in ["PASS"])):
                            start = max(0, rec.start - 1)
                            end = rec.INFO.get("END")
                            if not end:
                                end = start + 1
                            if isinstance(end, (list, tuple)):
                                end = end[0]
                            end = int(end)
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
    with open(out_file, "w") as out_handle:
        for feat in pybedtools.BedTool(in_file):
            out_handle.write("\t".join([feat.chrom, str(feat.start), str(feat.end),
                                        "cnv%s_%s" % (feat.score, caller)])
                             + "\n")

CALLER_TO_BED = {"lumpy": _vcf_to_bed,
                 "longranger": _vcf_to_bed,
                 "delly": _vcf_to_bed,
                 "manta": _vcf_to_bed,
                 "metasv": _vcf_to_bed,
                 "cnvkit": _vcf_to_bed,
                 "gridss": _vcf_to_bed,
                 "seq2c": _vcf_to_bed,
                 "cn_mops": _cnvbed_to_bed,
                 "wham": _vcf_to_bed}
SUBSET_BY_SUPPORT = {}

def to_bed(call, sample, work_dir, calls, data):
    """Create a simplified BED file from caller specific input.
    """
    out_file = os.path.join(work_dir, "%s-%s-flat.bed" % (sample, call["variantcaller"]))
    if call.get("vrn_file") and not utils.file_uptodate(out_file, call["vrn_file"]):
        with file_transaction(data, out_file) as tx_out_file:
            convert_fn = CALLER_TO_BED.get(call["variantcaller"])
            if convert_fn:
                vrn_file = call["vrn_file"]
                if call["variantcaller"] in SUBSET_BY_SUPPORT:
                    ecalls = [x for x in calls if x["variantcaller"] in SUBSET_BY_SUPPORT[call["variantcaller"]]]
                    if len(ecalls) > 0:
                        vrn_file = _subset_by_support(call["vrn_file"], ecalls, data)
                convert_fn(vrn_file, call["variantcaller"], tx_out_file)
    if utils.file_exists(out_file):
        return out_file

def _subset_by_support(orig_vcf, cmp_calls, data):
    """Subset orig_vcf to calls also present in any of the comparison callers.
    """
    cmp_vcfs = [x["vrn_file"] for x in cmp_calls]
    out_file = "%s-inensemble.vcf.gz" % utils.splitext_plus(orig_vcf)[0]
    if not utils.file_uptodate(out_file, orig_vcf):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = "bedtools intersect -header -wa -f 0.5 -r -a {orig_vcf} -b "
            for cmp_vcf in cmp_vcfs:
                cmd += "<(bcftools view -f 'PASS,.' %s) " % cmp_vcf
            cmd += "| bgzip -c > {tx_out_file}"
            do.run(cmd.format(**locals()), "Subset calls by those present in Ensemble output")
    return vcfutils.bgzip_and_index(out_file, data["config"])
