"""Annotate with potential DNA damage artifacts by examining strand/read bias.

Uses DKFZBiasFilter to identify strand and PCR bias and converts these into
INFO level annotations of low frequency variants:

https://github.com/bcbio/bcbio.github.io/blob/master/_posts/2017-01-31-damage-filters.md
"""
import io
import os
import shutil

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.variation import vcfutils

def run_filter(vrn_file, align_bam, ref_file, data, items):
    """Filter and annotate somatic VCFs with damage/bias artifacts on low frequency variants.

    Moves damage estimation to INFO field, instead of leaving in FILTER.
    """
    if not should_filter(items) or not vcfutils.vcf_has_variants(vrn_file):
        return data
    else:
        raw_file = "%s-damage.vcf" % utils.splitext_plus(vrn_file)[0]
        out_plot_files = ["%s%s" % (utils.splitext_plus(raw_file)[0], ext)
                          for ext in ["_seq_bias_simplified.pdf", "_pcr_bias_simplified.pdf"]]
        if not utils.file_uptodate(raw_file, vrn_file) and not utils.file_uptodate(raw_file + ".gz", vrn_file):
            with file_transaction(items[0], raw_file) as tx_out_file:
                # Does not apply --qcSummary plotting due to slow runtimes
                dkfzbiasfilter = utils.which(config_utils.get_program("dkfzbiasfilter.py", data))
                cmd = [dkfzbiasfilter, "--filterCycles", "1", "--passOnly",
                       "--tempFolder", os.path.dirname(tx_out_file),
                       vrn_file, align_bam, ref_file, tx_out_file]
                do.run(cmd, "Filter low frequency variants for DNA damage and strand bias")
                for out_plot in out_plot_files:
                    tx_plot_file = os.path.join("%s_qcSummary" % utils.splitext_plus(tx_out_file)[0], "plots",
                                                os.path.basename(out_plot))
                    if utils.file_exists(tx_plot_file):
                        shutil.move(tx_plot_file, out_plot)
        raw_file = vcfutils.bgzip_and_index(raw_file, items[0]["config"])
        data["vrn_file"] = _filter_to_info(raw_file, items[0])
        out_plot_files = [x for x in out_plot_files if utils.file_exists(x)]
        data["damage_plots"] = out_plot_files
        return data

def _filter_to_info(in_file, data):
    """Move DKFZ filter information into INFO field.
    """
    header = ("""##INFO=<ID=DKFZBias,Number=.,Type=String,"""
              """Description="Bias estimation based on unequal read support from DKFZBiasFilterVariant Depth">\n""")
    out_file = "%s-ann.vcf" % utils.splitext_plus(in_file)[0]
    if not utils.file_uptodate(out_file, in_file) and not utils.file_uptodate(out_file + ".gz", in_file):
        with file_transaction(data, out_file) as tx_out_file:
            with utils.open_gzipsafe(in_file) as in_handle:
                with io.open(tx_out_file, "w", encoding="utf-8") as out_handle:
                    for line in in_handle:
                        if line.startswith("#CHROM"):
                            out_handle.write(header + line)
                        elif line.startswith("#"):
                            out_handle.write(line)
                        else:
                            out_handle.write(_rec_filter_to_info(line))
    return vcfutils.bgzip_and_index(out_file, data["config"])

def _rec_filter_to_info(line):
    """Move a DKFZBias filter to the INFO field, for a record.
    """
    parts = line.rstrip().split("\t")
    move_filters = {"bSeq": "strand", "bPcr": "damage"}
    new_filters = []
    bias_info = []
    for f in parts[6].split(";"):
        if f in move_filters:
            bias_info.append(move_filters[f])
        elif f not in ["."]:
            new_filters.append(f)
    if bias_info:
        parts[7] += ";DKFZBias=%s" % ",".join(bias_info)
    parts[6] = ";".join(new_filters or ["PASS"])
    return "\t".join(parts) + "\n"

def should_filter(items):
    """Check if we should do damage filtering on somatic calling with low frequency events.
    """
    return (vcfutils.get_paired(items) is not None and
            any("damage_filter" in dd.get_tools_on(d) for d in items))
