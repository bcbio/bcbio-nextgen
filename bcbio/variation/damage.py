"""Annotate with potential DNA damage artifacts by examining strand/read bias.

Uses DKFZBiasFilter to identify strand and PCR bias and converts these into
INFO level annotations of low frequency variants:

https://github.com/bcbio/bcbio.github.io/blob/master/_posts/2017-01-31-damage-filters.md
"""
import os
import shutil

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import vcfutils

def run_filter(vrn_file, align_bam, ref_file, data, items):
    """Filter and annotate somatic VCFs with damage/bias artifacts on low frequency variants.
    """
    if not should_filter(items):
        return data
    else:
        out_file = "%s-damage.vcf" % utils.splitext_plus(vrn_file)[0]
        out_plot_files = ["%s%s" % (utils.splitext_plus(out_file)[0], ext)
                          for ext in ["_seq_bias_simplified.pdf", "_pcr_bias_simplified.pdf"]]
        if not utils.file_uptodate(out_file, vrn_file) and not utils.file_uptodate(out_file + ".gz", vrn_file):
            with file_transaction(items[0], out_file) as tx_out_file:
                # Does not apply --qcSummary plotting due to slow runtimes
                cmd = ["dkfzbiasfilter.py", "--filterCycles", "1", "--passOnly",
                       "--tempFolder", os.path.dirname(tx_out_file),
                       vrn_file, align_bam, ref_file, tx_out_file]
                do.run(cmd, "Filter low frequency variants for DNA damage and strand bias")
                for out_plot in out_plot_files:
                    tx_plot_file = os.path.join("%s_qcSummary" % utils.splitext_plus(tx_out_file)[0], "plots",
                                                os.path.basename(out_plot))
                    if utils.file_exists(tx_plot_file):
                        shutil.move(tx_plot_file, out_plot)
        # work around issue with vcfanno/vcfgo
        prep_cmd = ("sed 's/ACGTNacgtnPLUS,Number=10,/ACGTNacgtnPLUS,Number=.,/' |"
                    "sed 's/ACGTNacgtnMINUS,Number=10,/ACGTNacgtnMINUS,Number=.,/' ")
        out_file = vcfutils.bgzip_and_index(out_file, items[0]["config"], prep_cmd=prep_cmd)
        data["vrn_file"] = out_file
        out_plot_files = [x for x in out_plot_files if utils.file_exists(x)]
        data["damage_plots"] = out_plot_files
        return data

def should_filter(items):
    """Check if we should do damage filtering on somatic calling with low frequency events.
    """
    return (vcfutils.get_paired(items) is not None and
            any("damage_filter" in dd.get_tools_on(d) for d in items))
