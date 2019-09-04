import os
import shutil
import bcbio.bam as bam
from bcbio.utils import (file_exists, safe_makedir, append_stem)
from bcbio.pipeline import config_utils
from bcbio.bam import is_paired
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
import bcbio.pipeline.datadict as dd

try:
    import pandas as pd
except ImportError:
    pd = None

def count(data):
    """
    count reads mapping to genes using featureCounts
    http://subread.sourceforge.net
    """
    in_bam = dd.get_work_bam(data) or dd.get_align_bam(data)
    out_dir = os.path.join(dd.get_work_dir(data), "align", dd.get_sample_name(data))
    if dd.get_aligner(data) == "star":
        out_dir = os.path.join(out_dir, "%s_%s" % (dd.get_sample_name(data), dd.get_aligner(data)))
    sorted_bam = bam.sort(in_bam, dd.get_config(data), order="queryname", out_dir=safe_makedir(out_dir))
    gtf_file = dd.get_gtf_file(data)
    work_dir = dd.get_work_dir(data)
    out_dir = os.path.join(work_dir, "htseq-count")
    safe_makedir(out_dir)
    count_file = os.path.join(out_dir, dd.get_sample_name(data)) + ".counts"
    summary_file = os.path.join(out_dir, dd.get_sample_name(data)) + ".counts.summary"
    if file_exists(count_file) and _is_fixed_count_file(count_file):
        return count_file

    featureCounts = config_utils.get_program("featureCounts", dd.get_config(data))
    paired_flag = _paired_flag(in_bam)
    strand_flag = _strand_flag(data)

    filtered_bam = bam.filter_primary(sorted_bam, data)

    cmd = ("{featureCounts} -a {gtf_file} -o {tx_count_file} -s {strand_flag} "
           "{paired_flag} {filtered_bam}")

    message = ("Count reads in {tx_count_file} mapping to {gtf_file} using "
               "featureCounts")
    with file_transaction(data, [count_file, summary_file]) as tx_files:
        tx_count_file, tx_summary_file = tx_files
        do.run(cmd.format(**locals()), message.format(**locals()))
    fixed_count_file = _format_count_file(count_file, data)
    fixed_summary_file = _change_sample_name(
        summary_file, dd.get_sample_name(data), data=data)
    shutil.move(fixed_count_file, count_file)
    shutil.move(fixed_summary_file, summary_file)

    return count_file

def _change_sample_name(in_file, sample_name, data=None):
    """Fix name in feature counts log file to get the same
       name in multiqc report.
    """
    out_file = append_stem(in_file, "_fixed")
    with file_transaction(data, out_file) as tx_out:
        with open(tx_out, "w") as out_handle:
            with open(in_file) as in_handle:
                for line in in_handle:
                    if line.startswith("Status"):
                        line = "Status\t%s.bam" % sample_name
                    out_handle.write("%s\n" % line.strip())
    return out_file

def _is_fixed_count_file(count_file):
    if os.path.exists(count_file):
        with open(count_file) as in_handle:
            line = in_handle.readline().split("\t")
            return len(line) == 2

def _format_count_file(count_file, data):
    """
    this cuts the count file produced from featureCounts down to
    a two column file of gene ids and number of reads mapping to
    each gene
    """
    COUNT_COLUMN = 5
    out_file = os.path.splitext(count_file)[0] + ".fixed.counts"
    if file_exists(out_file) and _is_fixed_count_file(out_file):
        return out_file

    df = pd.io.parsers.read_csv(count_file, sep="\t", index_col=0, header=1)
    df_sub = df.iloc[:, COUNT_COLUMN]
    with file_transaction(data, out_file) as tx_out_file:
        df_sub.to_csv(tx_out_file, sep="\t", index_label="id", header=False)
    return out_file

def _strand_flag(data):
    """
    0: unstranded 1: stranded 2: reverse stranded
    """
    strand_flag = {"unstranded": "0",
                   "firststrand": "2",
                   "secondstrand": "1"}
    stranded = dd.get_strandedness(data)

    assert stranded in strand_flag, ("%s is not a valid strandedness value. "
                                     "Valid values are 'firststrand', 'secondstrand', "
                                     "and 'unstranded")
    return strand_flag[stranded]

def _paired_flag(bam_file):
    """
    sets flags to handle paired-end BAM files
    """
    if is_paired(bam_file):
        return "-p -B -C"
    else:
        return ""
