import os
import toolz as tz

from bcbio import bam, utils
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction

pybedtools = utils.LazyImport("pybedtools")

mapped_filter_query = (
    "not unmapped and "
    "not mate_is_unmapped and "
    "not secondary_alignment and "
    "not failed_quality_control")

def work_dir(data):
    return utils.safe_makedir(os.path.join(dd.get_work_dir(data), "coverage", dd.get_sample_name(data), "sambamba"))

def number_of_mapped_reads(data, bam_file, keep_dups=True, bed_file=None, target_name=None):
    """Count mapped reads, allow adjustment for duplicates and BED regions.

    Since samtools view does not use indexes for BED files
    (https://github.com/samtools/samtools/issues/88)
    we loop over regions in a BED file and add the counts together.
    """
    # Flag explainer https://broadinstitute.github.io/picard/explain-flags.html
    if keep_dups:
        query_name = mapped_filter_query
        flag = 780  # not (read unmapped or mate unmapped or fails QC or secondary alignment)
    else:
        query_name = mapped_filter_query + " and not duplicate"
        flag = 1804  # as above plus not duplicate
    cmd_id = "num_" + query_name.replace(" ", "_")
    if bed_file is not None:
        target_name = target_name or os.path.basename(bed_file)
        cmd_id += "_on_" + target_name
    output_file = os.path.join(work_dir(data), cmd_id)
    if not utils.file_uptodate(output_file, bam_file):
        bam.index(bam_file, data["config"], check_timestamp=False)
        num_cores = dd.get_num_cores(data)
        with file_transaction(data, output_file) as tx_out_file:
            count = 0
            # Covert to samtools regions (they are 1-based, BED is 0-based)
            regions = (["%s:%s-%s" % (r.chrom, r.start + 1, r.end) for r in pybedtools.BedTool(bed_file)]
                       if bed_file else [None])
            logger.debug("Count mapped reads with samtools view: %s" % (dd.get_sample_name(data)))
            for i, region_group in enumerate(tz.partition_all(10000, regions)):
                if len(region_group) == 1 and not region_group[0]:
                    region_str = ""
                else:
                    region_in = "%s-regions-%s.bed" % (utils.splitext_plus(tx_out_file)[0], i)
                    with open(region_in, "w") as out_handle:
                        out_handle.write(" ".join(region_group))
                    region_str = " `cat %s`" % region_in
                count_out = "%s-count-%s.txt" % (utils.splitext_plus(tx_out_file)[0], i)
                cmd = "samtools view -c -F {flag} -@ {num_cores} {bam_file}{region_str} > {count_out}"
                do.run(cmd.format(**locals()))
                with open(count_out) as in_handle:
                    count += int(in_handle.read().strip())
            with open(tx_out_file, "w") as out_handle:
                out_handle.write(str(count))
    with open(output_file) as f:
        return int(f.read().strip())
