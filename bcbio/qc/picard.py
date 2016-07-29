import os

from bcbio import utils
from bcbio import broad
from bcbio.broad.metrics import PicardMetrics
from bcbio import bam
from bcbio.distributed.transaction import tx_tmpdir
from bcbio.provenance import do
from bcbio.pipeline import datadict as dd

def run(bam_file, data, out_dir):
    config = data["config"]
    if "picard" not in dd.get_tools_on(data):
        return {}
    ref_file = dd.get_ref_file(data)
    sample = dd.get_sample_name(data)
    target_file = dd.get_variant_regions(data)
    broad_runner = broad.PicardCmdRunner("picard", config)
    bam_fname = os.path.abspath(bam_file)
    path = os.path.dirname(bam_fname)
    utils.safe_makedir(out_dir)
    hsmetric_file = os.path.join(out_dir, "%s-sort.hs_metrics" % sample)
    hsinsert_file = os.path.join(out_dir, "%s-sort.insert_metrics" % sample)
    if utils.file_exists(hsmetric_file):
        return hsmetric_file
    with utils.chdir(out_dir):
        with tx_tmpdir() as tmp_dir:
            cur_bam = os.path.basename(bam_fname)
            if not os.path.exists(cur_bam):
                os.symlink(bam_fname, cur_bam)
            gen_metrics = PicardMetrics(broad_runner, tmp_dir)
            gen_metrics.report(cur_bam, ref_file,
                               bam.is_paired(bam_fname),
                               target_file, target_file, None, config)
    do.run("sed -i 's/-sort.bam//g' %s" % hsmetric_file, "")
    do.run("sed -i 's/-sort.bam//g' %s" % hsinsert_file, "")
    return hsmetric_file


