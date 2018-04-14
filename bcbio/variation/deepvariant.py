"""DeepVariant calling: https://github.com/google/deepvariant
"""
import os
import glob

from bcbio import utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import strelka2, vcfutils

def run(align_bams, items, ref_file, assoc_files, region, out_file):
    """Return DeepVariant calling on germline samples.

    region can be a single region or list of multiple regions for multicore calling.
    """
    assert not vcfutils.is_paired_analysis(align_bams, items), \
        ("DeepVariant currently only supports germline calling: %s" %
         (", ".join([dd.get_sample_name(d) for d in items])))
    assert len(items) == 1, \
        ("DeepVariant currently only supports single sample calling: %s" %
         (", ".join([dd.get_sample_name(d) for d in items])))
    out_file = _run_germline(align_bams[0], items[0], ref_file,
                             region, out_file)
    return vcfutils.bgzip_and_index(out_file, items[0]["config"])

def _run_germline(bam_file, data, ref_file, region, out_file):
    """Single sample germline variant calling.
    """
    work_dir = utils.safe_makedir("%s-work" % utils.splitext_plus(out_file)[0])
    region_bed = strelka2.get_region_bed(region, [data], out_file, want_gzip=False)
    example_dir = _make_examples(bam_file, data, ref_file, region_bed, out_file, work_dir)
    if _has_candidate_variants(example_dir):
        tfrecord_file = _call_variants(example_dir, region_bed, data, out_file)
        return _postprocess_variants(tfrecord_file, data, ref_file, out_file)
    else:
        return vcfutils.write_empty_vcf(out_file, data["config"], [dd.get_sample_name(data)])

def _has_candidate_variants(example_dir):
    return all(utils.is_empty_gzipsafe(f) for f in glob.glob(os.path.join(example_dir, "*tfrecord*gz")))

def _make_examples(bam_file, data, ref_file, region_bed, out_file, work_dir):
    """Create example pileup images to feed into variant calling.
    """
    log_dir = utils.safe_makedir(os.path.join(work_dir, "log"))
    example_dir = utils.safe_makedir(os.path.join(work_dir, "examples"))
    if len(glob.glob(os.path.join(example_dir, "%s.tfrecord*.gz" % dd.get_sample_name(data)))) == 0:
        with tx_tmpdir(data) as tx_example_dir:
            cmd = ["dv_make_examples.py", "--cores", dd.get_num_cores(data), "--ref", ref_file,
                   "--reads", bam_file, "--regions", region_bed, "--logdir", log_dir,
                   "--examples", tx_example_dir, "--sample", dd.get_sample_name(data)]
            do.run(cmd, "DeepVariant make_examples %s" % dd.get_sample_name(data))
            for fname in glob.glob(os.path.join(tx_example_dir, "%s.tfrecord*.gz" % dd.get_sample_name(data))):
                utils.copy_plus(fname, os.path.join(example_dir, os.path.basename(fname)))
    return example_dir

def _call_variants(example_dir, region_bed, data, out_file):
    """Call variants from prepared pileup examples, creating tensorflow record file.
    """
    tf_out_file = "%s-tfrecord.gz" % utils.splitext_plus(out_file)[0]
    if not utils.file_exists(tf_out_file):
        with file_transaction(data, tf_out_file) as tx_out_file:
            model = "wes" if strelka2.coverage_interval_from_bed(region_bed) == "targeted" else "wgs"
            cmd = ["dv_call_variants.py", "--cores", dd.get_num_cores(data),
                   "--outfile", tx_out_file, "--examples", example_dir,
                   "--sample", dd.get_sample_name(data), "--model", model]
            do.run(cmd, "DeepVariant call_variants %s" % dd.get_sample_name(data))
    return tf_out_file

def _postprocess_variants(record_file, data, ref_file, out_file):
    """Post-process variants, converting into standard VCF file.
    """
    if not utils.file_uptodate(out_file, record_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = ["dv_postprocess_variants.py", "--ref", ref_file,
                   "--infile", record_file, "--outfile", tx_out_file]
            do.run(cmd, "DeepVariant postprocess_variants %s" % dd.get_sample_name(data))
    return out_file
