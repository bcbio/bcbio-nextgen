"""calling using Pindel

http://gmt.genome.wustl.edu/packages/pindel/
"""

from __future__ import print_function
import os

try:
    import vcf
except ImportError:
    vcf = None

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions, remove_lcr_regions
from bcbio.provenance import do
#from bcbio.variation import annotation
#from bcbio.variation.vcfutils import get_paired_bams, is_paired_analysis, bgzip_and_index


def _pindel_options_from_config(items, config, out_file, region, tmp_path):
    variant_regions = utils.get_in(config, ("algorithm", "variant_regions"))
    target = subset_variant_regions(variant_regions, region, out_file, items)
    opts = ""
    if target:
        if isinstance(target, basestring) and os.path.isfile(target):
            target_bed = target
        else:
            target_bed = os.path.join(tmp_path, "tmp.bed")
            with file_transaction(config, target_bed) as tx_tmp_bed:
                if not isinstance(region, (list, tuple)):
                    message = ("Region must be a tuple - something odd just happened")
                    raise ValueError(message)
                chrom, start, end = region
                with open(tx_tmp_bed, "w") as out_handle:
                    print("%s\t%s\t%s" % (chrom, start, end), file=out_handle)
        opts = "-j" + remove_lcr_regions(target_bed, items)
    return opts


def _pindel_regions(bed_file):
    """get regions from bed file"""
    regions = set()
    with open(bed_file) as in_handle:
        for line in in_handle:
            cols = line.strip().split("\t")
            regions.add((cols[0], int(cols[1]), int(cols[2])))
    return regions


def is_installed(config):
    """Check for pindel installation on machine.
    """
    try:
        config_utils.get_program("pindel", config)
        return True
    except config_utils.CmdNotFound:
        return False


def run(items):
    """in case go inside structural variation"""
    if not all(utils.get_in(data, ("config", "algorithm", "aligner")) == "bwa" for data in items):
        raise ValueError("Require bwa-mem alignment input for lumpy structural variation detection")
    work_dir = utils.safe_makedir(os.path.join(items[0]["dirs"]["work"], "structural", items[0]["name"][-1],
                                               "pindel"))
    work_bams = [x["align_bam"] for x in items]
    ref_file = items[0]["sam_ref"]
    out_file = os.path.join(work_dir, items[0]["name"][-1])
    out_file = _run_pindel_caller(work_bams, items, ref_file, out_file=out_file)
    return items


def run_pindel(align_bams, items, ref_file, assoc_files, region=None,
               out_file):
    """Run pindel calling, either paired tumor/normal or germline calling.
    """
    config = items[0]["config"]
    out_file = run_pindel_caller(align_bams, items, ref_file, out_file)
    return merge_call_file


def _run_pindel_caller(align_bams, items, ref_file,
                       out_file=None):
    """Detect SV with pindel.

    Single sample mode.
    """
    print("parameters{align_bams} {ref_file} {region}".format(**locals()))
    config = items[0]["config"]
    if out_file is None:
        out_file = "%s-pindelroot" % os.path.splitext(align_bams[0])[0]
    if not utils.file_exists(out_file + "_SI"):
        with file_transaction(config, out_file) as tx_out_file:
            for align_bam in align_bams:
                bam.index(align_bam, config)
            pindel = config_utils.get_program("pindel", config)
            tmp_path = os.path.dirname(tx_out_file)
            opts = _pindel_options_from_config(items, config, out_file, region, tmp_path)
            sample_name_str = items[0]["name"][1]
            tmp_input = _create_tmp_input(align_bams, tmp_path, config)
            cmd = ("{pindel} -f {ref_file} -i {tmp_input} -o {out_file} {opts} -r false -l false -k false -t false ")
            print("final command for pindel %s" % cmd.format(**locals()))
            do.run(cmd.format(**locals()), "Genotyping with pindel", {})
    out_file = _create_vcf(out_file)
    return out_file


def _create_tmp_input(input_bams, tmp_path, config):
    """Create input file for pindel. tab file: bam file, insert size, name"""
    tmp_input = os.path.join(tmp_path, "pindel.txt")
    with open(tmp_input, 'w') as out_handle:
        for bam_file in input_bams:
            name_bam = utils.splitext_plus(os.path.basename(bam_file))[0]
            print("%s\t%s\t%s\n" % (bam_file, 250, name_bam), file=out_handle)
            print("%s\t%s\t%s\n" % (bam_file, 250, name_bam))
    return tmp_input


def _create_vcf(root_file, reference, name_ref, data_ref):
    """use pindel2vcf to create vcf format"""
    in_file = root_file + "_SI"
    out_file = root_file + "-SI-pindel.vcf"
    if not out_file:
        pindel2vcf = config_utils.get_program("pindel2vcf", config)
        cmd = ("{pindel2vcf} -p {in_file} -r {reference} -R {name_ref} -d {date_ref} -v {out_file}")
        print("final command for pindel %s" % cmd.format(**locals()))
        do.run(cmd.format(**locals()), "Converting to vcf", {})
    return out_file


def _merge_output(prefix, data):
    """merge outputs files as one single bed file"""
    exts = ["D", "INV", "LI", "SI", "TD"]
    out_file = prefix + "-pindel.bed"
    for ext in exts:
        single_file = prefix + "_" + ext
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, 'w') as out_handle:
                for ext in exts:
                    single_file = prefix + "_" + ext
                    with open(single_file) as in_handle:
                        line = in_handle.readline()
                        while line:
                            if line.startswith('#'):
                                line = in_handle.readline()
                                cols = line.strip().split("\t")
                                out_handle.write("%s\t%s\t%s\t%s\n" % (cols[3].split(" ")[1], cols[4].split(" ")[1], cols[5], cols[1]))
                            line = in_handle.readline()
    return out_file

