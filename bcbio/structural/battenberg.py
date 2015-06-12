"""Copy number variant calling with cgpBattenberg from Sanger

https://github.com/cancerit/cgpBattenberg
"""
import os

import toolz as tz

from bcbio import install, utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import vcfutils

def run(items, background=None):
    """Detect copy number variations from tumor/normal samples using Battenberg.
    """
    paired = vcfutils.get_paired_bams([x["align_bam"] for x in items], items)
    assert paired and paired.normal_bam, "Battenberg only works on paired tumor/normal inputs"
    assert tz.get_in(["genome_resources", "aliases", "human"], paired.tumor_data), \
        "Battenberg only works on human data"
    batout = _do_run(paired)
    batout["variantcaller"] = "battenberg"
    out = []
    for data in items:
        if "sv" not in data:
            data["sv"] = []
        data["sv"].append(batout)
        out.append(data)
    return out

def _do_run(paired):
    """Perform Battenberg caling with the paired dataset.

    This purposely does not use a temporary directory for the output
    since Battenberg does smart restarts.
    """
    work_dir = _sv_workdir(paired.tumor_data)
    out = _get_battenberg_out(paired, work_dir)
    if len(_missing_files(out)) > 0:
        ref_file = dd.get_ref_file(paired.tumor_data)
        bat_datadir = os.path.normpath(os.path.join(os.path.dirname(ref_file), os.pardir, "battenberg"))
        local_sitelib = os.path.join(install.get_defaults().get("tooldir", "/usr/local"),
                                     "lib", "R", "site-library")
        tumor_bam = paired.tumor_bam
        normal_bam = paired.normal_bam
        platform = dd.get_platform(paired.tumor_data)
        genome_build = paired.tumor_data["genome_build"]
        cores = dd.get_num_cores(paired.tumor_data)
        cmd = ("export R_LIBS_USER={local_sitelib} && "
               "battenberg.pl -t {cores} -o {work_dir} -r {ref_file}.fai "
               "-tb {tumor_bam} -nb {normal_bam} -e {bat_datadir}/impute/imput_info.txt "
               "-u {bat_datadir}/1000genomesloci -c {bat_datadir}/probloci.txt "
               "-assembly {genome_build} -species Human -platform {platform}")
        do.run(cmd.format(**locals()), "Battenberg CNV calling")
    assert len(_missing_files(out)) == 0, "Missing Battenberg output: %s" % _missing_files(out)
    return out

def _missing_files(out):
    missing_files = []
    for key, fname in out.items():
        if not os.path.exists(fname):
            missing_files.append((key, fname))
    return missing_files

def _get_battenberg_out(paired, work_dir):
    out = {}
    for key, ext in [("vrn_file", "battenberg_cn.vcf.gz"),
                     ("subclones", "subclones.txt"),
                     ("contamination", "normal_contamination.txt")]:
        out[key] = os.path.join(work_dir, "%s_%s" % (paired.tumor_name, ext))
    return out

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "battenberg"))
