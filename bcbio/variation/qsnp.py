"""identify somatic variants in cancer samples

https://qcmg.org/bioinformatics/tiki-index.php?page=qSNP#EXAMPLES
"""

from __future__ import print_function
import os
import shutil
from re import sub

from bcbio import utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.provenance import do
from bcbio.variation import annotation, bedutils
from bcbio.variation.vcfutils import get_paired_bams, bgzip_and_index, combine_variant_files, PairedData

import six


def is_installed(config):
    """Check for qsnp installation on machine.
    """
    try:
        config_utils.get_program("qsnp", config)
        return True
    except config_utils.CmdNotFound:
        return False

def run_qsnp(align_bams, items, ref_file, assoc_files, region=None,
             out_file=None):
    """Run qSNP calling on paired tumor/normal.
    """
    if utils.file_exists(out_file):
        return out_file
    paired = get_paired_bams(align_bams, items)
    if paired.normal_bam:
        region_files = []
        regions = _clean_regions(items, region)
        if regions:
            for region in regions:
                out_region_file = out_file.replace(".vcf.gz", _to_str(region) + ".vcf.gz")
                region_file = _run_qsnp_paired(align_bams, items, ref_file,
                                               assoc_files, region, out_region_file)
                region_files.append(region_file)
            out_file = combine_variant_files(region_files, out_file, ref_file, items[0]["config"])
        if not region:
            out_file = _run_qsnp_paired(align_bams, items, ref_file,
                                        assoc_files, region, out_file)
        return out_file
    else:
        raise ValueError("qSNP only works on paired samples")

def _run_qsnp_paired(align_bams, items, ref_file, assoc_files,
                     region=None, out_file=None):
    """Detect somatic mutations with qSNP.

    This is used for paired tumor / normal samples.
    """
    config = items[0]["config"]
    if out_file is None:
        out_file = "%s-paired-variants.vcf" % os.path.splitext(align_bams[0])[0]
    if not utils.file_exists(out_file):
        out_file = out_file.replace(".gz", "")
        with file_transaction(config, out_file) as tx_out_file:
            with tx_tmpdir(config) as tmpdir:
                with utils.chdir(tmpdir):
                    paired = get_paired_bams(align_bams, items)
                    qsnp = config_utils.get_program("qsnp", config)
                    resources = config_utils.get_resources("qsnp", config)
                    mem = " ".join(resources.get("jvm_opts", ["-Xms750m -Xmx4g"]))
                    qsnp_log = os.path.join(tmpdir, "qsnp.log")
                    qsnp_init = os.path.join(tmpdir, "qsnp.ini")
                    if region:
                        paired = _create_bam_region(paired, region, tmpdir)
                    _create_input(paired, tx_out_file, ref_file, assoc_files['dbsnp'], qsnp_init)
                    cl = ("{qsnp} {mem} -i {qsnp_init} -log {qsnp_log}")
                    do.run(cl.format(**locals()), "Genotyping paired variants with Qsnp", {})
        out_file = _filter_vcf(out_file)
        out_file = bgzip_and_index(out_file, config)
    return out_file

def _clean_regions(items, region):
    """Intersect region with target file if it exists"""
    variant_regions = bedutils.population_variant_regions(items, merged=True)
    with utils.tmpfile() as tx_out_file:
        target = subset_variant_regions(variant_regions, region, tx_out_file, items)
        if target:
            if isinstance(target, six.string_types) and os.path.isfile(target):
                target = _load_regions(target)
            else:
                target = [target]
            return target

def _load_regions(target):
    """Get list of tupples from bed file"""
    regions = []
    with open(target) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                c, s, e = line.strip().split("\t")
                regions.append((c, s, e))
    return regions

def _create_bam_region(paired, region, tmp_dir):
    """create temporal normal/tumor bam_file only with reads on that region"""
    tumor_name, normal_name = paired.tumor_name, paired.normal_name
    normal_bam = _slice_bam(paired.normal_bam, region, tmp_dir, paired.tumor_config)
    tumor_bam = _slice_bam(paired.tumor_bam, region, tmp_dir, paired.tumor_config)
    paired = PairedData(tumor_bam, tumor_name, normal_bam, normal_name, None, None, None)
    return paired

def _slice_bam(in_bam, region, tmp_dir, config):
    """Use sambamba to slice a bam region"""
    name_file = os.path.splitext(os.path.basename(in_bam))[0]
    out_file = os.path.join(tmp_dir, os.path.join(tmp_dir, name_file + _to_str(region) + ".bam"))
    sambamba = config_utils.get_program("sambamba", config)
    region = _to_sambamba(region)
    with file_transaction(out_file) as tx_out_file:
        cmd = ("{sambamba} slice {in_bam} {region} -o {tx_out_file}")
        do.run(cmd.format(**locals()), "Slice region", {})
    return out_file

def _create_input(paired, out_file, ref_file, snp_file, qsnp_file):
    """Create INI input for qSNP"""
    ini_file["[inputFiles]"]["dbSNP"] = snp_file
    ini_file["[inputFiles]"]["ref"] = ref_file
    ini_file["[inputFiles]"]["normalBam"] = paired.normal_bam
    ini_file["[inputFiles]"]["tumourBam"] = paired.tumor_bam
    ini_file["[ids]"]["normalSample"] = paired.normal_name
    ini_file["[ids]"]["tumourSample"] = paired.tumor_name
    ini_file["[ids]"]["donor"] = paired.tumor_name
    ini_file["[outputFiles]"]["vcf"] = out_file
    with open(qsnp_file, "w") as out_handle:
        for k, v in ini_file.items():
            out_handle.write("%s\n" % k)
            for opt, value in v.items():
                if value != "":
                    out_handle.write("%s = %s\n" % (opt, value))

def _has_ambiguous_ref_allele(line):
    if not line.startswith("#"):
        parts = line.split("\t")
        return len(parts) > 4 and parts[3].upper() not in set(["C", "A", "T", "G"])

def _filter_vcf(out_file):
    """Fix sample names, FILTER and FORMAT fields. Remove lines with ambiguous reference.
    """
    in_file = out_file.replace(".vcf", "-ori.vcf")
    FILTER_line = ('##FILTER=<ID=SBIAS,Description="Due to bias">\n'
                   '##FILTER=<ID=5BP,Description="Due to 5BP">\n'
                   '##FILTER=<ID=REJECT,Description="Not somatic due to qSNP filters">\n')
    SOMATIC_line = '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="somatic event">\n'
    if not utils.file_exists(in_file):
        shutil.move(out_file, in_file)
    with file_transaction(out_file) as tx_out_file:
        with open(in_file) as in_handle, open(tx_out_file, "w") as out_handle:
            for line in in_handle:
                if line.startswith("##normalSample="):
                    normal_name = line.strip().split("=")[1]
                if line.startswith("##patient_id="):
                    tumor_name = line.strip().split("=")[1]
                if line.startswith("#CHROM"):
                    line = line.replace("Normal", normal_name)
                    line = line.replace("Tumour", tumor_name)
                if line.startswith("##INFO=<ID=FS"):
                    line = line.replace("ID=FS", "ID=RNT")
                if line.find("FS=") > -1:
                    line = line.replace("FS=", "RNT=")
                if "5BP" in line:
                    line = sub("5BP[0-9]+", "5BP", line)
                if line.find("PASS") == -1:
                    line = _set_reject(line)
                if line.find("PASS") > - 1 and line.find("SOMATIC") == -1:
                    line = _set_reject(line)
                if not _has_ambiguous_ref_allele(line):
                    out_handle.write(line)
                if line.startswith("##FILTER") and FILTER_line:
                    out_handle.write("%s" % FILTER_line)
                    FILTER_line = ""
                if line.startswith("##INFO") and SOMATIC_line:
                    out_handle.write("%s" % SOMATIC_line)
                    SOMATIC_line = ""
    return out_file

def _to_str(region):
    return "_" + "_".join(map(str, list(region)))

def _to_sambamba(region):
    return "%s:%s-%s" % (region[0], region[1]+1, region[2]+1)

def _set_reject(line):
    """Set REJECT in VCF line, or add it if there is something else."""
    if line.startswith("#"):
        return line
    parts = line.split("\t")
    if parts[6] == "PASS":
        parts[6] = "REJECT"
    else:
        parts[6] += ";REJECT"
    return "\t".join(parts)


ini_file = {"[inputFiles]":{
          "dbSNP":"",
          "ref":"",
          "normalBam":"",
          "tumourBam":""},

          "[parameters]":{
          "annotateMode":"vcf",
          "runMode":"standard",
          "minimumBaseQuality":"10",
          "includeIndels": "true"},

          "[ids]":{
          "donor":"",
          "normalSample":"",
          "tumourSample":"",
          "somaticAnalysis":"",
          "germlineAnalysis":""},

          "[outputFiles]":{
          "vcf":""},

          "[rules]":{
          "normal1":"0,20,3",
          "normal2":"21,50,4",
          "normal3":"51,,10",
          "tumour1":"0,20,3",
          "tumour2":"21,50,4",
          "tumour3":"51,,5"}}
