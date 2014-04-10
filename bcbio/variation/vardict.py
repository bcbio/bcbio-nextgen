"""Sensitive variant calling using VarDict.

"""
import os
import itertools

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.provenance import do
from bcbio.variation import annotation
from bcbio.variation.vcfutils import get_paired_bams, is_paired_analysis, merge_variant_files
from bcbio.variation import bamprep


def _vardict_options_from_config(items, config, out_file, region=None):
    opts = []

    resources = config_utils.get_resources("vardict", config)
    if resources.get("options"):
        opts += resources["options"]

    variant_regions = utils.get_in(config, ("algorithm", "variant_regions"))
    target = subset_variant_regions(variant_regions, region, out_file)
    if target:
        if isinstance(target, basestring) and os.path.isfile(target):
            opts += [target] # this must be the last option
        else:
            # one-based, end-inclusive coordinates as for Gatk
            opts += ["-R", bamprep.region_to_gatk(target)]
    return opts

def run_vardict(align_bams, items, ref_file, assoc_files, region=None,
                  out_file=None):
    """Run VarDict variant calling.
    """
    call_file = _run_vardict_caller(align_bams, items, ref_file,
                                          assoc_files, region, out_file)
    return call_file

def _run_vardict_caller(align_bams, items, ref_file, assoc_files,
                          region=None, out_file=None):
    """Detect SNPs and indels with VarDict.
    """
    config = items[0]["config"]
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % os.path.splitext(align_bams[0])[0]
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            for align_bam in align_bams:
                bam.index(align_bam, config)
            
            #vardict = config_utils.get_program("vardict", config)
            #input_bams = " ".join("-b %s" % x for x in align_bams)
            num_bams = len(align_bams)
            sample_vcf_names = [] # for individual sample names, given batch calling may be required
            for bamfile, item in itertools.izip(align_bams, items):
                temp_file_prefix = out_file.replace(".gz","").replace(".vcf","") + item["name"][1]
                tmp_out1 = temp_file_prefix + ".temp1.txt"
                tmp_out2 = temp_file_prefix + ".temp2.txt"
                tmp_out3 = temp_file_prefix + ".temp3.vcf"
                if out_file.endswith("gz"):
                    tmp_out3 += ".gz"
                # prepare commands
                vardict_dir = os.path.dirname(os.path.realpath(config_utils.get_program("vardict", config)))
                # assumption: the following scripts are all in the same path
                vardict = os.path.join(vardict_dir, "vardict")
                if not utils.file_exists(vardict):
                    vardict = "vardict.pl" # desperado attempt
                strandbias = os.path.join(vardict_dir, "teststrandbias.R")
                if not utils.file_exists(strandbias):
                    strandbias = "teststrandbias.R" # desperado attempt
                var2vcf = os.path.join(vardict_dir, "var2vcf_valid.pl")
                if not utils.file_exists(var2vcf):
                    var2vcf = "var2vcf_valid.pl" # desperado attempt
                sample_vcf_names.append(tmp_out3)
                opts = " ".join(_vardict_options_from_config(items, config, out_file, region))            
                vcfstreamsort = config_utils.get_program("vcfstreamsort", config)
                compress_cmd = "| bgzip -c" if out_file.endswith("gz") else ""            
                freq = 0.01
                sample = item["name"][1]
                # 3 steps to produce vcf
                # 1. scan and sift BAM
                with file_transaction(tmp_out1) as tx_tmp1_file:
                    cmd = ("{vardict} -z -c 1 -S 2 -E 3 -g 4 -x 0 -G {ref_file} -f {freq} -k 3 " 
                           "-X 5 -N {sample} -b {bamfile} {opts} > {tx_tmp1_file}")
                    do.run(cmd.format(**locals()), "Genotyping with VarDict: Inference", {})
                # 2. filter based on strand bias
                with file_transaction(tmp_out2) as tx_tmp2_file:
                    cmd = ("{strandbias} {tmp_out1} > {tx_tmp2_file}")
                    do.run(cmd.format(**locals()), "Genotyping with VarDict: Strand bias", {})
                # 3. produce vcf
                if num_bams > 1:
                    with file_transaction(tmp_out3) as tx_tmp3_file:
                        cmd = ("{var2vcf} -S -f {freq} {tmp_out2} | {vcfstreamsort} {compress_cmd} > {tmp_out3}")
                        do.run(cmd.format(**locals()), "Genotyping with VarDict: VCF output", {})
                else:
                    cmd = ("var2vcf_valid.pl -S -f {freq} -N {sample} {tmp_out2} | {vcfstreamsort} {compress_cmd} > {tx_out_file}")
                    do.run(cmd.format(**locals()), "Genotyping with VarDict: VCF output", {})
            if num_bams > 1:
                # N.B. merge_variant_files wants region in 1-based end-inclusive 
                # coordinates. Thus use bamprep.region_to_gatk
                merge_variant_files(orig_files=sample_vcf_names, 
                                    out_file=tx_out_file, ref_file=ref_file, 
                                    config=config, region=bamprep.region_to_gatk(region))            
    ann_file = annotation.annotate_nongatk_vcf(out_file, align_bams,
                                               assoc_files["dbsnp"],
                                               ref_file, config)
    return ann_file

