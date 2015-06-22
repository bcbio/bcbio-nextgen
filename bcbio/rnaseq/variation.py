import os
from bcbio.utils import file_exists
import bcbio.pipeline.datadict as dd
from bcbio.ngsalign.postalign import dedup_bam
from bcbio.distributed.transaction import file_transaction
from bcbio import broad, bam


def rnaseq_gatk_variant_calling(data):
    data = dd.set_deduped_bam(data, dedup_bam(dd.get_work_bam(data), data))
    data = gatk_splitreads(data)
    data = gatk_rnaseq_calling(data)
    return data

def gatk_splitreads(data):
    """
    use GATK to split reads with Ns in the CIGAR string, hard clipping regions
    that end up in introns
    """
    broad_runner = broad.runner_from_config(dd.get_config(data))
    ref_file = dd.get_ref_file(data)
    deduped_bam = dd.get_deduped_bam(data)
    base, ext = os.path.splitext(deduped_bam)
    split_bam = base + ".splitN" + ext
    if dd.get_quality_format(data) == "illumina":
        quality_flag = ["--fix_misencoded_quality_scores", "-fixMisencodedQuals"]
    else:
        quality_flag = []
    if file_exists(split_bam):
        data = dd.set_split_bam(data, split_bam)
        return data
    with file_transaction(split_bam) as tx_split_bam:
        params = ["-T", "SplitNCigarReads",
                  "-R", ref_file,
                  "-I", deduped_bam,
                  "-o", tx_split_bam,
                  "-rf", "ReassignOneMappingQuality",
                  "-RMQF", "255",
                  "-RMQT", "60",
                  "-rf", "UnmappedRead",
                  "-U", "ALLOW_N_CIGAR_READS"] + quality_flag
        broad_runner.run_gatk(params)
    bam.index(split_bam, dd.get_config(data))
    data = dd.set_split_bam(data, split_bam)
    return data

def gatk_rnaseq_calling(data):
    """
    use GATK to perform variant calling on RNA-seq data
    """
    broad_runner = broad.runner_from_config(dd.get_config(data))
    ref_file = dd.get_ref_file(data)
    split_bam = dd.get_split_bam(data)
    out_file = os.path.splitext(split_bam)[0] + ".gvcf"
    num_cores = dd.get_num_cores(data)
    if file_exists(out_file):
        data = dd.set_vrn_file(data, out_file)
        return data
    with file_transaction(out_file) as tx_out_file:
        params = ["-T", "HaplotypeCaller",
                  "-R", ref_file,
                  "-I", split_bam,
                  "-o", tx_out_file,
                  "-nct", str(num_cores),
                  "--emitRefConfidence", "GVCF",
                  "--variant_index_type", "LINEAR",
                  "--variant_index_parameter", "128000",
                  "-dontUseSoftClippedBases",
                  "-stand_call_conf", "20.0",
                  "-stand_emit_conf", "20.0"]
        broad_runner.run_gatk(params)
    data = dd.set_vrn_file(data, out_file)
    return data

def gatk_joint_calling(data, vrn_files, ref_file, out_file=None):
    if out_file is None:
        out_file = os.path.join("variation", "combined.vcf")
    if not file_exists(out_file):
        out_file = _run_genotype_gvcfs(data, vrn_files, ref_file, out_file)
    return out_file

def _run_genotype_gvcfs(data, vrn_files, ref_file, out_file):
    if not file_exists(out_file):
        broad_runner = broad.runner_from_config(data["config"])
        with file_transaction(data, out_file) as tx_out_file:
            params = ["-T", "GenotypeGVCFs",
                      "-R", ref_file, "-o", tx_out_file]
            for vrn_file in vrn_files:
                params += ["--variant", vrn_file]
            broad_runner.new_resources("gatk-haplotype")
            cores = dd.get_cores(data)
            if cores > 1:
                params += ["-nt", str(cores)]
                memscale = {"magnitude": 0.9 * cores, "direction": "increase"}
            else:
                memscale = None
            broad_runner.run_gatk(params, memscale=memscale)
    return out_file
