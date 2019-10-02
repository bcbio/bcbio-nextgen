"""
functions to access the data dictionary in a clearer way
"""

import os
import toolz as tz
from bcbio.utils import file_exists, to_single_data
from bcbio.log import logger
import sys

LOOKUPS = {
    "config": {"keys": ['config']},
    "tmp_dir": {"keys": ['config', 'resources', 'tmp', 'dir']},
    "num_cores": {"keys": ['config', 'algorithm', 'num_cores'],
                  "default": 1},
    "svprioritize": {"keys": ['config', 'algorithm', 'svprioritize']},
    "effects_transcripts": {"keys": ["config", "algorithm", "effects_transcripts"], "default": "all"},
    "genome_build": {"keys": ["genome_build"]},
    "gtf_file": {"keys": ['genome_resources', 'rnaseq', 'transcripts'],
                 "checker": file_exists},
    "transcriptome_fasta": {"keys": ["config", "algorithm", "transcriptome_fasta"],
                            "default": None},
    "transcriptome_gtf": {"keys": ["config", "algorithm", "transcriptome_gtf"],
                            "default": None},
    "singlecell_quantifier": {"keys": ["config", "algorithm",
                                       "singlecell_quantifier"],
                            "default": "rapmap"},
    "positional_umi": {"keys": ["config", "algorithm", "positional_umi"]},
    "tx2gene": {"keys": ["tx2gene"]},
    "ref_file": {"keys": ["reference", "fasta", "base"]},
    "ref_file_compressed": {"keys": ["reference", "fastagz", "base"]},
    "srna_gtf_file": {"keys": ['genome_resources', 'srnaseq', 'srna_transcripts'],
                      "checker": file_exists},
    "srna_trna_file": {"keys": ['genome_resources', 'srnaseq', 'trna_fasta'],
                       "checker": file_exists},
    "srna_mint_lookup": {"keys": ['genome_resources', 'srnaseq', 'mint_lookup'],
                         "checker": file_exists},
    "srna_mint_space": {"keys": ['genome_resources', 'srnaseq', 'mint_space'],
                        "checker": file_exists},
    "srna_mint_other": {"keys": ['genome_resources', 'srnaseq', 'mint_other'],
                        "checker": file_exists},
    "mirdeep2_file": {"keys": ['genome_resources', 'srnaseq', 'mirdeep2_fasta'],
                      "checker": file_exists},
    "mirbase_hairpin": {"keys": ['genome_resources', 'srnaseq', 'mirbase_hairpin'],
                        "checker": file_exists},
    "mirbase_mature": {"keys": ['genome_resources', 'srnaseq', 'mirbase_mature'],
                       "checker": file_exists},
    "gene_bed": {"keys": ['genome_resources', 'rnaseq', 'gene_bed'],
                 "checker": file_exists},
    "work_dir": {"keys": ['dirs', 'work']},
    "sam_ref": {"keys": ["sam_ref"]},
    "disambiguate": {"keys": ["config", "algorithm", "disambiguate"],
                     "default": False},
    "lane": {"keys": ["rgnames", "lane"]},
    "cores": {"keys": ["config", "algorithm", "num_cores"], "default": 1},
    "sample_name": {"keys": ['rgnames', 'sample']},
    "strandedness": {"keys": ['config', 'algorithm', 'strandedness'],
                     "default": "unstranded"},
    "vcfanno": {"keys": ['config', 'algorithm', 'vcfanno'], "default": [], "always_list": True},
    "analysis": {"keys": ["analysis"], "default": ""},
    "square_vcf": {"keys": ['square_vcf']},
    "ploidy": {"keys": ['config', 'algorithm', 'ploidy'], "default": 2},
    "gender": {"keys": ["metadata", "sex"], "default": ""},
    "batch": {"keys": ["metadata", "batch"]},
    "bcbiornaseq": {"keys": ["config", "algorithm", "bcbiornaseq"], "default": {}},
    "mark_duplicates": {"keys": ["config", "algorithm", "mark_duplicates"], "default": True},
    "phenotype": {"keys": ["metadata", "phenotype"], "default": ""},
    "svclass": {"keys": ["metadata", "svclass"], "default": ""},
    "prep_method": {"keys": ["metadata", "prep_method"], "default": ""},
    "disease": {"keys": ["metadata", "disease"], "default": ""},
    "hetcaller": {"keys": ["config", "algorithm", "hetcaller"]},
    "variantcaller": {"keys": ['config', 'algorithm', 'variantcaller']},
    "variantcaller_order": {"keys": ['config', 'algorithm', 'variantcaller_order'], "default": 0},
    "svcaller": {"keys": ['config', 'algorithm', 'svcaller'], "default": [], "always_list": True},
    "jointcaller": {"keys": ['config', 'algorithm', 'jointcaller']},
    "hlacaller": {"keys": ['config', 'algorithm', 'hlacaller']},
    "recalibrate": {"keys": ['config', 'algorithm', 'recalibrate'], "default": False},
    "realign": {"keys": ['config', 'algorithm', 'realign'], "default": False},
    "ensemble": {"keys": ["config", "algorithm", "ensemble"], "default": {}},
    "background_variant": {"keys": ["config", "algorithm", "background", "variant"]},
    "peakcaller": {"keys": ['config', 'algorithm', 'peakcaller'], "default": []},
    "chip_method": {"keys": ['config', 'algorithm', 'chip_method'], "default": "chip"},
    "spikein_counts": {"keys": ["spikein_counts"]},
    "count_file": {"keys": ["count_file"]},
    "mirna_counts": {"keys": ["mirna_counts"]},
    "isomir_counts": {"keys": ["isomir_counts"]},
    "novel_mirna_counts": {"keys": ["novel_mirna_counts"]},
    "novel_isomir_counts": {"keys": ["novel_isomir_counts"]},
    "combined_counts": {"keys": ["combined_counts"]},
    "combined_histogram": {"keys": ["combined_histogram"]},
    "annotated_combined_counts": {"keys": ["annotated_combined_counts"]},
    "genome_context_files": {"keys": ["reference", "genome_context"], "default": [], "always_list": True},
    "viral_files": {"keys": ["reference", "viral"], "default": [], "always_list": True},
    "archive": {"keys": ["config", "algorithm", "archive"], "default": [], "always_list": True},
    "dexseq_gff": {"keys": ['genome_resources', 'rnaseq', 'dexseq']},
    "combined_fpkm": {"keys": ['combined_fpkm']},
    "combined_fpkm_isoform": {"keys": ['combined_fpkm_isoform']},
    "express_fpkm": {"keys": ['express_fpkm']},
    "express_tpm": {"keys": ['express_tpm']},
    "express_counts": {"keys": ['express_counts']},
    "histogram_counts": {"keys": ['histogram_counts']},
    "isoform_to_gene": {"keys": ['isoform_to_gene']},
    "fusion_mode": {"keys": ['config', 'algorithm', 'fusion_mode']},
    "fusion_caller": {"keys": ['config', 'algorithm', 'fusion_caller']},
    "dexseq_counts": {"keys": ['dexseq_counts']},
    "description": {"keys": ['description']},
    "aligner": {"keys": ['config', 'algorithm', 'aligner']},
    "align_split_size": {"keys": ['config', 'algorithm', 'align_split_size']},
    "bam_clean": {"keys": ['config', 'algorithm', 'bam_clean']},
    "platform": {"keys": ['config', 'algorithm', 'platform'],
                 "default": "illumina"},
    "quality_format": {"keys": ['config', 'algorithm', 'quality_format'],
                       "default": "standard"},
    "algorithm_qc": {"keys": ['config', 'algorithm', 'qc'], "default": [], "always_list": True},
    "summary_qc": {"keys": ['summary', 'qc'], "default": {}},
    "summary_metrics": {"keys": ['summary', 'metrics'], "default": {}},
    "adapters": {"keys": ['config', 'algorithm', 'adapters'], "default": [], "always_list": True},
    "custom_trim": {"keys": ['config', 'algorithm', 'custom_trim'],
                 "default": []},
    "species": {"keys": ['config', 'algorithm', 'species'],
                 "default": None},
    "trim_reads": {"keys": ['config', 'algorithm', 'trim_reads'],
                 "default": None},
    "trim_ends": {"keys": ['config', 'algorithm', 'trim_ends'],
                 "default": []},
    "min_read_length": {"keys": ['config', 'algorithm', 'min_read_length'],
                        "default": 25},
    "variation_resources": {"keys": ["genome_resources", "variation"], "default": {}},
    "qsig_file": {"keys": ['genome_resources', 'variation', 'qsignature'],
                  "checker": file_exists},
    "mixup_check": {"keys": ["config", "algorithm", "mixup_check"],
                    "default": False},
    "cufflinks_dir": {"keys": ['cufflinks_dir']},
    "stringtie_dir": {"keys": ['stringtie_dir']},
    "rsem": {"keys": ["config", "algorithm", "rsem"], "default": False},
    "transcriptome_align": {"keys": ["config", "algorithm", "transcriptome_align"],
                            "default": False},
    "quantify_genome_alignments": {"keys": ["config", "algorithm", "quantify_genome_alignments"],
                             "default": False},
    "expression_caller": {"keys": ["config", "algorithm", "expression_caller"],
                          "default": [], "always_list": True},
    "fusion_caller": {"keys": ["config", "algorithm", "fusion_caller"], "default": []},
    "spikein_fasta" : {"keys": ["config", "algorithm", "spikein_fasta"], "default": None},
    "transcriptome_bam": {"keys": ["transcriptome_bam"]},
    "junction_bed": {"keys": ["junction_bed"]},
    "starjunction": {"keys": ["starjunction"]},
    "chimericjunction": {"keys": ["chimericjunction"]},
    "fpkm_isoform": {"keys": ["fpkm_isoform"]},
    "fpkm": {"keys": ["fpkm"]},
    "galaxy_dir": {"keys": ["dirs", "galaxy"]},
    "assembled_gtf": {"keys": ["assembled_gtf"], "default": []},
    "merged_gtf": {"keys": ["merged_gtf"], "default": None},
    "transcript_assembler": {"keys": ["config", "algorithm", "transcript_assembler"],
                             "default": []},
    "oncofuse_file": {"keys": ["oncofuse_file"]},
    "pizzly_dir": {"keys": ["pizzly_dir"]},
    "split_bam": {"keys": ["split_bam"]},
    "vrn_file": {"keys": ["vrn_file"]},
    "exclude_regions": {"keys": ["config", "algorithm", "exclude_regions"], "default": [],
                        "always_list": True},
    "variant_regions": {"keys": ["config", "algorithm", "variant_regions"]},
    "variant_regions_merged": {"keys": ["config", "algorithm", "variant_regions_merged"]},
    "variant_regions_orig": {"keys": ["config", "algorithm", "variant_regions_orig"]},
    "sv_regions": {"keys": ["config", "algorithm", "sv_regions"]},
    "coverage": {"keys": ["config", "algorithm", "coverage"]},
    "coverage_merged": {"keys": ["config", "algorithm", "coverage_merged"]},
    "coverage_orig": {"keys": ["config", "algorithm", "coverage_orig"]},
    "callable_regions": {"keys": ["regions", "callable"]},
    "avg_coverage": {"keys": ["regions", "avg_coverage"]},
    "callable_min_size": {"keys": ["config", "algorithm", "callable_min_size"],
                          "default": 1000000},
    "min_allele_fraction": {"keys": ["config", "algorithm", "min_allele_fraction"], "default": 10},
    "normalized_depth": {"keys": ["depth", "bins", "normalized"]},

    "save_diskspace": {"keys": ["config", "algorithm", "save_diskspace"]},
    "salmon": {"keys": ["salmon"]},
    "umi_type": {"keys": ["config", "algorithm", "umi_type"]},
    "correct_umis": {"keys": ["config", "algorithm", "correct_umis"]},
    "sample_barcodes": {"keys": ["config", "algorithm", "sample_barcodes"]},
    "cellular_barcodes": {"keys": ["config", "algorithm", "cellular_barcodes"],
                          "default": []},
    "minimum_barcode_depth": {"keys": ["config", "algorithm", "minimum_barcode_depth"],
                              "default": 100000},
    "cellular_barcode_correction": {"keys": ["config", "algorithm",
                                             "cellular_barcode_correction"],
                                    "default": 1},
    "demultiplexed": {"keys": ["config", "algorithm", "demultiplexed"]},
    "kallisto_quant": {"keys": ["kallisto_quant"]},
    "salmon_dir": {"keys": ["salmon_dir"]},
    "salmon_fraglen_file": {"keys": ["salmon_fraglen_file"]},
    "sailfish": {"keys": ["sailfish"]},
    "sailfish_dir": {"keys": ["sailfish_dir"]},
    "sailfish_transcript_tpm": {"keys": ["sailfish_transcript_tpm"]},
    "sailfish_gene_tpm": {"keys": ["sailfish_gene_tpm"]},
    "sample_callable": {"keys": ["regions", "sample_callable"]},
    "coverage_interval": {"keys": ["config", "algorithm", "coverage_interval"]},
    "coverage_depth_min": {"keys": ["config", "algorithm", "coverage_depth_min"],
                           "default": 4},
    "maxcov_downsample": {"keys": ["config", "algorithm", "maxcov_downsample"],
                          "default": False},
    "joint_group_size": {"keys": ["config", "algorithm", "joint_group_size"],
                         "default": 200},
    "arriba": {"keys": ["arriba"], "default": {}},
    "report": {"keys": ["config", "algorithm", "report"]},
    "work_bam": {"keys": ["work_bam"]},
    "deduped_bam": {"keys": ["deduped_bam"]},
    "align_bam": {"keys": ["align_bam"]},
    "disc_bam": {"keys": ["work_bam_plus", "disc"]},
    "sr_bam": {"keys": ["work_bam_plus", "sr"]},
    "peddy_report": {"keys": ["peddy_report"]},
    "tools_off": {"keys": ["config", "algorithm", "tools_off"], "default": [], "always_list": True},
    "tools_on": {"keys": ["config", "algorithm", "tools_on"], "default": [], "always_list": True},
    "cwl_reporting": {"keys": ["config", "algorithm", "cwl_reporting"]},
}

def get_background_cnv_reference(data, caller):
    out = tz.get_in(["config", "algorithm", "background", "cnv_reference"], data)
    if out:
        return out.get(caller) if isinstance(out, dict) else out

def get_batches(data):
    batches = get_batch(data)
    if batches:
        if not isinstance(batches, (list, tuple)):
            batches = [batches]
        return batches
    return []

def get_input_sequence_files(data, default=None):
    """
    returns the input sequencing files, these can be single or paired FASTQ
    files or BAM files
    """
    if "files" not in data or data.get("files") is None:
        file1, file2 = None, None
    elif len(data["files"]) == 2:
        file1, file2 = data["files"]
    else:
        assert len(data["files"]) == 1, data["files"]
        file1, file2 = data["files"][0], None
    return file1, file2

def get_umi_consensus(data):
    """Retrieve UMI for consensus based preparation.

    We specify this either as a separate fastq file or embedded
    in the read name as `fastq_name`.`
    """
    consensus_choices = (["fastq_name"])
    umi = tz.get_in(["config", "algorithm", "umi_type"], data)
    # don't run consensus UMI calling for scrna-seq
    if tz.get_in(["analysis"], data, "").lower() == "scrna-seq":
        return False
    if umi and (umi in consensus_choices or os.path.exists(umi)):
        assert tz.get_in(["config", "algorithm", "mark_duplicates"], data, True), \
            "Using consensus UMI inputs requires marking duplicates"
        return umi

def get_correct_umis(data):
    """
    Do we need to correct UMIs with a whitelist?
    """
    umi_whitelist = tz.get_in(["config", "algorithm", "correct_umis"], data)
    if umi_whitelist:
        return True
    else:
        return False

def get_dexseq_gff(config, default=None):
    """
    some older versions of the genomes have the DEXseq gff file as
    gff instead of gff3, so this handles that by looking for either one
    """
    dexseq_gff = tz.get_in(tz.get_in(['dexseq_gff', 'keys'], LOOKUPS, {}),
                           config, None)
    if not dexseq_gff:
        return None
    gtf_file = get_gtf_file(config)
    if gtf_file:
        base_dir = os.path.dirname(gtf_file)
    else:
        base_dir = os.path.dirname(dexseq_gff)
    base, _ = os.path.splitext(dexseq_gff)
    gff_file = os.path.join(base_dir, base + ".gff")
    if file_exists(gff_file):
        return gff_file
    gtf_file = os.path.join(base_dir, base + ".gff3")
    if file_exists(gtf_file):
        return gtf_file
    else:
        return None

def getter(keys, global_default=None, always_list=False):
    def lookup(config, default=None):
        default = global_default if not default else default
        val = tz.get_in(keys, config, default)
        if always_list:
            if not val:
                val = []
            elif not isinstance(val, (list, tuple)): val = [val]
        return val
    return lookup

def setter(keys, checker):
    def update(config, value):
        if checker and not checker(value):
            logger.error("%s fails check %s." % (value, checker))
            sys.exit(1)
        return tz.update_in(config, keys, lambda x: value, default=value)
    return update

def is_setter(keys):
    def present(config):
        try:
            value = tz.get_in(keys, config, no_default=True)
        except:
            value = False
        return True if value else False
    return present

"""
generate the getter and setter functions but don't override any explicitly
defined
"""
_g = globals()
for k, v in LOOKUPS.items():
    keys = v['keys']
    getter_fn = 'get_' + k
    if getter_fn not in _g:
        _g["get_" + k] = getter(keys, v.get('default', None), v.get("always_list", False))
    setter_fn = 'set_' + k
    if setter_fn not in _g:
        _g["set_" + k] = setter(keys, v.get('checker', None))
    is_setter_fn = "is_set" + k
    if is_setter_fn not in _g:
        _g["is_set_" + k] = is_setter(keys)

def sample_data_iterator(samples):
    """
    for a list of samples, return the data dictionary of each sample
    """
    for sample in samples:
        yield to_single_data(sample)

def get_in_samples(samples, fn):
    """
    for a list of samples, return the value of a global option
    """
    for sample in samples:
        sample = to_single_data(sample)
        if fn(sample, None):
            return fn(sample)
    return None

def set_in_samples(samples, fn, value):
    """
    update a list of samples with a given value
    """
    return [fn(x, value) for x in samples]

def get_keys(lookup):
    """
    return the keys used to look up a function in the datadict
    """
    return tz.get_in((lookup, "keys"), LOOKUPS, None)

def update_summary_qc(data, key, base=None, secondary=None):
    """
    updates summary_qc with a new section, keyed by key.
    stick files into summary_qc if you want them propagated forward
    and available for multiqc
    """
    summary = get_summary_qc(data, {})
    if base and secondary:
        summary[key] = {"base": base, "secondary": secondary}
    elif base:
        summary[key] = {"base": base}
    elif secondary:
        summary[key] = {"secondary": secondary}
    data = set_summary_qc(data, summary)
    return data

def has_variantcalls(data):
    """
    returns True if the data dictionary is configured for variant calling
    """
    analysis = get_analysis(data).lower()
    variant_pipeline = analysis.startswith(("standard", "variant", "variant2"))
    variantcaller = get_variantcaller(data)
    return variant_pipeline or variantcaller

def get_algorithm_keys():
    """
    returns a list of all defined keys under the algorithm section
    """
    keys = []
    for k, v in LOOKUPS.items():
        if not "algorithm" in v["keys"]:
            continue
        if k == v["keys"][2]:
            keys.append(k)
    return keys
