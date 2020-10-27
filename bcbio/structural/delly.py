"""Structural variant calling with Delly

https://github.com/tobiasrausch/delly
"""
import collections
import copy
import os
import subprocess

import vcf

from bcbio import bam, utils
from bcbio.distributed.multi import run_multicore, zeromq_aware_logging
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.structural import shared as sshared
from bcbio.variation import vcfutils

def _get_full_exclude_file(items, work_bams, work_dir):
    base_file = os.path.join(work_dir, "%s-svs" % (os.path.splitext(os.path.basename(work_bams[0]))[0]))
    return sshared.prepare_exclude_file(items, base_file)

def _delly_exclude_file(items, base_file, chrom):
    """Prepare a delly-specific exclude file eliminating chromosomes.
    Delly wants excluded chromosomes listed as just the chromosome, with no coordinates.
    """
    base_exclude = sshared.prepare_exclude_file(items, base_file, chrom)
    out_file = "%s-delly%s" % utils.splitext_plus(base_exclude)
    with file_transaction(items[0], out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            with open(base_exclude) as in_handle:
                for line in in_handle:
                    parts = line.split("\t")
                    if parts[0] == chrom:
                        out_handle.write(line)
                    else:
                        out_handle.write("%s\n" % parts[0])
    return out_file

@utils.map_wrap
@zeromq_aware_logging
def _run_delly(bam_files, chrom, ref_file, work_dir, items):
    """Run delly, calling structural variations for the specified type.
    """
    batch = sshared.get_cur_batch(items)
    ext = "-%s-svs" % batch if batch else "-svs"
    out_file = os.path.join(work_dir, "%s%s-%s.bcf"
                            % (os.path.splitext(os.path.basename(bam_files[0]))[0], ext, chrom))
    final_file = "%s.vcf.gz" % (utils.splitext_plus(out_file)[0])
    cores = min(utils.get_in(items[0], ("config", "algorithm", "num_cores"), 1),
                len(bam_files))
    if not utils.file_exists(out_file) and not utils.file_exists(final_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            if sshared.has_variant_regions(items, out_file, chrom):
                exclude = ["-x", _delly_exclude_file(items, out_file, chrom)]
                cmd = ["delly", "call", "-g", ref_file, "-o", tx_out_file] + exclude + bam_files
                locale_to_use = utils.get_locale()
                multi_cmd = "export OMP_NUM_THREADS=%s && export LC_ALL=%s && " % (cores, locale_to_use)
                try:
                    do.run(multi_cmd + " ".join(cmd), "delly structural variant")
                except subprocess.CalledProcessError as msg:
                    # Small input samples, write an empty vcf
                    if "Sample has not enough data to estimate library parameters" in str(msg):
                        pass
                    # delly returns an error exit code if there are no variants
                    elif "No structural variants found" not in str(msg):
                        raise
    return [_bgzip_and_clean(out_file, items)]

def _bgzip_and_clean(bcf_file, items):
    """Create a clean bgzipped VCF output file from bcf for downstream processing.

    Also corrects problems with missing likelihoods: https://github.com/dellytools/delly/issues/37
    GATK does not like missing GLs like '.,.,.'. This converts them to the recognized '.'
    """
    out_file = "%s.vcf.gz" % (utils.splitext_plus(bcf_file)[0])
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            if not utils.file_exists(bcf_file):
                vcfutils.write_empty_vcf(tx_out_file, samples=[dd.get_sample_name(d) for d in items])
            else:
                cmd = ("bcftools view {bcf_file} | sed 's/\.,\.,\././' | bgzip -c > {tx_out_file}")
                do.run(cmd.format(**locals()), "Convert and clean delly output")
    return vcfutils.bgzip_and_index(out_file, items[0]["config"])

@utils.map_wrap
@zeromq_aware_logging
def _prep_subsampled_bams(data, work_dir):
    """Prepare a subsampled BAM file with discordants from samblaster and minimal correct pairs.

    This attempts to minimize run times by pre-extracting useful reads mixed
    with subsampled normal pairs to estimate paired end distributions:

    https://groups.google.com/d/msg/delly-users/xmia4lwOd1Q/uaajoBkahAIJ

    Subsamples correctly aligned reads to 100 million based on speedseq defaults and
    evaluations on NA12878 whole genome data:

    https://github.com/cc2qe/speedseq/blob/ca624ba9affb0bd0fb88834ca896e9122639ec94/bin/speedseq#L1102

    XXX Currently not used as new versions of delly do not get good sensitivity
    with downsampled BAMs.
    """
    sr_bam, disc_bam = sshared.get_split_discordants(data, work_dir)
    ds_bam = bam.downsample(dd.get_align_bam(data), data, 1e8,
                            read_filter="-F 'not secondary_alignment and proper_pair'",
                            always_run=True, work_dir=work_dir)
    out_bam = "%s-final%s" % utils.splitext_plus(ds_bam)
    if not utils.file_exists(out_bam):
        bam.merge([ds_bam, sr_bam, disc_bam], out_bam, data["config"])
    bam.index(out_bam, data["config"])
    return [out_bam]

def _delly_count_evidence_filter(in_file, data):
    """Filter delly outputs based on read support (DV) and evidence (split and paired).

    We require DV > 4 and either both paired end and split read evidence or
    5 or more evidence for either individually.
    """
    filtname = "DVSupport"
    filtdoc = "FMT/DV < 4 || (SR < 1 && PE < 5) || (SR < 5 && PE < 1)"
    out_file = "%s-filter%s" % utils.splitext_plus(in_file)
    cur_out_file = out_file.replace(".vcf.gz", ".vcf")
    if not utils.file_exists(out_file):
        with file_transaction(data, cur_out_file) as tx_out_file:
            with utils.open_gzipsafe(in_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    inp = vcf.Reader(in_handle, in_file)
                    inp.filters["DVSupport"] = vcf.parser._Filter(filtname, filtdoc)
                    outp = vcf.Writer(out_handle, inp)
                    for rec in inp:
                        sr = rec.INFO.get("SR", 0)
                        pe = rec.INFO.get("PE", 0)
                        call = rec.samples[0].data
                        dv = call.DV if hasattr(call, "DV") else 0
                        if dv < 4 or (sr < 1 and pe < 5) or (sr < 5 and pe < 1):
                            rec.add_filter(filtname)
                        outp.write_record(rec)
    if out_file.endswith(".vcf.gz"):
        out_file = vcfutils.bgzip_and_index(cur_out_file, data["config"])
    return out_file

def run(items):
    """Perform detection of structural variations with delly.

    Performs post-call filtering with a custom filter tuned based
    on NA12878 Moleculo and PacBio data, using calls prepared by
    @ryanlayer and @cc2qe

    Filters using the high quality variant pairs (DV) compared with
    high quality reference pairs (DR).
    """
    work_dir = utils.safe_makedir(os.path.join(items[0]["dirs"]["work"], "structural",
                                               dd.get_sample_name(items[0]), "delly"))
    # Add core request for delly
    config = copy.deepcopy(items[0]["config"])
    delly_config = utils.get_in(config, ("resources", "delly"), {})
    delly_config["cores"] = 1
    config["resources"]["delly"] = delly_config
    parallel = {"type": "local", "cores": config["algorithm"].get("num_cores", 1),
                "progs": ["delly"]}
    work_bams = [dd.get_align_bam(d) for d in items]
    ref_file = dd.get_ref_file(items[0])
    exclude_file = _get_full_exclude_file(items, work_bams, work_dir)
    bytype_vcfs = run_multicore(_run_delly,
                                [(work_bams, chrom, ref_file, work_dir, items)
                                 for chrom in sshared.get_sv_chroms(items, exclude_file)],
                                config, parallel)
    out_file = "%s.vcf.gz" % sshared.outname_from_inputs(bytype_vcfs)
    combo_vcf = vcfutils.combine_variant_files(bytype_vcfs, out_file, ref_file, config)
    out = []
    upload_counts = collections.defaultdict(int)
    for data in items:
        if "sv" not in data:
            data["sv"] = []
        base, ext = utils.splitext_plus(combo_vcf)
        final_vcf = sshared.finalize_sv(combo_vcf, data, items)
        if final_vcf:
            delly_vcf = _delly_count_evidence_filter(final_vcf, data)
            data["sv"].append({"variantcaller": "delly", "vrn_file": delly_vcf,
                               "do_upload": upload_counts[final_vcf] == 0,  # only upload a single file per batch
                               "exclude": exclude_file})
            upload_counts[final_vcf] += 1
        out.append(data)
    return out
