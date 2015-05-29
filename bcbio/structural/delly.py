"""Structural variant calling with Delly

https://github.com/tobiasrausch/delly
"""
import copy
import itertools
import os
import re
import subprocess

import toolz as tz
import vcf

from bcbio import bam, utils
from bcbio.distributed.multi import run_multicore, zeromq_aware_logging
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.structural import shared as sshared
from bcbio.variation import vcfutils

def _get_full_exclude_file(items, work_dir):
    base_file = os.path.join(work_dir, "%s-svs" % (os.path.splitext(os.path.basename(items[0]["work_bam"]))[0]))
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
def _run_delly(bam_files, chrom, sv_type, ref_file, work_dir, items):
    """Run delly, calling structural variations for the specified type.
    """
    batch = sshared.get_cur_batch(items)
    ext = "-%s-svs" % batch if batch else "-svs"
    out_file = os.path.join(work_dir, "%s%s-%s-%s.vcf"
                            % (os.path.splitext(os.path.basename(bam_files[0]))[0], ext, chrom, sv_type))
    cores = min(utils.get_in(items[0], ("config", "algorithm", "num_cores"), 1),
                len(bam_files))
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            names = [tz.get_in(["rgnames", "sample"], x) for x in items]
            if not sshared.has_variant_regions(items, out_file, chrom):
                vcfutils.write_empty_vcf(tx_out_file, samples=names)
            else:
                exclude = ["-x", _delly_exclude_file(items, out_file, chrom)]
                map_qual = ["--map-qual", "1"]
                cmd = ["delly", "-t", sv_type, "-g", ref_file, "-o", tx_out_file] + exclude + bam_files
                multi_cmd = "export OMP_NUM_THREADS=%s && export LC_ALL=C && " % cores
                try:
                    do.run(multi_cmd + " ".join(cmd), "delly structural variant")
                    # Delly will write nothing if no variants found
                    if not utils.file_exists(tx_out_file):
                        vcfutils.write_empty_vcf(tx_out_file, samples=names)
                except subprocess.CalledProcessError, msg:
                    # delly returns an error exit code if there are no variants
                    if "No structural variants found" in str(msg):
                        vcfutils.write_empty_vcf(tx_out_file, samples=names)
                    else:
                        raise
    return [vcfutils.bgzip_and_index(_clean_delly_output(out_file, items), items[0]["config"])]

def _clean_delly_output(in_file, items):
    """Clean delly output, fixing sample names and removing problem GL specifications from output.

    GATK does not like missing GLs like '.,.,.'. This converts them to the recognized '.'
    """
    pat = re.compile(r"\.,\.,\.")
    out_file = "%s-clean.vcf" % utils.splitext_plus(in_file)[0]
    if not utils.file_exists(out_file) and not utils.file_exists(out_file + ".gz"):
        with file_transaction(items[0], out_file) as tx_out_file:
            with open(in_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    for line in in_handle:
                        if line.startswith("#"):
                            if line.startswith("#CHROM"):
                                line = _fix_sample_names(line, items)
                        else:
                            line = pat.sub(".", line)
                        out_handle.write(line)
    return out_file

def _fix_sample_names(line, items):
    """Substitute Delly output sample names (filenames) with actual sample names.
    """
    names = [tz.get_in(["rgnames", "sample"], x) for x in items]
    parts = line.split("\t")
    # If we're not empty and actually have genotype information
    if "FORMAT" in parts:
        format_i = parts.index("FORMAT") + 1
        assert len(parts[format_i:]) == len(names), (parts[format_i:], names)
        return "\t".join(parts[:format_i] + names) + "\n"
    else:
        return line

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

    XXX Currently does not downsample as new versions of delly do not get good sensitivity
    with downsampled BAMs.
    """
    return [data["align_bam"]]

    # full_bam, sr_bam, disc_bam = sshared.get_split_discordants(data, work_dir)
    # ds_bam = bam.downsample(full_bam, data, 1e8, read_filter="-F 'not secondary_alignment and proper_pair'",
    #                         always_run=True, work_dir=work_dir)
    # out_bam = "%s-final%s" % utils.splitext_plus(ds_bam)
    # if not utils.file_exists(out_bam):
    #     bam.merge([ds_bam, sr_bam, disc_bam], out_bam, data["config"])
    # bam.index(out_bam, data["config"])
    # return [out_bam]

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
                                               items[0]["name"][-1], "delly"))
    # Add core request for delly
    config = copy.deepcopy(items[0]["config"])
    delly_config = utils.get_in(config, ("resources", "delly"), {})
    delly_config["cores"] = 1
    config["resources"]["delly"] = delly_config
    parallel = {"type": "local", "cores": config["algorithm"].get("num_cores", 1),
                "progs": ["delly"]}
    work_bams = run_multicore(_prep_subsampled_bams,
                              [(data, work_dir) for data in items],
                              config, parallel)
    ref_file = utils.get_in(items[0], ("reference", "fasta", "base"))
    sv_types = ["DEL", "DUP"]  # "TRA" has invalid VCF END specifications that GATK doesn't like, "INV" very slow
    exclude_file = _get_full_exclude_file(items, work_dir)
    bytype_vcfs = run_multicore(_run_delly,
                                [(work_bams, chrom, sv_type, ref_file, work_dir, items)
                                 for (chrom, sv_type)
                                 in itertools.product(sshared.get_sv_chroms(items, exclude_file), sv_types)],
                                config, parallel)
    out_file = "%s.vcf.gz" % sshared.outname_from_inputs(bytype_vcfs)
    combo_vcf = vcfutils.combine_variant_files(bytype_vcfs, out_file, ref_file, config)
    out = []
    for data in items:
        if "sv" not in data:
            data["sv"] = []
        base, ext = utils.splitext_plus(combo_vcf)
        sample = tz.get_in(["rgnames", "sample"], data)
        delly_sample_vcf = vcfutils.select_sample(combo_vcf, sample,
                                                  "%s-%s%s" % (base, sample, ext), data["config"])
        delly_vcf = _delly_count_evidence_filter(delly_sample_vcf, data)
        data["sv"].append({"variantcaller": "delly", "vrn_file": delly_vcf,
                           "exclude": exclude_file})
        out.append(data)
    return out
