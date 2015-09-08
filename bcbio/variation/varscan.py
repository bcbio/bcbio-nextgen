"""Provide variant calling with VarScan from TGI at Wash U.

http://varscan.sourceforge.net/
"""

from collections import namedtuple
import contextlib
from distutils.version import LooseVersion
import os
import shutil

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.pipeline import config_utils
from bcbio.provenance import do, programs
from bcbio.utils import file_exists, append_stem
from bcbio.variation import freebayes, samtools, vcfutils
from bcbio.variation.vcfutils import (combine_variant_files, write_empty_vcf,
                                      get_paired_bams, is_paired_analysis,
                                      bgzip_and_index, move_vcf)

import pysam
import vcf


def run_varscan(align_bams, items, ref_file, assoc_files,
                region=None, out_file=None):
    paired = get_paired_bams(align_bams, items)
    if paired and paired.normal_bam and paired.tumor_bam:
        call_file = samtools.shared_variantcall(_varscan_paired, "varscan",
                                                align_bams, ref_file, items,
                                                assoc_files, region, out_file)
    else:
        vcfutils.check_paired_problems(items)
        call_file = samtools.shared_variantcall(_varscan_work, "varscan",
                                                align_bams, ref_file,
                                                items, assoc_files,
                                                region, out_file)
    return call_file


def _get_varscan_opts(config, tmp_dir):
    """Retrieve common options for running VarScan.
    Handles jvm_opts, setting user and country to English to avoid issues
    with different locales producing non-compliant VCF.
    """
    resources = config_utils.get_resources("varscan", config)
    jvm_opts = resources.get("jvm_opts", ["-Xmx750m", "-Xmx2g"])
    jvm_opts = config_utils.adjust_opts(jvm_opts,
                                        {"algorithm": {"memory_adjust":
                                                       {"magnitude": 1.1, "direction": "decrease"}}})
    jvm_opts += ["-Duser.language=en", "-Duser.country=US"]
    jvm_opts += broad.get_default_jvm_opts(tmp_dir)
    return " ".join(jvm_opts)


def _varscan_paired(align_bams, ref_file, items, target_regions, out_file):

    """Run a paired VarScan analysis, also known as "somatic". """

    max_read_depth = "1000"
    config = items[0]["config"]

    version = programs.jar_versioner("varscan", "VarScan")(config)
    if LooseVersion(version) < LooseVersion("v2.3.6"):
        raise IOError(
            "Please install version 2.3.6 or better of VarScan with support "
            "for multisample calling and indels in VCF format.")
    varscan_jar = config_utils.get_jar(
        "VarScan",
        config_utils.get_program("varscan", config, "dir"))

    remove_zerocoverage = "grep -v -P '\t0\t\t$'"

    # No need for names in VarScan, hence the "_"

    paired = get_paired_bams(align_bams, items)
    if not paired.normal_bam:
        affected_batch = items[0]["metadata"]["batch"]
        message = ("Batch {} requires both tumor and normal BAM files for"
            " VarScan cancer calling").format(affected_batch)
        raise ValueError(message)

    if not file_exists(out_file):
        orig_out_file = out_file
        out_file = orig_out_file.replace(".vcf.gz", ".vcf")
        base, ext = utils.splitext_plus(out_file)
        cleanup_files = []
        for fname, mpext in [(paired.normal_bam, "normal"), (paired.tumor_bam, "tumor")]:
            mpfile = "%s-%s.mpileup" % (base, mpext)
            cleanup_files.append(mpfile)
            with file_transaction(config, mpfile) as mpfile_tx:
                mpileup = samtools.prep_mpileup([fname], ref_file,
                                                config, max_read_depth,
                                                target_regions=target_regions,
                                                want_bcf=False)
                cmd = "{mpileup} > {mpfile_tx}"
                cmd = cmd.format(**locals())
                do.run(cmd, "samtools mpileup".format(**locals()), None,
                       [do.file_exists(mpfile_tx)])

        # Sometimes mpileup writes an empty file: in this case we
        # just skip the rest of the analysis (VarScan will hang otherwise)

        if any(os.stat(filename).st_size == 0 for filename in cleanup_files):
            write_empty_vcf(orig_out_file, config)
            return

        # First index is normal, second is tumor
        normal_tmp_mpileup = cleanup_files[0]
        tumor_tmp_mpileup = cleanup_files[1]

        indel_file = base + ".indel.vcf"
        snp_file = base + ".snp.vcf"
        cleanup_files.append(indel_file)
        cleanup_files.append(snp_file)
        with file_transaction(config, indel_file, snp_file) as (tx_indel, tx_snp):
            with tx_tmpdir(items[0]) as tmp_dir:
                jvm_opts = _get_varscan_opts(config, tmp_dir)
                fix_ambig = vcfutils.fix_ambiguous_cl()
                tx_snp_in = "%s-orig" % os.path.splitext(tx_snp)[0]
                tx_indel_in = "%s-orig" % os.path.splitext(tx_indel)[0]
                varscan_cmd = ("java {jvm_opts} -jar {varscan_jar} somatic"
                               " {normal_tmp_mpileup} {tumor_tmp_mpileup} "
                               "--output-snp {tx_snp_in} --output-indel {tx_indel_in} "
                               " --output-vcf --min-coverage 5 --p-value 0.98 "
                               "--strand-filter 1 ")
                # add minimum AF
                if "--min-var-freq" not in varscan_cmd:
                    min_af = float(utils.get_in(paired.tumor_config, ("algorithm",
                                                                      "min_allele_fraction"),10)) / 100.0
                    varscan_cmd += "--min-var-freq {min_af} "
                do.run(varscan_cmd.format(**locals()), "Varscan", None, None)
                for orig_fname, fname in [(tx_snp_in, tx_snp), (tx_indel_in, tx_indel)]:
                    cmd = "vcfuniqalleles {orig_fname}.vcf | {fix_ambig} > {fname}"
                    do.run(cmd.format(**locals()), "Varscan paired fix")

        # VarScan files need to be corrected to match the VCF specification
        # We do this before combining them otherwise merging may fail
        # if there are invalid records
        to_combine = []
        if do.file_exists(snp_file):
            to_combine.append(snp_file)
            _fix_varscan_vcf(snp_file, paired.normal_name, paired.tumor_name, config)

        if do.file_exists(indel_file):
            to_combine.append(indel_file)
            _fix_varscan_vcf(indel_file, paired.normal_name, paired.tumor_name, config)

        if not to_combine:
            write_empty_vcf(orig_out_file, config)
            return

        out_file = combine_variant_files([snp_file, indel_file],
                                         out_file, ref_file, config,
                                         region=target_regions)

        # Remove cleanup files

        for extra_file in cleanup_files:
            for ext in ["", ".gz", ".gz.tbi"]:
                if os.path.exists(extra_file + ext):
                    os.remove(extra_file + ext)

        if os.path.getsize(out_file) == 0:
            write_empty_vcf(out_file)

        if orig_out_file.endswith(".gz"):
            out_file = bgzip_and_index(out_file, config)

        _add_reject_flag(out_file, config)

def _fix_varscan_vcf(orig_file, normal_name, tumor_name, config):
    """Fixes issues with the standard VarScan VCF output.

    - Remap sample names back to those defined in the input BAM file.
    - Convert indels into correct VCF representation.
    """
    tmp_file = append_stem(orig_file, "-origsample")

    if not file_exists(tmp_file):
        shutil.move(orig_file, tmp_file)

        with file_transaction(config, orig_file) as tx_out_file:
            with open(tmp_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:

                    for line in in_handle:
                        line = _clean_varscan_line(_fix_varscan_output(line, normal_name, tumor_name))
                        if not line:
                            continue
                        out_handle.write(line)

def _add_reject_flag(in_file, config):

    """Add REJECT flag to all records that aren't flagged somatic
    (SS=2)"""

    Filter = namedtuple('Filter', ['id', 'desc'])
    reject_filter = Filter(id='REJECT',
                           desc='Rejected as non-SOMATIC or by quality')
    # NOTE: PyVCF will write an uncompressed VCF
    base, ext = utils.splitext_plus(in_file)
    name = "rejectfix"
    out_file = "{0}-{1}{2}".format(base, name, ".vcf")

    if utils.file_exists(in_file):
        reader = vcf.VCFReader(filename=in_file)
        # Add info to the header of the reader
        reader.filters["REJECT"] = reject_filter
        with file_transaction(config, out_file) as tx_out_file:
            with open(tx_out_file, "wb") as handle:
                writer = vcf.VCFWriter(handle, template=reader)
                for record in reader:
                    if "SS" in record.INFO:
                        # VarScan encodes it as a string
                        # TODO: Set it as integer when cleaning

                        if record.INFO["SS"] != "2":
                            record.add_filter("REJECT")
                    writer.write_record(record)

        # Re-compress the file
        out_file = bgzip_and_index(out_file, config)
        move_vcf(in_file, "{0}.orig".format(in_file))
        move_vcf(out_file, in_file)
        with open(out_file, "w") as out_handle:
            out_handle.write("Moved to {0}".format(in_file))


def _fix_varscan_output(line, normal_name, tumor_name):
    """Fix a varscan VCF line

    Fixes the ALT column and also fixes floating point values
    output as strings to by Floats: FREQ, SSC.

    :param line: a pre-split and stripped varscan line

    This function was contributed by Sean Davis <sdavis2@mail.nih.gov>,
    with minor modifications by Luca Beltrame <luca.beltrame@marionegri.it>.
    """
    line = line.strip()

    tofix = ("##INFO=<ID=SSC", "##FORMAT=<ID=FREQ")
    if(line.startswith("##")):
        if line.startswith(tofix):
            line = line.replace('Number=1,Type=String',
                                'Number=1,Type=Float')
        return line + "\n"

    line = line.split("\t")

    mapping = {"NORMAL": normal_name, "TUMOR": tumor_name}

    if(line[0].startswith("#CHROM")):

        base_header = line[:9]
        old_samples = line[9:]

        if len(old_samples) == 0:
            return "\t".join(line) + "\n"

        samples = [mapping[sample_name] for sample_name in old_samples]

        assert len(old_samples) == len(samples)
        return "\t".join(base_header + samples) + "\n"

    try:
        REF, ALT = line[3:5]
    except ValueError:
        return "\t".join(line) + "\n"

    Ifreq = line[8].split(":").index("FREQ")
    ndat = line[9].split(":")
    tdat = line[10].split(":")
    somatic_status = line[7].split(";")  # SS=<number>
    # HACK: The position of the SS= changes, so we just search for it
    somatic_status = [item for item in somatic_status
                      if item.startswith("SS=")][0]
    somatic_status = int(somatic_status.split("=")[1])  # Get the number

    ndat[Ifreq] = str(float(ndat[Ifreq].rstrip("%")) / 100)
    tdat[Ifreq] = str(float(tdat[Ifreq].rstrip("%")) / 100)
    line[9] = ":".join(ndat)
    line[10] = ":".join(tdat)

    #FIXME: VarScan also produces invalid REF records (e.g. CAA/A)
    # This is not handled yet.

    if somatic_status == 5:

        # "Unknown" states are broken in current versions of VarScan
        # so we just bail out here for now

        return

    if "+" in ALT or "-" in ALT:
        if "/" not in ALT:
            if ALT[0] == "+":
                R = REF
                A = REF + ALT[1:]
            elif ALT[0] == "-":
                R = REF + ALT[1:]
                A = REF
        else:
            Ins = [p[1:] for p in ALT.split("/") if p[0] == "+"]
            Del = [p[1:] for p in ALT.split("/") if p[0] == "-"]

            if len(Del):
                REF += sorted(Del, key=lambda x: len(x))[-1]

            A = ",".join([REF[::-1].replace(p[::-1], "", 1)[::-1]
                          for p in Del] + [REF + p for p in Ins])
            R = REF

        REF = R
        ALT = A
    else:
        ALT = ALT.replace('/', ',')

    line[3] = REF
    line[4] = ALT
    return "\t".join(line) + "\n"


def _create_sample_list(in_bams, vcf_file):
    """Pull sample names from input BAMs and create input sample list.
    """
    out_file = "%s-sample_list.txt" % os.path.splitext(vcf_file)[0]
    with open(out_file, "w") as out_handle:
        for in_bam in in_bams:
            with contextlib.closing(pysam.Samfile(in_bam, "rb")) as work_bam:
                for rg in work_bam.header.get("RG", []):
                    out_handle.write("%s\n" % rg["SM"])
    return out_file


def _varscan_work(align_bams, ref_file, items, target_regions, out_file):
    """Perform SNP and indel genotyping with VarScan.
    """
    config = items[0]["config"]

    orig_out_file = out_file
    out_file = orig_out_file.replace(".vcf.gz", ".vcf")

    max_read_depth = "1000"
    version = programs.jar_versioner("varscan", "VarScan")(config)
    if version < "v2.3.6":
        raise IOError("Please install version 2.3.6 or better of VarScan"
                      " with support for multisample calling and indels"
                      " in VCF format.")
    varscan_jar = config_utils.get_jar("VarScan",
                                       config_utils.get_program("varscan", config, "dir"))
    sample_list = _create_sample_list(align_bams, out_file)
    mpileup = samtools.prep_mpileup(align_bams, ref_file, config, max_read_depth,
                                    target_regions=target_regions, want_bcf=False)
    # VarScan fails to generate a header on files that start with
    # zerocoverage calls; strip these with grep, we're not going to
    # call on them
    remove_zerocoverage = "grep -v -P '\t0\t\t$'"
    # write a temporary mpileup file so we can check if empty
    mpfile = "%s.mpileup" % os.path.splitext(out_file)[0]
    with file_transaction(config, mpfile) as mpfile_tx:
        cmd = ("{mpileup} | {remove_zerocoverage} > {mpfile_tx}")
        do.run(cmd.format(**locals()), "mpileup for Varscan")
    if os.path.getsize(mpfile) == 0:
        write_empty_vcf(out_file)
    else:
        with tx_tmpdir(items[0]) as tmp_dir:
            jvm_opts = _get_varscan_opts(config, tmp_dir)
            fix_ambig = vcfutils.fix_ambiguous_cl()
            cmd = ("cat {mpfile} "
                   "| java {jvm_opts} -jar {varscan_jar} mpileup2cns --min-coverage 5 --p-value 0.98 "
                   "  --vcf-sample-list {sample_list} --output-vcf --variants "
                   "| {fix_ambig} | vcfuniqalleles > {out_file}")
            do.run(cmd.format(**locals()), "Varscan", None,
                   [do.file_exists(out_file)])
    os.remove(sample_list)
    os.remove(mpfile)
    # VarScan can create completely empty files in regions without
    # variants, so we create a correctly formatted empty file
    if os.path.getsize(out_file) == 0:
        write_empty_vcf(out_file)
    else:
        freebayes.clean_vcf_output(out_file, _clean_varscan_line, config)

    if orig_out_file.endswith(".gz"):
        vcfutils.bgzip_and_index(out_file, config)


def _clean_varscan_line(line):
    """Avoid lines with non-GATC bases, ambiguous output bases make GATK unhappy.
    """
    if line and not line.startswith("#"):
        parts = line.split("\t")
        alleles = [x.strip() for x in parts[4].split(",")] + [parts[3].strip()]
        for a in alleles:
            if len(set(a) - set("GATCgatc")) > 0:
                return None
    return line
