"""Provide variant calling with VarScan from TGI at Wash U.

http://varscan.sourceforge.net/
"""

import contextlib
import os
import sys

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.variation import samtools, vcfutils
from bcbio.variation.vcfutils import (combine_variant_files, write_empty_vcf,
                                      get_paired_bams, bgzip_and_index)

import pysam


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

def _safe_to_float(x):
    if x is None:
        return None
    else:
        try:
            return float(x)
        except ValueError:
            return None

def spv_freq_filter(line, tumor_index):
    """Filter VarScan calls based on the SPV value and frequency.

    Removes calls with SPV < 0.05 and a tumor FREQ > 0.35.

    False positives dominate these higher frequency, low SPV calls. They appear
    to be primarily non-somatic/germline variants not removed by other filters.
    """
    if line.startswith("#CHROM"):
        headers = [('##FILTER=<ID=SpvFreq,Description="High frequency (tumor FREQ > 0.35) '
                    'and low p-value for somatic (SPV < 0.05)">')]
        return "\n".join(headers) + "\n" + line
    elif line.startswith("#"):
        return line
    else:
        parts = line.split("\t")
        sample_ft = {a: v for (a, v) in zip(parts[8].split(":"), parts[9 + tumor_index].split(":"))}
        freq = _safe_to_float(sample_ft.get("FREQ"))
        spvs = [x for x in parts[7].split(";") if x.startswith("SPV=")]
        spv = _safe_to_float(spvs[0].split("=")[-1] if spvs else None)
        fname = None
        if spv is not None and freq is not None:
            if spv < 0.05 and freq > 0.35:
                fname = "SpvFreq"
        if fname:
            if parts[6] in set([".", "PASS"]):
                parts[6] = fname
            else:
                parts[6] += ";%s" % fname
        line = "\t".join(parts)
        return line

def _varscan_paired(align_bams, ref_file, items, target_regions, out_file):

    """Run a paired VarScan analysis, also known as "somatic". """

    max_read_depth = "1000"
    config = items[0]["config"]
    paired = get_paired_bams(align_bams, items)
    if not paired.normal_bam:
        affected_batch = items[0]["metadata"]["batch"]
        message = ("Batch {} requires both tumor and normal BAM files for"
                   " VarScan cancer calling").format(affected_batch)
        raise ValueError(message)

    if not utils.file_exists(out_file):
        assert out_file.endswith(".vcf.gz"), "Expect bgzipped output to VarScan"
        normal_mpileup_cl = samtools.prep_mpileup([paired.normal_bam], ref_file,
                                                  config, max_read_depth,
                                                  target_regions=target_regions,
                                                  want_bcf=False)
        tumor_mpileup_cl = samtools.prep_mpileup([paired.tumor_bam], ref_file,
                                                 config, max_read_depth,
                                                 target_regions=target_regions,
                                                 want_bcf=False)
        base, ext = utils.splitext_plus(out_file)
        indel_file = base + "-indel.vcf"
        snp_file = base + "-snp.vcf"
        with file_transaction(config, indel_file, snp_file) as (tx_indel, tx_snp):
            with tx_tmpdir(items[0]) as tmp_dir:
                jvm_opts = _get_varscan_opts(config, tmp_dir)
                remove_zerocoverage = r"{ ifne grep -v -P '\t0\t\t$' || true; }"
                varscan_cmd = ("varscan {jvm_opts} somatic "
                               " <({normal_mpileup_cl} | {remove_zerocoverage}) "
                               "<({tumor_mpileup_cl} | {remove_zerocoverage}) "
                               "--output-snp {tx_snp} --output-indel {tx_indel} "
                               " --output-vcf --min-coverage 5 --p-value 0.98 "
                               "--strand-filter 1 ")
                # add minimum AF
                if "--min-var-freq" not in varscan_cmd:
                    min_af = float(utils.get_in(paired.tumor_config, ("algorithm",
                                                                      "min_allele_fraction"), 10)) / 100.0
                    varscan_cmd += "--min-var-freq {min_af} "
                do.run(varscan_cmd.format(**locals()), "Varscan", None, None)

        to_combine = []
        for fname in [snp_file, indel_file]:
            if utils.file_exists(fname):
                fix_file = "%s-fix.vcf.gz" % (utils.splitext_plus(fname)[0])
                with file_transaction(config, fix_file) as tx_fix_file:
                    fix_ambig_ref = vcfutils.fix_ambiguous_cl()
                    fix_ambig_alt = vcfutils.fix_ambiguous_cl(5)
                    py_cl = os.path.join(os.path.dirname(sys.executable), "py")
                    normal_name = paired.normal_name
                    tumor_name = paired.tumor_name
                    cmd = ("cat {fname} | "
                           "{py_cl} -x 'bcbio.variation.varscan.fix_varscan_output(x,"
                            """ "{normal_name}", "{tumor_name}")' | """
                           "{fix_ambig_ref} | {fix_ambig_alt} | ifne vcfuniqalleles | "
                           """bcftools filter -m + -s REJECT -e "SS != '.' && SS != '2'" 2> /dev/null | """
                           "{py_cl} -x 'bcbio.variation.varscan.spv_freq_filter(x, 1)' | "
                           "bgzip -c > {tx_fix_file}")
                    do.run(cmd.format(**locals()), "Varscan paired fix")
                to_combine.append(fix_file)

        if not to_combine:
            out_file = write_empty_vcf(out_file, config)
        else:
            out_file = combine_variant_files(to_combine,
                                             out_file, ref_file, config,
                                             region=target_regions)
        if os.path.getsize(out_file) == 0:
            write_empty_vcf(out_file)
        if out_file.endswith(".gz"):
            out_file = bgzip_and_index(out_file, config)

def fix_varscan_output(line, normal_name="", tumor_name=""):
    """Fix a varscan VCF line.

    Fixes the ALT column and also fixes floating point values
    output as strings to by Floats: FREQ, SSC.

    This function was contributed by Sean Davis <sdavis2@mail.nih.gov>,
    with minor modifications by Luca Beltrame <luca.beltrame@marionegri.it>.
    """
    line = line.strip()

    tofix = ("##INFO=<ID=SSC", "##FORMAT=<ID=FREQ")
    if(line.startswith("##")):
        if line.startswith(tofix):
            line = line.replace('Number=1,Type=String',
                                'Number=1,Type=Float')
        return line
    line = line.split("\t")

    if line[0].startswith("#CHROM"):
        if tumor_name and normal_name:
            mapping = {"NORMAL": normal_name, "TUMOR": tumor_name}
            base_header = line[:9]
            old_samples = line[9:]

            if len(old_samples) == 0:
                return "\t".join(line)

            samples = [mapping[sample_name] for sample_name in old_samples]

            assert len(old_samples) == len(samples)
            return "\t".join(base_header + samples)
        else:
            return "\t".join(line)

    try:
        REF, ALT = line[3:5]
    except ValueError:
        return "\t".join(line)

    if len(line) > 10:
        #print line
        Ifreq = line[8].split(":").index("FREQ")
        #print repr(Ifreq)
        ndat = line[9].split(":")
        tdat = line[10].split(":")
        #print ndat
        #print tdat
        somatic_status = line[7].split(";")  # SS=<number>
        # HACK: The position of the SS= changes, so we just search for it
        somatic_status = [item for item in somatic_status
                          if item.startswith("SS=")][0]
        somatic_status = int(somatic_status.split("=")[1])  # Get the number

        try:
            ndat[Ifreq] = str(float(ndat[Ifreq].rstrip("%")) / 100)
        except ValueError:  # illegal binary characters -- set frequency to zero
            ndat[Ifreq] = "0.0"
        try:
            tdat[Ifreq] = str(float(tdat[Ifreq].rstrip("%")) / 100)
        except ValueError:  # illegal binary characters -- set frequency to zero
            tdat[Ifreq] = "0.0"
        line[9] = ":".join(ndat)
        line[10] = ":".join(tdat)
        if somatic_status == 5:
            # "Unknown" states are broken in current versions of VarScan
            # so we just bail out here for now
            return

    #FIXME: VarScan also produces invalid REF records (e.g. CAA/A)
    # This is not handled yet.

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
    return "\t".join(line)


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
    sample_list = _create_sample_list(align_bams, out_file)
    mpileup = samtools.prep_mpileup(align_bams, ref_file, config, max_read_depth,
                                    target_regions=target_regions, want_bcf=False)
    # VarScan fails to generate a header on files that start with
    # zerocoverage calls; strip these with grep, we're not going to
    # call on them
    remove_zerocoverage = r"{ ifne grep -v -P '\t0\t\t$' || true; }"
    # we use ifne from moreutils to ensure we process only on files with input, skipping otherwise
    # http://manpages.ubuntu.com/manpages/natty/man1/ifne.1.html
    with tx_tmpdir(items[0]) as tmp_dir:
        jvm_opts = _get_varscan_opts(config, tmp_dir)
        fix_ambig_ref = vcfutils.fix_ambiguous_cl()
        fix_ambig_alt = vcfutils.fix_ambiguous_cl(5)
        py_cl = os.path.join(os.path.dirname(sys.executable), "py")
        cmd = ("{mpileup} | {remove_zerocoverage} | "
                "ifne varscan {jvm_opts} mpileup2cns --min-coverage 5 --p-value 0.98 "
                "  --vcf-sample-list {sample_list} --output-vcf --variants | "
               "{py_cl} -x 'bcbio.variation.varscan.fix_varscan_output(x)' | "
                "{fix_ambig_ref} | {fix_ambig_alt} | ifne vcfuniqalleles > {out_file}")
        do.run(cmd.format(**locals()), "Varscan", None,
                [do.file_exists(out_file)])
    os.remove(sample_list)
    # VarScan can create completely empty files in regions without
    # variants, so we create a correctly formatted empty file
    if os.path.getsize(out_file) == 0:
        write_empty_vcf(out_file)

    if orig_out_file.endswith(".gz"):
        vcfutils.bgzip_and_index(out_file, config)
