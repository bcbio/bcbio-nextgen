"""Normalization of VCF files for GEMINI and ensemble calling.

1. Split multi-allelic variants into single sample records.

   Downstream tools like GEMINI don't handle multi-allelic variants.
   `normalize` splits these into multiple single allele variants.

   More generally, multi-allelic variants are tricky to represent in
   comparisons and storage. There are useful discussions on-going in the
   GA4GH about allele/count based approaches to handling this more
   generally:

   https://github.com/ga4gh/schemas/issues/169

   There are multiple approaches to do the decomposition:

      - vt decompose: https://github.com/atks/vt
      - vcfbreakmulti: https://github.com/ekg/vcflib
      - bcftools norm -m '-any': https://samtools.github.io/bcftools/bcftools.html
      - vcf_parser --split: https://github.com/moonso/vcf_parser https://github.com/moonso/genmod

   vt handles the most cases cleanly and correctly reconstructs PLs for
   each genotype, removing other FORMAT items which are not changed so
   the resulting VCF is still valid.

2. Decompose biallelic block substitutions

   As part of normalization, decompose MNPs into SNPs, e.g. AG>CT -> A>C and G>T

   There are a few approaches:

     - vcfallelicprimitives -t DECOMPOSED --keep-geno: https://github.com/vcflib/vcflib#vcfallelicprimitives
     - vt decompose_blocksub: https://genome.sph.umich.edu/wiki/Vt#Decompose_biallelic_block_substitutions
          with `-a`, align/aggressive mode is turned on, which runs 3-4 longer but resolves
          additional clumped SNPs

   We are using `vcfallelicprimitives` in this module.

3. Normalizing indels

   Left-align and normalize indels, check if REF alleles match the reference.

     - bcftools norm: https://samtools.github.io/bcftools/bcftools.html
     - vt normalize: https://github.com/atks/vt

   We are using `vt normalize` in this module.

"""
import os

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import effects, vcfutils

cyvcf2 = utils.LazyImport("cyvcf2")

def normalize(in_file, data, passonly=False, normalize_indels=True, split_biallelic=True,
              rerun_effects=True, remove_oldeffects=False, nonrefonly=False, work_dir=None):
    """Normalizes variants and reruns SnpEFF for resulting VCF
    """
    if remove_oldeffects:
        out_file = "%s-noeff-nomultiallelic%s" % utils.splitext_plus(in_file)
    else:
        out_file = "%s-nomultiallelic%s" % utils.splitext_plus(in_file)
    if work_dir:
        out_file = os.path.join(work_dir, os.path.basename(out_file))
    if not utils.file_exists(out_file):
        if vcfutils.vcf_has_variants(in_file):
            ready_ma_file = _normalize(in_file, data, passonly=passonly,
                                       normalize_indels=normalize_indels,
                                       split_biallelic=split_biallelic,
                                       remove_oldeffects=remove_oldeffects,
                                       nonrefonly=nonrefonly,
                                       work_dir=work_dir)
            if rerun_effects:
                ann_ma_file, _ = effects.add_to_vcf(ready_ma_file, data)
                if ann_ma_file:
                    ready_ma_file = ann_ma_file
            utils.symlink_plus(ready_ma_file, out_file)
        else:
            utils.symlink_plus(in_file, out_file)
    return vcfutils.bgzip_and_index(out_file, data["config"])

def _normalize(in_file, data, passonly=False, normalize_indels=True, split_biallelic=True,
               remove_oldeffects=False, nonrefonly=False, work_dir=None):
    """Convert multi-allelic variants into single allelic.

    `vt normalize` has the -n flag passed (skipping reference checks) because
    of errors where the reference genome has non GATCN ambiguous bases. These
    are not supported in VCF, so you'll have a mismatch of N in VCF versus R
    (or other ambiguous bases) in the genome.
    """
    if remove_oldeffects:
        out_file = "%s-noeff-decompose%s" % utils.splitext_plus(in_file)
        old_effects = [a for a in ["CSQ", "ANN"] if a in cyvcf2.VCF(in_file)]
        if old_effects:
            clean_effects_cmd = " | bcftools annotate -x %s " % (",".join(["INFO/%s" % x for x in old_effects]))
        else:
            clean_effects_cmd = ""
    else:
        clean_effects_cmd = ""
        out_file = "%s-decompose%s" % utils.splitext_plus(in_file)
    if passonly or nonrefonly:
        subset_vcf_cmd = " | bcftools view "
        if passonly:
            subset_vcf_cmd += "-f 'PASS,.' "
        if nonrefonly:
            subset_vcf_cmd += "--min-ac 1:nref "
    else:
        subset_vcf_cmd = ""
    if work_dir:
        out_file = os.path.join(work_dir, os.path.basename(out_file))
    if not utils.file_exists(out_file):
        ref_file = dd.get_ref_file(data)
        assert out_file.endswith(".vcf.gz")
        with file_transaction(data, out_file) as tx_out_file:
            cmd = ("gunzip -c " + in_file +
                   subset_vcf_cmd + clean_effects_cmd +
                   (" | vcfallelicprimitives -t DECOMPOSED --keep-geno" if split_biallelic else "") +
                   " | sed 's/ID=AD,Number=./ID=AD,Number=R/'" +
                   " | vt decompose -s - " +
                   ((" | vt normalize -n -r " + ref_file + " - ") if normalize_indels else "") +
                   " | awk '{ gsub(\"./-65\", \"./.\"); print $0 }'" +
                   " | sed -e 's/Number=A/Number=1/g'" +
                   " | bgzip -c > " + tx_out_file
                   )
            do.run(cmd, "Multi-allelic to single allele")
    return vcfutils.bgzip_and_index(out_file, data["config"])
