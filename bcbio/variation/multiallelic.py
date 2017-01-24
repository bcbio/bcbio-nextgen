"""Handle multi-allelic variants by splitting into single sample records.

Downstream tools like GEMINI don't handle multi-allelic variants. This
splits these into multiple single allele variants and re-annotates
effects.

More generally, multi-allelic variants are tricky to represent in
comparisons and storage.  There are useful discussions on-going in the
GA4GH about allele/count based approaches to handling this more
generally:

https://github.com/ga4gh/schemas/issues/169

There are multiple approaches to do the decomposition:

   - vt decompose: https://github.com/atks/vt
   - vcfbreakmulti: https://github.com/ekg/vcflib
   - bcftools norm -m '-any': http://samtools.github.io/bcftools/
   - vcf_parser --split: https://github.com/moonso/vcf_parser https://github.com/moonso/genmod

vt handles the most cases cleanly and correctly reconstructs PLs for
each genotype, removing other FORMAT items which are not changed so
the resulting VCF is still valid.
"""

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import effects, vcfutils

def to_single(in_file, data, passonly=False):
    """Convert multi-allelic inputs in the original VCF file into single alleles.
    """
    out_file = "%s-nomultiallelic%s" % utils.splitext_plus(in_file)
    if not utils.file_exists(out_file):
        if vcfutils.vcf_has_variants(in_file):
            ready_ma_file = _decompose(in_file, data, passonly=passonly)
            ann_ma_file, _ = effects.add_to_vcf(ready_ma_file, data)
            if ann_ma_file:
                ready_ma_file = ann_ma_file
            out_file = ready_ma_file
        else:
            utils.symlink_plus(in_file, out_file)
    return vcfutils.bgzip_and_index(out_file, data["config"])

def _decompose(in_file, data, passonly=False):
    """Convert multi-allelic variants into single allelic.

    vt normalize has the -n flag passed (skipping reference checks) because
    of errors where the reference genome has non GATCN ambiguous bases. These
    are not supported in VCF, so you'll have a mismatch of N in VCF versus R
    (or other ambiguous bases) in the genome.
    """
    out_file = "%s-decompose%s" % utils.splitext_plus(in_file)
    if not utils.file_exists(out_file):
        ref_file = dd.get_ref_file(data)
        assert out_file.endswith(".vcf.gz")
        sample_filter = "| bcftools view -f 'PASS,.' " if passonly else ""
        with file_transaction(data, out_file) as tx_out_file:
            cmd = ("gunzip -c %s | "
                   "sed 's/ID=AD,Number=./ID=AD,Number=R/' | "
                   "vt decompose -s - "
                   "| vt normalize -n -r %s - "
                   """| awk '{ gsub("./-65", "./."); print $0 }'"""
                   "%s"
                   "| sed -e 's/Number=A/Number=1/g' | bgzip -c > %s")
            do.run(cmd % (in_file, ref_file, sample_filter, tx_out_file), "Multi-allelic to single allele")
    return vcfutils.bgzip_and_index(out_file, data["config"])
