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

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import effects, vcfutils

def to_single(in_file, data):
    """Convert multi-allelic inputs in the original VCF file into single alleles.
    """
    out_file = "%s-nomultiallelic%s" % utils.splitext_plus(in_file)
    if not utils.file_exists(out_file):
        ba_file, ma_file = _split_mulitallelic(in_file, data)
        if vcfutils.vcf_has_variants(ma_file):
            ready_ma_file = _decompose(ma_file, data)
            ann_ma_file = effects.add_to_vcf(ready_ma_file, data)
            if ann_ma_file:
                ready_ma_file = ann_ma_file
            out_file = vcfutils.merge_sorted([ready_ma_file, ba_file], out_file, data)
        else:
            utils.symlink_plus(in_file, out_file)
    return vcfutils.bgzip_and_index(out_file, data["config"])

def _decompose(in_file, data):
    """Convert multi-allelic variants into single allelic.
    """
    out_file = "%s-decompose%s" % utils.splitext_plus(in_file)
    if not utils.file_exists(out_file):
        assert out_file.endswith(".vcf.gz")
        with file_transaction(data, out_file) as tx_out_file:
            cmd = ("gunzip -c %s | "
                   "sed 's/ID=AD,Number=./ID=AD,Number=R/' | "
                   "vt decompose -s - "
                   """| awk '{ gsub("./-65", "./."); print $0 }'"""
                   "| bgzip -c > %s")
            do.run(cmd % (in_file, tx_out_file), "Multi-allelic to single allele")
    return vcfutils.bgzip_and_index(out_file, data["config"])

def _split_mulitallelic(in_file, data):
    """Split input into biallelic and multiallelic files.
    """
    ba_out = "%s-biallelic%s" % utils.splitext_plus(in_file)
    ma_out = "%s-multiallelic%s" % utils.splitext_plus(in_file)
    for out_file, select_type in [(ba_out, "BIALLELIC"), (ma_out, "MULTIALLELIC")]:
        if not utils.file_exists(out_file):
            with file_transaction(data, out_file) as tx_out_file:
                params = ["-T", "SelectVariants", "-R", dd.get_ref_file(data),
                          "--variant", in_file, "--out", tx_out_file,
                          "-restrictAllelesTo", select_type]
                jvm_opts = broad.get_gatk_framework_opts(data["config"])
                cmd = ["gatk-framework"] + jvm_opts + params
                do.run(cmd, "Select %s variants" % select_type)
        vcfutils.bgzip_and_index(out_file, data["config"])
    return ba_out, ma_out
