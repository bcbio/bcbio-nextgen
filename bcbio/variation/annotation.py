"""Annotated variant VCF files with additional information.

- GATK variant annotation with snpEff predicted effects.
"""
import glob
import os

import pybedtools
import toolz as tz

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import vcfutils

def get_gatk_annotations(config, include_depth=True, include_baseqranksum=True,
                         gatk_input=True):
    """Retrieve annotations to use for GATK VariantAnnotator.

    If include_depth is false, we'll skip annotating DP. Since GATK downsamples
    this will undercount on high depth sequencing and the standard outputs
    from the original callers may be preferable.

    BaseQRankSum can cause issues with some MuTect2 and other runs, so we
    provide option to skip it.
    """
    broad_runner = broad.runner_from_config(config)
    anns = ["MappingQualityRankSumTest", "MappingQualityZero",
            "QualByDepth", "ReadPosRankSumTest", "RMSMappingQuality"]
    if include_baseqranksum:
        anns += ["BaseQualityRankSumTest"]
    # Some annotations not working correctly with external datasets and GATK 3
    if gatk_input or broad_runner.gatk_type() == "gatk4":
        anns += ["FisherStrand"]
    if broad_runner.gatk_type() == "gatk4":
        anns += ["MappingQuality"]
    else:
        anns += ["GCContent", "HaplotypeScore", "HomopolymerRun"]
    if include_depth:
        anns += ["DepthPerAlleleBySample"]
        if broad_runner.gatk_type() in ["restricted", "gatk4"]:
            anns += ["Coverage"]
        else:
            anns += ["DepthOfCoverage"]
    return anns

def finalize_vcf(in_file, variantcaller, items):
    """Perform cleanup and dbSNP annotation of the final VCF.

    - Adds contigs to header for bcftools compatibility
    - adds sample information for tumor/normal
    """
    out_file = "%s-annotated%s" % utils.splitext_plus(in_file)
    if not utils.file_uptodate(out_file, in_file):
        header_cl = _add_vcf_header_sample_cl(in_file, items, out_file)
        contig_cl = _add_contig_cl(in_file, items, out_file)
        cls = [x for x in (contig_cl, header_cl) if x]
        if cls:
            post_cl = " | ".join(cls) + " | "
        else:
            post_cl = None
        dbsnp_file = tz.get_in(("genome_resources", "variation", "dbsnp"), items[0])
        if dbsnp_file:
            out_file = _add_dbsnp(in_file, dbsnp_file, items[0], out_file, post_cl)
    if utils.file_exists(out_file):
        return vcfutils.bgzip_and_index(out_file, items[0]["config"])
    else:
        return in_file

def _add_contig_cl(in_file, items, out_file):
    has_contigs = False
    with utils.open_gzipsafe(in_file) as in_handle:
        for line in in_handle:
            if line.startswith("##contig"):
                has_contigs = True
                break
            elif not line.startswith("##"):
                break
    if not has_contigs:
        return vcfutils.add_contig_to_header_cl(dd.get_ref_file(items[0]), out_file)

def _fix_generic_tn_names(paired):
    """Convert TUMOR/NORMAL names in output into sample IDs.
    """
    def run(line):
        parts = line.rstrip("\n\r").split("\t")
        if "TUMOR" in parts:
            parts[parts.index("TUMOR")] = paired.tumor_name
        if "TUMOUR" in parts:
            parts[parts.index("TUMOUR")] = paired.tumor_name
        if "NORMAL" in parts:
            assert paired.normal_name
            parts[parts.index("NORMAL")] = paired.normal_name
        return "\t".join(parts) + "\n"
    return run

def _add_vcf_header_sample_cl(in_file, items, base_file):
    """Add phenotype information to a VCF header.

    Encode tumor/normal relationships in VCF header.
    Could also eventually handle more complicated pedigree information if useful.
    """
    paired = vcfutils.get_paired(items)
    if paired:
        toadd = ["##SAMPLE=<ID=%s,Genomes=Tumor>" % paired.tumor_name]
        if paired.normal_name:
            toadd.append("##SAMPLE=<ID=%s,Genomes=Germline>" % paired.normal_name)
            toadd.append("##PEDIGREE=<Derived=%s,Original=%s>" % (paired.tumor_name, paired.normal_name))
        new_header = _update_header(in_file, base_file, toadd, _fix_generic_tn_names(paired))
        if vcfutils.vcf_has_variants(in_file):
            cmd = "bcftools reheader -h {new_header} | bcftools view "
            return cmd.format(**locals())

def _update_header(orig_vcf, base_file, new_lines, chrom_process_fn=None):
    """Fix header with additional lines and remapping of generic sample names.
    """
    new_header = "%s-sample_header.txt" % utils.splitext_plus(base_file)[0]
    with open(new_header, "w") as out_handle:
        chrom_line = None
        with utils.open_gzipsafe(orig_vcf) as in_handle:
            for line in in_handle:
                if line.startswith("##"):
                    out_handle.write(line)
                else:
                    chrom_line = line
                    break
        assert chrom_line is not None
        for line in new_lines:
            out_handle.write(line + "\n")
        if chrom_process_fn:
            chrom_line = chrom_process_fn(chrom_line)
        out_handle.write(chrom_line)
    return new_header

_DBSNP_TEMPLATE = """
[[annotation]]
file="%s"
fields=["ID"]
names=["rs_ids"]
ops=["concat"]

[[postannotation]]
name="ID"
fields=["rs_ids"]
op="setid"
type="String"

[[postannotation]]
fields=["rs_ids"]
op="delete"
"""

def _add_dbsnp(orig_file, dbsnp_file, data, out_file=None, post_cl=None):
    """Annotate a VCF file with dbSNP.

    vcfanno has flexible matching for NON_REF gVCF positions, matching
    at position and REF allele, matching ALT NON_REF as a wildcard.
    """
    orig_file = vcfutils.bgzip_and_index(orig_file, data["config"])
    if out_file is None:
        out_file = "%s-wdbsnp.vcf.gz" % utils.splitext_plus(orig_file)[0]
    if not utils.file_uptodate(out_file, orig_file):
        with file_transaction(data, out_file) as tx_out_file:
            conf_file = os.path.join(os.path.dirname(out_file), "dbsnp.conf")
            with open(conf_file, "w") as out_handle:
                out_handle.write(_DBSNP_TEMPLATE % os.path.normpath(os.path.join(dd.get_work_dir(data), dbsnp_file)))
            if not post_cl: post_cl = ""
            cores = dd.get_num_cores(data)
            cmd = ("vcfanno -p {cores} {conf_file} {orig_file} | {post_cl} "
                   "bgzip -c > {tx_out_file}")
            do.run(cmd.format(**locals()), "Annotate with dbSNP")
    return vcfutils.bgzip_and_index(out_file, data["config"])

def get_context_files(data):
    """Retrieve pre-installed annotation files for annotating genome context.
    """
    ref_file = dd.get_ref_file(data)
    all_files = []
    for ext in [".bed.gz"]:
        all_files += sorted(glob.glob(os.path.normpath(os.path.join(os.path.dirname(ref_file), os.pardir,
                                                                    "coverage", "problem_regions", "*",
                                                                    "*%s" % ext))))
    return sorted(all_files)

def add_genome_context(orig_file, data):
    """Annotate a file with annotations of genome context using vcfanno.
    """
    out_file = "%s-context.vcf.gz" % utils.splitext_plus(orig_file)[0]
    if not utils.file_uptodate(out_file, orig_file):
        with file_transaction(data, out_file) as tx_out_file:
            config_file = "%s.toml" % (utils.splitext_plus(tx_out_file)[0])
            with open(config_file, "w") as out_handle:
                all_names = []
                for fname in dd.get_genome_context_files(data):
                    bt = pybedtools.BedTool(fname)
                    if bt.field_count() >= 4:
                        d, base = os.path.split(fname)
                        _, prefix = os.path.split(d)
                        name = "%s_%s" % (prefix, utils.splitext_plus(base)[0])
                        out_handle.write("[[annotation]]\n")
                        out_handle.write('file = "%s"\n' % fname)
                        out_handle.write("columns = [4]\n")
                        out_handle.write('names = ["%s"]\n' % name)
                        out_handle.write('ops = ["uniq"]\n')
                        all_names.append(name)
                out_handle.write("[[postannotation]]\n")
                out_handle.write("fields = [%s]\n" % (", ".join(['"%s"' % n for n in all_names])))
                out_handle.write('name = "genome_context"\n')
                out_handle.write('op = "concat"\n')
                out_handle.write('type = "String"\n')
            cmd = "vcfanno {config_file} {orig_file} | bgzip -c > {tx_out_file}"
            do.run(cmd.format(**locals()), "Annotate with problem annotations", data)
    return vcfutils.bgzip_and_index(out_file, data["config"])
