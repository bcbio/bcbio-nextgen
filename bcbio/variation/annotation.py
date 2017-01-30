"""Annotated variant VCF files with additional information.

- GATK variant annotation with snpEff predicted effects.
"""
import glob
import gzip
import os
import sys

import pybedtools

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import vcfutils

def get_gatk_annotations(config, include_depth=True):
    """Retrieve annotations to use for GATK VariantAnnotator.

    If include_depth is false, we'll skip annotating DP. Since GATK downsamples
    this will undercount on high depth sequencing and the standard outputs
    from the original callers may be preferable.
    """
    broad_runner = broad.runner_from_config(config)
    anns = ["BaseQualityRankSumTest", "FisherStrand",
            "GCContent", "HaplotypeScore", "HomopolymerRun",
            "MappingQualityRankSumTest", "MappingQualityZero",
            "QualByDepth", "ReadPosRankSumTest", "RMSMappingQuality"]
    if include_depth:
        anns += ["DepthPerAlleleBySample"]
        if broad_runner.gatk_type() == "restricted":
            anns += ["Coverage"]
        else:
            anns += ["DepthOfCoverage"]
    return anns

def finalize_vcf(in_file, variantcaller, items):
    """Perform cleanup and annotation of the final VCF.
    """
    out_file = "%s-annotated%s" % utils.splitext_plus(in_file)
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            cl = _add_vcf_header_sample_cl(in_file, items, out_file)
            if cl:
                cmd = "{cl} | bgzip -c > {tx_out_file}"
                do.run(cmd.format(**locals()), "Annotate")
    if utils.file_exists(out_file):
        return vcfutils.bgzip_and_index(out_file, items[0]["config"])
    else:
        return in_file

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
            cmd = "bcftools reheader -h {new_header} {in_file} | bcftools view "
        # bcftools reheader does not work with empty VCF files as of samtools 1.3
        else:
            cmd = "cat {new_header}"
        return cmd.format(**locals())

def _update_header(orig_vcf, base_file, new_lines, chrom_process_fn=None):
    """Fix header with additional lines and remapping of generic sample names.
    """
    new_header = "%s-header.txt" % utils.splitext_plus(base_file)[0]
    with open(new_header, "w") as out_handle:
        chrom_line = None
        with gzip.open(orig_vcf) as in_handle:
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

def add_dbsnp(orig_file, dbsnp_file, data, out_file=None):
    """Annotate a VCF file with dbSNP.
    """
    orig_file = vcfutils.bgzip_and_index(orig_file, data["config"])
    if out_file is None:
        out_file = "%s-wdbsnp.vcf.gz" % utils.splitext_plus(orig_file)[0]
    if not utils.file_uptodate(out_file, orig_file):
        with file_transaction(data, out_file) as tx_out_file:
            conf_file = os.path.join(os.path.dirname(tx_out_file), "dbsnp.conf")
            with open(conf_file, "w") as out_handle:
                out_handle.write('[[annotation]]\n')
                out_handle.write('file="%s"\n' % os.path.normpath(os.path.join(dd.get_work_dir(data), dbsnp_file)))
                out_handle.write('fields=["ID"]\n')
                out_handle.write('names=["rs_ids"]\n')
                out_handle.write('ops=["concat"]\n')
            ref_file = dd.get_ref_file(data)
            cmd = ("vcfanno {conf_file} {orig_file} | "
                   "bcftools annotate --set-id +'%INFO/rs_ids' -o {tx_out_file} -O z")
            do.run(cmd.format(**locals()), "Annotate with dbSNP")
    return vcfutils.bgzip_and_index(out_file, data["config"])

def annotate_nongatk_vcf(orig_file, bam_files, dbsnp_file, ref_file, config):
    """Annotate a VCF file with dbSNP and standard GATK called annotations.
    """
    orig_file = vcfutils.bgzip_and_index(orig_file, config)
    broad_runner = broad.runner_from_config_safe(config)
    if not broad_runner or not broad_runner.has_gatk():
        return orig_file
    else:
        out_file = "%s-gatkann%s" % utils.splitext_plus(orig_file)
        if not utils.file_exists(out_file):
            with file_transaction(config, out_file) as tx_out_file:
                # Avoid issues with incorrectly created empty GATK index files.
                # Occurs when GATK cannot lock shared dbSNP database on previous run
                idx_file = orig_file + ".idx"
                if os.path.exists(idx_file) and not utils.file_exists(idx_file):
                    os.remove(idx_file)
                annotations = get_gatk_annotations(config, include_depth=False)
                params = ["-T", "VariantAnnotator",
                          "-R", ref_file,
                          "--variant", orig_file,
                          "--out", tx_out_file,
                          "-L", orig_file]
                if dbsnp_file:
                    params += ["--dbsnp", dbsnp_file]
                for bam_file in bam_files:
                    params += ["-I", bam_file]
                for x in annotations:
                    params += ["-A", x]
                if ("--allow_potentially_misencoded_quality_scores" not in params
                      and "-allowPotentiallyMisencodedQuals" not in params):
                    params += ["--allow_potentially_misencoded_quality_scores"]
                # be less stringent about BAM and VCF files (esp. N in CIGAR for RNA-seq)
                # start by removing existing -U or --unsafe opts
                # (if another option is added to Gatk that starts with -U... this may create a bug)
                unsafe_options = [x for x in params if x.startswith(("-U", "--unsafe"))]
                for my_opt in unsafe_options:
                    ind_to_rem = params.index(my_opt)
                    # are the options given as separate strings or in one?
                    if my_opt.strip() == "-U" or my_opt.strip() == "--unsafe":
                        params.pop(ind_to_rem + 1)
                    params.pop(ind_to_rem)
                params.extend(["-U", "ALL"])
                broad_runner = broad.runner_from_config(config)
                broad_runner.run_gatk(params)
        vcfutils.bgzip_and_index(out_file, config)
        return out_file

def get_context_files(data):
    """Retrieve pre-installed annotation files for annotating genome context.
    """
    ref_file = dd.get_ref_file(data)
    all_files = []
    for ext in [".bed.gz"]:
        all_files += glob.glob(os.path.normpath(os.path.join(os.path.dirname(ref_file), os.pardir,
                                                             "coverage", "problem_regions", "*",
                                                             "*%s" % ext)))
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
