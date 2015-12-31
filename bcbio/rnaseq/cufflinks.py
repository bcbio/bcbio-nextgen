"""Assess transcript abundance in RNA-seq experiments using Cufflinks.

http://cufflinks.cbcb.umd.edu/manual.html
"""
import os
import shutil
import tempfile
import pandas as pd
from bcbio.utils import get_in, file_exists, safe_makedir
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.rnaseq import gtf, annotate_gtf


def run(align_file, ref_file, data):
    config = data["config"]
    cmd = _get_general_options(align_file, config)
    cmd.extend(_get_no_assembly_options(ref_file, data))
    out_dir = _get_output_dir(align_file, data)
    tracking_file = os.path.join(out_dir, "genes.fpkm_tracking")
    fpkm_file = os.path.join(out_dir, data['rgnames']['sample']) + ".fpkm"
    tracking_file_isoform = os.path.join(out_dir, "isoforms.fpkm_tracking")
    fpkm_file_isoform = os.path.join(out_dir, data['rgnames']['sample']) + ".isoform.fpkm"
    if not file_exists(fpkm_file):
        with file_transaction(data, out_dir) as tmp_out_dir:
            safe_makedir(tmp_out_dir)
            cmd.extend(["--output-dir", tmp_out_dir])
            cmd.extend([align_file])
            cmd = map(str, cmd)
            do.run(cmd, "Cufflinks on %s." % (align_file))
        fpkm_file = gene_tracking_to_fpkm(tracking_file, fpkm_file)
        fpkm_file_isoform = gene_tracking_to_fpkm(tracking_file_isoform, fpkm_file_isoform)
    return out_dir, fpkm_file, fpkm_file_isoform


def gene_tracking_to_fpkm(tracking_file, out_file):
    """
    take a gene-level tracking file from Cufflinks and output a two column
    table with the first column as IDs and the second column as FPKM for the
    sample. combines FPKM from the same genes into one FPKM value to fix
    this bug: http://seqanswers.com/forums/showthread.php?t=5224&page=2
    """
    if file_exists(out_file):
        return out_file
    df = pd.io.parsers.read_table(tracking_file, sep="\t", header=0)
    df = df[['tracking_id', 'FPKM']]
    df = df.groupby(['tracking_id']).sum()
    df.to_csv(out_file, sep="\t", header=False, index_label=False)
    return out_file


def _get_general_options(align_file, config):
    options = []
    cufflinks = config_utils.get_program("cufflinks", config)
    options.extend([cufflinks])
    options.extend(["--num-threads", config["algorithm"].get("num_cores", 1)])
    options.extend(["--quiet"])
    options.extend(["--no-update-check"])
    options.extend(["--max-bundle-frags", 2000000])
    options.extend(_get_stranded_flag(config))
    return options


def _get_no_assembly_options(ref_file, data):
    options = []
    options.extend(["--frag-bias-correct", ref_file])
    options.extend(["--multi-read-correct"])
    options.extend(["--upper-quartile-norm"])
    gtf_file = data["genome_resources"]["rnaseq"].get("transcripts", "")
    if gtf_file:
        options.extend(["--GTF", gtf_file])
    mask_file = data["genome_resources"]["rnaseq"].get("transcripts_mask", "")
    if mask_file:
        options.extend(["--mask-file", mask_file])

    return options


def _get_stranded_flag(config):
    strand_flag = {"unstranded": "fr-unstranded",
                   "firststrand": "fr-firststrand",
                   "secondstrand": "fr-secondstrand"}
    stranded = get_in(config, ("algorithm", "strandedness"), "unstranded").lower()
    assert stranded in strand_flag, ("%s is not a valid strandedness value. "
                                     "Valid values are 'firststrand', "
                                     "'secondstrand' and 'unstranded" % (stranded))
    flag = strand_flag[stranded]
    return ["--library-type", flag]


def _get_output_dir(align_file, data, sample_dir=True):
    config = data["config"]
    name = data["rgnames"]["sample"] if sample_dir else ""
    return os.path.join(get_in(data, ("dirs", "work")), "cufflinks", name)


def assemble(bam_file, ref_file, num_cores, out_dir, data):
    out_dir = os.path.join(out_dir, data["rgnames"]["sample"])
    safe_makedir(out_dir)
    out_file = os.path.join(out_dir, "cufflinks-assembly.gtf")
    cufflinks_out_file = os.path.join(out_dir, "transcripts.gtf")
    library_type = " ".join(_get_stranded_flag(data["config"]))
    if file_exists(out_file):
        return out_file
    with file_transaction(data, out_dir) as tmp_out_dir:
        cmd = ("cufflinks --output-dir {tmp_out_dir} --num-threads {num_cores} "
               "--frag-bias-correct {ref_file} "
               "{library_type} --multi-read-correct --upper-quartile-norm {bam_file}")
        cmd = cmd.format(**locals())
        do.run(cmd, "Assembling transcripts with Cufflinks using %s." % bam_file)
    shutil.move(cufflinks_out_file, out_file)
    return out_file


def clean_assembly(gtf_file, clean=None, dirty=None):
    """
    clean the likely garbage transcripts from the GTF file including:
    1. any novel single-exon transcripts
    2. any features with an unknown strand
    """
    base, ext = os.path.splitext(gtf_file)
    db = gtf.get_gtf_db(gtf_file, in_memory=True)
    clean = clean if clean else base + ".clean" + ext
    dirty = dirty if dirty else base + ".dirty" + ext
    if file_exists(clean):
        return clean, dirty
    logger.info("Cleaning features with an unknown strand from the assembly.")
    with open(clean, "w") as clean_handle, open(dirty, "w") as dirty_handle:
        for gene in db.features_of_type('gene'):
            for transcript in db.children(gene, level=1):
                if is_likely_noise(db, transcript):
                    write_transcript(db, dirty_handle, transcript)
                else:
                    write_transcript(db, clean_handle, transcript)
    return clean, dirty


def write_transcript(db, handle, transcript):
    for feature in db.children(transcript):
        handle.write(str(feature) + "\n")


def is_likely_noise(db, transcript):
    if is_novel_single_exon(db, transcript):
        return True
    if strand_unknown(db, transcript):
        return True


def strand_unknown(db, transcript):
    """
    for unstranded data with novel transcripts single exon genes
    will have no strand information. single exon novel genes are also
    a source of noise in the Cufflinks assembly so this removes them
    """
    features = list(db.children(transcript))
    strand = features[0].strand
    if strand == ".":
        return True
    else:
        return False


def is_novel_single_exon(db, transcript):
    features = list(db.children(transcript))
    exons = [x for x in features if x.featuretype == "exon"]
    class_code = features[0].attributes.get("class_code", None)[0]
    if len(exons) == 1 and class_code == "u":
        return True
    return False


def fix_cufflinks_attributes(ref_gtf, merged_gtf, data, out_file=None):
    """
    replace the cufflinks gene_id and transcript_id with the
    gene_id and transcript_id from ref_gtf, where available

    """
    base, ext = os.path.splitext(merged_gtf)
    fixed = out_file if out_file else base + ".clean.fixed" + ext
    if file_exists(fixed):
        return fixed
    ref_db = gtf.get_gtf_db(ref_gtf)
    merged_db = gtf.get_gtf_db(merged_gtf, in_memory=True)

    ref_tid_to_gid = {}
    for gene in ref_db.features_of_type('gene'):
        for transcript in ref_db.children(gene, level=1):
            ref_tid_to_gid[transcript.id] = gene.id

    ctid_to_cgid = {}
    ctid_to_oid = {}
    for gene in merged_db.features_of_type('gene'):
        for transcript in merged_db.children(gene, level=1):
            ctid_to_cgid[transcript.id] = gene.id
            feature = list(merged_db.children(transcript))[0]
            oid = feature.attributes.get("oId", [None])[0]
            if oid:
                ctid_to_oid[transcript.id] = oid
    cgid_to_gid = {}
    for ctid, oid in ctid_to_oid.items():
        cgid = ctid_to_cgid.get(ctid, None)
        oid = ctid_to_oid.get(ctid, None)
        gid = ref_tid_to_gid.get(oid, None) if oid else None
        if cgid and gid:
            cgid_to_gid[cgid] = gid

    with file_transaction(data, fixed) as tmp_fixed_file:
        with open(tmp_fixed_file, "w") as out_handle:
            for gene in merged_db.features_of_type('gene'):
                for transcript in merged_db.children(gene, level=1):
                    for feature in merged_db.children(transcript):
                        cgid = feature.attributes.get("gene_id", [None])[0]
                        gid = cgid_to_gid.get(cgid, None)
                        ctid = None
                        if gid:
                            feature.attributes["gene_id"][0] = gid
                            ctid = feature.attributes.get("transcript_id",
                                                          [None])[0]
                        tid = ctid_to_oid.get(ctid, None)
                        if tid:
                            feature.attributes["transcript_id"][0] = tid
                        if "nearest_ref" in feature.attributes:
                            del feature.attributes["nearest_ref"]
                        if "oId" in feature.attributes:
                            del feature.attributes["oId"]
                        out_handle.write(str(feature) + "\n")
    return fixed


def merge(assembled_gtfs, ref_file, gtf_file, num_cores, data):
    """
    run cuffmerge on a set of assembled GTF files
    """
    assembled_file = tempfile.NamedTemporaryFile(delete=False).name
    with open(assembled_file, "w") as temp_handle:
        for assembled in assembled_gtfs:
            temp_handle.write(assembled + "\n")
    out_dir = os.path.join("assembly", "cuffmerge")
    merged_file = os.path.join(out_dir, "merged.gtf")
    out_file = os.path.join(out_dir, "assembled.gtf")
    if file_exists(out_file):
        return out_file
    if not file_exists(merged_file):
        with file_transaction(data, out_dir) as tmp_out_dir:
            cmd = ("cuffmerge -o {tmp_out_dir} --ref-gtf {gtf_file} "
                   "--num-threads {num_cores} --ref-sequence {ref_file} "
                   "{assembled_file}")
            cmd = cmd.format(**locals())
            message = ("Merging the following transcript assemblies with "
                       "Cuffmerge: %s" % ", ".join(assembled_gtfs))
            do.run(cmd, message)
    clean, _ = clean_assembly(merged_file)
    fixed = fix_cufflinks_attributes(gtf_file, clean, data)
    classified = annotate_gtf.annotate_novel_coding(fixed, gtf_file, ref_file,
                                                    data)
    filtered = annotate_gtf.cleanup_transcripts(classified, gtf_file, ref_file)
    shutil.move(filtered, out_file)
    return out_file
