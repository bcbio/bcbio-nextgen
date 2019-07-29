import gffutils
import tempfile
import os
import random
import re

from bcbio import utils
from bcbio.utils import file_exists, open_gzipsafe
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.log import logger

def guess_infer_extent(gtf_file):
    """
    guess if we need to use the gene extent option when making a gffutils
    database by making a tiny database of 1000 lines from the original
    GTF and looking for all of the features
    """
    _, ext = os.path.splitext(gtf_file)
    tmp_out = tempfile.NamedTemporaryFile(suffix=".gtf", delete=False).name
    with open(tmp_out, "w") as out_handle:
        count = 0
        in_handle = utils.open_gzipsafe(gtf_file)
        for line in in_handle:
            if count > 1000:
                break
            out_handle.write(line)
            count += 1
        in_handle.close()
    db = gffutils.create_db(tmp_out, dbfn=":memory:", infer_gene_extent=False)
    os.remove(tmp_out)
    features = [x for x in db.featuretypes()]
    if "gene" in features and "transcript" in features:
        return False
    else:
        return True

def get_gtf_db(gtf, in_memory=False):
    """
    create a gffutils DB, in memory if we don't have write permissions
    """
    db_file = gtf + ".db"
    if file_exists(db_file):
        return gffutils.FeatureDB(db_file)
    if not os.access(os.path.dirname(db_file), os.W_OK | os.X_OK):
        in_memory = True
    db_file = ":memory:" if in_memory else db_file
    if in_memory or not file_exists(db_file):
        infer_extent = guess_infer_extent(gtf)
        disable_extent = not infer_extent
        db = gffutils.create_db(gtf, dbfn=db_file,
                                disable_infer_genes=disable_extent,
                                disable_infer_transcripts=disable_extent)
    if in_memory:
        return db
    else:
        return gffutils.FeatureDB(db_file)

def gtf_to_bed(gtf, alt_out_dir=None):
    """
    create a BED file of transcript-level features with attached gene name
    or gene ids
    """
    out_file = os.path.splitext(gtf)[0] + ".bed"
    if file_exists(out_file):
        return out_file
    if not os.access(os.path.dirname(out_file), os.W_OK | os.X_OK):
        if not alt_out_dir:
            raise IOError("Cannot write transcript BED output file %s" % out_file)
        else:
            out_file = os.path.join(utils.safe_makedir(alt_out_dir), os.path.basename(out_file))
    if file_exists(out_file):
        return out_file
    with open(out_file, "w") as out_handle:
        db = get_gtf_db(gtf)
        for feature in db.features_of_type('transcript', order_by=("seqid", "start", "end")):
            chrom = feature.chrom
            start = feature.start - 1
            end = feature.end
            attributes = feature.attributes.keys()
            strand = feature.strand
            name = (feature['gene_name'][0] if 'gene_name' in attributes else
                    feature['gene_id'][0])
            line = "\t".join([str(x) for x in [chrom, start, end, name, ".",
                                               strand]])
            out_handle.write(line + "\n")
    return out_file

def complete_features(db):
    """
    iterator returning features which are complete (have a 'gene_id' and a
    'transcript_id')
    """
    for feature in db.all_features():
        gene_id = feature.attributes.get('gene_id', [None])[0]
        transcript_id = feature.attributes.get('transcript_id', [None])[0]
        if gene_id and transcript_id and feature.featuretype != "transcript":
            yield feature

def gtf_to_fasta(gtf_file, ref_fasta, cds=False, out_file=None):
    """
    convert a GTF to FASTA format if cds=True, use the start/stop codons
    to output only the CDS
    handles malformed FASTA files where a single transcript is repeated multiple
    times by just using the first one
    """
    if out_file and file_exists(out_file):
        return out_file

    if not out_file:
        out_file = tempfile.NamedTemporaryFile(delete=False, suffix=".fa").name

    tmp_file = out_file + ".tmp"
    if cds:
        cmd = "gffread -g {ref_fasta} -x {tx_tmp_file} {gtf_file}"
    else:
        cmd = "gffread -g {ref_fasta} -w {tx_tmp_file} {gtf_file}"
    message = "Converting %s to FASTA format." % gtf_file
    with file_transaction(tmp_file) as tx_tmp_file:
        do.run(cmd.format(**locals()), message)

    transcript = ""
    skipping = False
    with file_transaction(out_file) as tx_out_file:
        with open(tmp_file) as in_handle, open(tx_out_file, "w") as out_handle:
            for line in in_handle:
                if line.startswith(">"):
                    cur_transcript = line.split(" ")[0][1:]
                    if transcript == cur_transcript:
                        logger.info("Transcript %s has already been seen, skipping this "
                                    "version." % cur_transcript)
                        skipping = True
                    else:
                        transcript = cur_transcript
                        skipping = False
                    line = ">" + transcript + "\n"
                if not skipping:
                    out_handle.write(line)
    os.remove(tmp_file)
    return out_file

def partition_gtf(gtf, coding=False, out_file=False):
    """
    return a GTF file of all non-coding or coding transcripts. the GTF must be annotated
    with gene_biotype = "protein_coding" or to have the source column set to the
    biotype for all coding transcripts. set coding to
    True to get only the coding, false to get only the non-coding
    """
    if out_file and file_exists(out_file):
        return out_file
    if not out_file:
        out_file = tempfile.NamedTemporaryFile(delete=False,
                                               suffix=".gtf").name

    if coding:
        pred = lambda biotype: biotype and biotype == "protein_coding"
    else:
        pred = lambda biotype: biotype and biotype != "protein_coding"

    biotype_lookup = _biotype_lookup_fn(gtf)

    db = get_gtf_db(gtf)
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            for feature in db.all_features():
                biotype = biotype_lookup(feature)
                if pred(biotype):
                    out_handle.write(str(feature) + "\n")
    return out_file

def split_gtf(gtf, sample_size=None, out_dir=None):
    """
    split a GTF file into two equal parts, randomly selecting genes.
    sample_size will select up to sample_size genes in total
    """
    if out_dir:
        part1_fn = os.path.basename(os.path.splitext(gtf)[0]) + ".part1.gtf"
        part2_fn = os.path.basename(os.path.splitext(gtf)[0]) + ".part2.gtf"
        part1 = os.path.join(out_dir, part1_fn)
        part2 = os.path.join(out_dir, part2_fn)
        if file_exists(part1) and file_exists(part2):
            return part1, part2
    else:
        part1 = tempfile.NamedTemporaryFile(delete=False, suffix=".part1.gtf").name
        part2 = tempfile.NamedTemporaryFile(delete=False, suffix=".part2.gtf").name

    db = get_gtf_db(gtf)
    gene_ids = set([x['gene_id'][0] for x in db.all_features()])
    if not sample_size or (sample_size and sample_size > len(gene_ids)):
        sample_size = len(gene_ids)
    gene_ids = set(random.sample(gene_ids, sample_size))
    part1_ids = set(random.sample(gene_ids, sample_size // 2))
    part2_ids = gene_ids.difference(part1_ids)
    with open(part1, "w") as part1_handle:
        for gene in part1_ids:
            for feature in db.children(gene):
                part1_handle.write(str(feature) + "\n")
    with open(part2, "w") as part2_handle:
        for gene in part2_ids:
            for feature in db.children(gene):
                part2_handle.write(str(feature) + "\n")
    return part1, part2

def get_coding_noncoding_transcript_ids(gtf):
    """
    return a set of coding and non-coding transcript_ids from a GTF
    """
    coding_gtf = partition_gtf(gtf, coding=True)
    coding_db = get_gtf_db(coding_gtf)
    coding_ids = set([x['transcript_id'][0] for x in coding_db.all_features()
                  if 'transcript_id' in x.attributes])
    noncoding_gtf = partition_gtf(gtf)
    noncoding_db = get_gtf_db(noncoding_gtf)
    noncoding_ids = set([x['transcript_id'][0] for x in noncoding_db.all_features()
                     if 'transcript_id' in x.attributes])
    return coding_ids, noncoding_ids

def get_gene_source_set(gtf):
    """
    get a dictionary of the set of all sources for a gene
    """
    gene_to_source = {}
    db = get_gtf_db(gtf)
    for feature in complete_features(db):
        gene_id = feature['gene_id'][0]
        sources = gene_to_source.get(gene_id, set([])).union(set([feature.source]))
        gene_to_source[gene_id] = sources
    return gene_to_source

def get_transcript_source_set(gtf):
    """
    get a dictionary of the set of all sources of the gene for a given
    transcript
    """
    gene_to_source = get_gene_source_set(gtf)
    transcript_to_source = {}
    db = get_gtf_db(gtf)
    for feature in complete_features(db):
        gene_id = feature['gene_id'][0]
        transcript_to_source[feature['transcript_id'][0]] = gene_to_source[gene_id]
    return transcript_to_source

def get_rRNA(gtf):
    """
    extract rRNA genes and transcripts from a gtf file
    """
    rRNA_biotypes = ["rRNA", "Mt_rRNA", "tRNA", "MT_tRNA"]
    features = set()
    with open_gzipsafe(gtf) as in_handle:
        for line in in_handle:
            if not "gene_id" in line or not "transcript_id" in line:
                continue
            if any(x in line for x in rRNA_biotypes):
                geneid = line.split("gene_id")[1].split(" ")[1]
                geneid = _strip_non_alphanumeric(geneid)
                geneid = _strip_feature_version(geneid)
                txid = line.split("transcript_id")[1].split(" ")[1]
                txid = _strip_non_alphanumeric(txid)
                txid = _strip_feature_version(txid)
                features.add((geneid, txid))
    return features

def _biotype_lookup_fn(gtf):
    """
    return a function that will look up the biotype of a feature
    this checks for either gene_biotype or biotype being set or for the source
    column to have biotype information
    """
    db = get_gtf_db(gtf)
    sources = set([feature.source for feature in db.all_features()])
    gene_biotypes = set([feature.attributes.get("gene_biotype", [None])[0]
                         for feature in db.all_features()])
    biotypes = set([feature.attributes.get("biotype", [None])[0]
                    for feature in db.all_features()])
    if "protein_coding" in sources:
        return lambda feature: feature.source
    elif "protein_coding" in biotypes:
        return lambda feature: feature.attributes.get("biotype", [None])[0]
    elif "protein_coding" in gene_biotypes:
        return lambda feature: feature.attributes.get("gene_biotype", [None])[0]
    else:
        return None

def tx2genedict(gtf, keep_version=False):
    """
    produce a tx2gene dictionary from a GTF file
    """
    d = {}
    with open_gzipsafe(gtf) as in_handle:
        for line in in_handle:
            if "gene_id" not in line or "transcript_id" not in line:
                continue
            geneid = line.split("gene_id")[1].split(" ")[1]
            geneid = _strip_non_alphanumeric(geneid)
            if not geneid:
                continue
            txid = line.split("transcript_id")[1].split(" ")[1]
            txid = _strip_non_alphanumeric(txid)
            if keep_version and "transcript_version" in line:
                txversion = line.split("transcript_version")[1].split(" ")[1]
                txversion = _strip_non_alphanumeric(txversion)
                txid  += "." + txversion
            if has_transcript_version(line) and not keep_version:
                txid = _strip_feature_version(txid)
                geneid = _strip_feature_version(geneid)
            d[txid] = geneid
    return d

def _strip_feature_version(featureid):
    """
    some feature versions are encoded as featureid.version, this strips those off, if they exist
    """
    version_detector = re.compile(r"(?P<featureid>.*)(?P<version>\.\d+)")
    match = version_detector.match(featureid)
    if match:
        return match.groupdict()["featureid"]
    else:
        return featureid

def _strip_non_alphanumeric(string):
    return string.replace('"', '').replace(';', '')

def has_transcript_version(line):
    version_detector = re.compile(r".*(?P<version>\.\d+)")
    if "transcript_version" in line:
        return True
    txid = line.split("transcript_id")[1].split(" ")[1]
    txid = _strip_non_alphanumeric(txid)
    if version_detector.match(txid):
        return True

def tx2genefile(gtf, out_file=None, data=None, tsv=True, keep_version=False):
    """
    write out a file of transcript->gene mappings.
    """
    if tsv:
        extension = ".tsv"
        sep = "\t"
    else:
        extension = ".csv"
        sep = ","
    if file_exists(out_file):
        return out_file
    with file_transaction(data, out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            for k, v in tx2genedict(gtf, keep_version).items():
                out_handle.write(sep.join([k, v]) + "\n")
    logger.info("tx2gene file %s created from %s." % (out_file, gtf))
    return out_file

def is_qualimap_compatible(gtf):
    """
    Qualimap needs a very specific GTF format or it fails, so skip it if
    the GTF is not in that format
    """
    if not gtf:
        return False
    db = get_gtf_db(gtf)
    def qualimap_compatible(feature):
        gene_id = feature.attributes.get('gene_id', [None])[0]
        transcript_id = feature.attributes.get('transcript_id', [None])[0]
        gene_biotype = feature.attributes.get('gene_biotype', [None])[0]
        return gene_id and transcript_id and gene_biotype
    for feature in db.all_features():
        if qualimap_compatible(feature):
            return True
    return False

def canonical_transcripts(gtf, out_file):
    """
    given a GTF file, produce a new GTF file with only the longest transcript
    for each gene
    function lifted from:
    https://pythonhosted.org/gffutils/_modules/gffutils/helpers.html
    """
    if file_exists(out_file):
        return out_file
    db = get_gtf_db(gtf)
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            for gene in db.features_of_type('gene'):
                exon_list = []
                for ti, transcript in enumerate(db.children(gene, level=1)):
                    cds_len = 0
                    total_len = 0
                    exons = list(db.children(transcript, level=1))
                    for exon in exons:
                        exon_length = len(exon)
                        if exon.featuretype == 'CDS':
                            cds_len += exon_length
                        total_len += exon_length

                    exon_list.append((cds_len, total_len, transcript, exons))

                # If we have CDS, then use the longest coding transcript
                if max(i[0] for i in exon_list) > 0:
                    best = sorted(exon_list)[0]
                # Otherwise, just choose the longest
                else:
                    best = sorted(exon_list, key=lambda x: x[1])[0]
                for exon in db.children(best[2], level=1):
                    out_handle.write(str(exon) + "\n")
    return out_file

def is_cpat_compatible(gtf):
    """
    CPAT needs some transcripts annotated with protein coding status to work
    properly
    """
    if not gtf:
        return False
    db = get_gtf_db(gtf)
    pred = lambda biotype: biotype and biotype == "protein_coding"
    biotype_lookup = _biotype_lookup_fn(gtf)
    if not biotype_lookup:
        return False
    db = get_gtf_db(gtf)
    for feature in db.all_features():
        biotype = biotype_lookup(feature)
        if pred(biotype):
            return True
    return False
