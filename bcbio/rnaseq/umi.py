"""
Unique Molecular Identifier (UMI) handling
Most of this either uses Valentine Svennson's umis repository or adapts
code written from it.
https://github.com/vals/umis
"""
import pandas as pd
import scipy.io
import os
import copy
import glob
import sys
import subprocess
from itertools import repeat, islice
from distutils.version import LooseVersion

import bcbio.pipeline.datadict as dd
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio import utils
from bcbio.utils import (file_exists, safe_makedir)
from bcbio.distributed.transaction import file_transaction
from bcbio.bam.fastq import open_fastq
from bcbio.log import logger
from bcbio.rnaseq import gtf

import six


class SparseMatrix(object):

    def __init__(self, matrix=None, rownames=None, colnames=None):
        self.matrix = matrix
        self.rownames = rownames
        self.colnames = colnames

    def __repr__(self):
        return "%d x %d matrix of class %s" %(len(self.rownames),
                                              len(self.colnames),
                                              type(self.matrix))

    def read(self, filename, rowprefix=None, colprefix=None, delim=":"):
        """read a sparse matrix, loading row and column name files. if
        specified, will add a prefix to the row or column names"""

        self.matrix = scipy.io.mmread(filename)
        with open(filename + ".rownames") as in_handle:
            self.rownames = [x.strip() for x in in_handle]
            if rowprefix:
                self.rownames = [rowprefix + delim + x for x in self.rownames]
        with open(filename + ".colnames") as in_handle:
            self.colnames = [x.strip() for x in in_handle]
            if colprefix:
                self.colnames = [colprefix + delim + x for x in self.colnames]

    def write(self, filename):
        """read a sparse matrix, loading row and column name files"""
        if file_exists(filename):
            return filename
        out_files = [filename, filename + ".rownames", filename + ".colnames"]
        with file_transaction(out_files) as tx_out_files:
            with open(tx_out_files[0], "wb") as out_handle:
                scipy.io.mmwrite(out_handle, scipy.sparse.csr_matrix(self.matrix))
            pd.Series(self.rownames).to_csv(tx_out_files[1], index=False)
            pd.Series(self.colnames).to_csv(tx_out_files[2], index=False)
        return filename

    def cat(self, newsm, byrow=False):
        if self.matrix is None:
            self.matrix = newsm.matrix
            self.colnames = newsm.colnames
            self.rownames = newsm.rownames
        else:
            catter = scipy.sparse.vstack if byrow else scipy.sparse.hstack
            self.matrix = catter((self.matrix, newsm.matrix))
            if byrow:
                self.rownames = self.rownames + newsm.rownames
            else:
                self.colnames = self.colnames + newsm.colnames

TRANSFORM_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "data",
                             "umis")
TRANSFORM_FILES = glob.glob(os.path.join(TRANSFORM_DIR, "*-transform.json"))
SUPPORTED_TRANSFORMS = set([os.path.basename(x).replace("-transform.json", "")
                            for x in TRANSFORM_FILES])

def is_supported_transform(data):
    return dd.get_umi_type(data) in SUPPORTED_TRANSFORMS

def get_transform_file(stem):
    transform_file = os.path.join(TRANSFORM_DIR, stem + "-transform.json")
    return transform_file

def get_cellular_barcodes(data):
    if dd.get_cellular_barcodes(data):
        return dd.get_cellular_barcodes(data)
    if is_supported_transform(data):
        stem = dd.get_umi_type(data)
        bc1 = os.path.join(TRANSFORM_DIR, stem + "-cb1.txt.gz")
        bc2 = os.path.join(TRANSFORM_DIR, stem + "-cb2.txt.gz")
        bc3 = os.path.join(TRANSFORM_DIR, stem + "-cb3.txt.gz")
        return list(filter(file_exists, [bc1, bc2, bc3]))
    else:
        return []

def get_sample_barcodes(fn, out_dir):
    if not fn:
        logger.error("Sample demultiplexing needs a list of known indexes provided "
                     "with via the sample_barcodes option in the algorithm section.")
        sys.exit(1)
    utils.safe_makedir(out_dir)
    out_fn = os.path.join(out_dir, "barcodes.csv")
    with open(fn) as inh:
        with open(out_fn, 'w') as outh:
            for line in inh:
                outh.write("%s\n" % (line.strip().split(",")[0]))
    return out_fn

def _umis_cmd(data):
    """Return umis command line argument, with correct python and locale.
    """
    return "%s %s %s" % (utils.locale_export(),
                         utils.get_program_python("umis"),
                         config_utils.get_program("umis", data["config"], default="umis"))

def umi_transform(data):
    """
    transform each read by identifying the barcode and UMI for each read
    and putting the information in the read name
    """
    fqfiles = data["files"]
    fqfiles.extend(list(repeat("", 4-len(fqfiles))))
    fq1, fq2, fq3, fq4 = fqfiles
    umi_dir = os.path.join(dd.get_work_dir(data), "umis")
    safe_makedir(umi_dir)
    transform = dd.get_umi_type(data)
    if not transform:
        logger.info("No UMI transform specified, assuming pre-transformed data.")
        if is_transformed(fq1):
            logger.info("%s detected as pre-transformed, passing it on unchanged." % fq1)
            data["files"] = [fq1]
            return [[data]]
        else:
            logger.error("No UMI transform was specified, but %s does not look "
                         "pre-transformed." % fq1)
            sys.exit(1)

    if file_exists(transform):
        transform_file = transform
    else:
        transform_file = get_transform_file(transform)
        if not file_exists(transform_file):
            logger.error(
                "The UMI transform can be specified as either a file or a "
                "bcbio-supported transform. Either the file %s does not exist "
                "or the transform is not supported by bcbio. Supported "
                "transforms are %s."
                %(dd.get_umi_type(data), ", ".join(SUPPORTED_TRANSFORMS)))
            sys.exit(1)
    out_base = dd.get_sample_name(data) + ".umitransformed.fq.gz"
    out_file = os.path.join(umi_dir, out_base)
    if file_exists(out_file):
        data["files"] = [out_file]
        return [[data]]
    cellular_barcodes = get_cellular_barcodes(data)
    if len(cellular_barcodes) > 1:
        split_option = "--separate_cb"
    else:
        split_option = ""
    if dd.get_demultiplexed(data):
        demuxed_option = "--demuxed_cb %s" % dd.get_sample_name(data)
        split_option = ""
    else:
        demuxed_option = ""
    cores = dd.get_num_cores(data)
    # skip transformation if the file already looks transformed
    with open_fastq(fq1) as in_handle:
        read = next(in_handle)
        if "UMI_" in read:
            data["files"] = [out_file]
            return [[data]]
    locale_export = utils.locale_export()
    umis = _umis_cmd(data)
    cmd = ("{umis} fastqtransform {split_option} {transform_file} "
           "--cores {cores} {demuxed_option} "
           "{fq1} {fq2} {fq3} {fq4}"
           "| seqtk seq -L 20 - | gzip > {tx_out_file}")
    message = ("Inserting UMI and barcode information into the read name of %s"
               % fq1)
    with file_transaction(out_file) as tx_out_file:
        do.run(cmd.format(**locals()), message)
    data["files"] = [out_file]
    return [[data]]

def filter_barcodes(data):
    # if data was pre-demultiplexed, there is no need to filter the barcodes
    if dd.get_demultiplexed(data):
        return [[data]]
    fq1 = dd.get_input_sequence_files(data)[0]
    umi_dir = os.path.join(dd.get_work_dir(data), "umis")
    correction = dd.get_cellular_barcode_correction(data)
    bc = get_cellular_barcodes(data)
    if not bc:
        logger.info("No cellular barcodes found, skipping filtering.")
        return [[data]]
    bc1 = None
    bc2 = None
    bc3 = None
    umi_dir = os.path.join(dd.get_work_dir(data), "umis")
    if isinstance(bc, six.string_types):
        bc1 = bc
    if len(bc) == 1:
        bc1 = bc[0]
    if len(bc) > 1:
        bc1 = bc[0]
        bc2 = bc[1]
    if len(bc) == 3:
        bc3 = bc[2]
    out_base = dd.get_sample_name(data) + ".filtered.fq.gz"
    out_file = os.path.join(umi_dir, out_base)
    if file_exists(out_file):
        data["files"] = [out_file]
        return [[data]]

    ncores = dd.get_num_cores(data)
    umis = _umis_cmd(data)
    cmd = "{umis} cb_filter --cores {ncores} "
    if bc1:
        cmd += "--bc1 {bc1} "
        if correction:
            cmd += "--nedit {correction} "
    if bc2:
        cmd += "--bc2 {bc2} "
    if bc3:
        cmd += "--bc3 {bc3} "

    fq1_cmd = "{fq1} "
    fq1_cmd = fq1_cmd.format(fq1=fq1)
    cmd += "{fq1_cmd} | gzip -c > {tx_out_file}"

    sample_dir = os.path.join(umi_dir, dd.get_sample_name(data))
    safe_makedir(sample_dir)
    with file_transaction(out_file) as tx_out_file:
        message = "Filtering by cellular barcode."
        do.run(cmd.format(**locals()), message)
    data["files"] = [out_file]
    return [[data]]

def barcode_histogram(data):
    fq1 = dd.get_input_sequence_files(data)[0]
    umi_dir = os.path.join(dd.get_work_dir(data), "umis")
    sample_dir = os.path.join(umi_dir, dd.get_sample_name(data))
    safe_makedir(sample_dir)
    out_file = os.path.join(sample_dir, "cb-histogram.txt")
    filtered_out_file = os.path.join(sample_dir, "cb-histogram-filtered.txt")
    fq1_cmd = fq1
    umis = _umis_cmd(data)
    cmd = "{umis} cb_histogram {fq1_cmd} > {tx_out_file}"
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            message = "Computing cellular barcode counts for %s." % fq1
            do.run(cmd.format(**locals()), message)
    cutoff = dd.get_minimum_barcode_depth(data)
    filter_barcode_histogram(filtered_out_file, out_file, cutoff)
    newdata = dd.set_histogram_counts(data, filtered_out_file)
    return [[newdata]]

def tagcount(data):
    bam = dd.get_transcriptome_bam(data)
    umi_dir = os.path.join(dd.get_work_dir(data), "umis")
    sample_dir = os.path.join(umi_dir, dd.get_sample_name(data))
    out_prefix = os.path.join(sample_dir, dd.get_sample_name(data))
    out_file = out_prefix + ".mtx"
    if file_exists(out_file):
        data = dd.set_count_file(data, out_file)
        return [[data]]
    safe_makedir(sample_dir)
    cutoff = dd.get_minimum_barcode_depth(data)
    cb_histogram = os.path.join(sample_dir, "cb-histogram.txt")
    positional = "--positional" if dd.get_positional_umi(data, False) else ""
    if use_installed_transcriptome(data):
        gtf_file = dd.get_gtf_file(data)
    else:
        gtf_file  = dd.get_transcriptome_gtf(data, None)

    if gtf_file:
        gene_map_file = os.path.join(dd.get_work_dir(data), "annotation",
                                     os.path.basename(os.path.splitext(gtf_file)[0]) + "-tx2gene.tsv")
        gene_map_file = gtf.tx2genefile(gtf_file, gene_map_file, tsv=True)
        gene_map_flag = " --genemap {0} ".format(gene_map_file)
    else:
        gene_map_flag = ""

    message = "Counting alignments of transcripts in %s." % bam
    umis = _umis_cmd(data)
    cmd = ("{umis} fasttagcount --cb_cutoff {cutoff} "
           "{gene_map_flag} "
           "{positional} "
           "--cb_histogram {cb_histogram}")
    out_files = [out_file, out_file + ".rownames", out_file + ".colnames"]
    umi_matrix_file = out_prefix + "-dupes.mtx"
    out_files += [umi_matrix_file, umi_matrix_file + ".rownames",
                  umi_matrix_file + ".colnames"]
    if has_umi_matrix(data):
        umi_matrix_flag = " --umi_matrix {tx_umi_matrix_full} "
    else:
        umi_matrix_flag = ""
    cmd += umi_matrix_flag
    cmd += " {bam} {tx_out_file_full}"
    with file_transaction(out_files) as tx_out_files:
        tx_out_file = tx_out_files[0]
        tx_out_file_full = tx_out_file + ".full"
        tx_umi_matrix = tx_out_files[3]
        tx_umi_matrix_full = tx_out_files[3] + ".full"
        do.run(cmd.format(**locals()), message)
        cmd = ("{umis} sparse {tx_out_file_full} {tx_out_file}")
        message = "Converting %s to sparse format." % tx_out_file_full
        do.run(cmd.format(**locals()), message)
        if has_umi_matrix(data):
            cmd = ("{umis} sparse {tx_umi_matrix_full} {tx_umi_matrix}")
            message = "Converting %s to sparse format." % tx_umi_matrix_full
        do.run(cmd.format(**locals()), message)
    data = dd.set_count_file(data, out_file)
    return [[data]]

def get_barcode_metadata(data):
    barcode_file = dd.get_barcode_file(data)
    df = pd.read_csv(barcode_file, sep=",", header=0)
    barcodes = df["barcodes"]

def convert_to_kallisto(data):
    files = dd.get_input_sequence_files(data)
    if len(files) == 2:
        fq1, fq2 = files
    else:
        fq1, fq2 = files[0], None
    samplename = dd.get_sample_name(data)
    work_dir = dd.get_work_dir(data)
    kallisto_dir = os.path.join(work_dir, "kallisto", samplename, "fastq")
    out_file = os.path.join(kallisto_dir, "barcodes.batch")
    umis = _umis_cmd(data)
    if file_exists(out_file):
        return out_file
    if dd.get_minimum_barcode_depth(data):
        cb_histogram = os.path.join(work_dir, "umis", samplename, "cb-histogram.txt")
        cb_cutoff = dd.get_minimum_barcode_depth(data)
        cb_options = "--cb_histogram {cb_histogram} --cb_cutoff {cb_cutoff}"
        cb_options = cb_options.format(**locals())
    else:
        cb_options = ""
    cmd = ("{umis} kallisto {cb_options} --out_dir {tx_kallisto_dir} {fq1}")
    with file_transaction(data, kallisto_dir) as tx_kallisto_dir:
        safe_makedir(tx_kallisto_dir)
        message = ("Transforming %s to Kallisto singlecell format. "
                   % fq1)
        do.run(cmd.format(**locals()), message)
    return out_file

def demultiplex_samples(data):
    """
    demultiplex a fastqtransformed FASTQ file into separate sample barcode files
    """
    work_dir = os.path.join(dd.get_work_dir(data), "umis")
    sample_dir = os.path.join(work_dir, dd.get_sample_name(data))
    demulti_dir = os.path.join(sample_dir, "demultiplexed")

    files = data["files"]
    if len(files) == 2:
        logger.error("Sample demultiplexing doesn't handle paired-end reads, but "
            "we can add it. Open an issue here https://github.com/bcbio/bcbio-nextgen/issues if you need this and we'll add it.")
        sys.exit(1)
    else:
        fq1 = files[0]
    # check if samples need to be demultiplexed
    with open_fastq(fq1) as in_handle:
        read = next(in_handle)
        if "SAMPLE_" not in read:
            return [[data]]

    bcfile = get_sample_barcodes(dd.get_sample_barcodes(data), sample_dir)
    demultiplexed = glob.glob(os.path.join(demulti_dir, "*.fq*"))
    if demultiplexed:
        return [split_demultiplexed_sampledata(data, demultiplexed)]
    umis = _umis_cmd(data)
    cmd = ("{umis} demultiplex_samples --nedit 1 --barcodes {bcfile} "
           "--out_dir {tx_dir} {fq1}")
    msg = "Demultiplexing {fq1}."
    with file_transaction(data, demulti_dir) as tx_dir:
        do.run(cmd.format(**locals()), msg.format(**locals()))
    demultiplexed = glob.glob(os.path.join(demulti_dir, "*.fq*"))
    return [split_demultiplexed_sampledata(data, demultiplexed)]

def split_demultiplexed_sampledata(data, demultiplexed):
    """
    splits demultiplexed samples into separate entries in the global sample
    datadict
    """
    datadicts = []
    samplename = dd.get_sample_name(data)
    for fastq in demultiplexed:
        barcode = os.path.basename(fastq).split(".")[0]
        datadict = copy.deepcopy(data)
        datadict = dd.set_sample_name(datadict, samplename + "-" + barcode)
        datadict = dd.set_description(datadict, samplename + "-" + barcode)
        datadict["rgnames"]["rg"] = samplename + "-" + barcode
        datadict["name"]= ["", samplename + "-" + barcode]
        datadict["files"] = [fastq]
        datadicts.append(datadict)
    return datadicts

def concatenate_sparse_counts(*samples):
    samples = concatenate_sparse_matrices(samples, deduped=True)
    samples = concatenate_sparse_matrices(samples, deduped=False)
    samples = concatenate_cb_histograms(samples)
    return samples

def concatenate_sparse_matrices(samples, deduped=True):
    work_dir = dd.get_in_samples(samples, dd.get_work_dir)
    umi_dir = os.path.join(work_dir, "umis")
    if deduped:
        out_file = os.path.join(umi_dir, "tagcounts.mtx")
    else:
        out_file = os.path.join(umi_dir, "tagcounts-dupes.mtx")
    if file_exists(out_file):
        if deduped:
            newsamples = []
            for data in dd.sample_data_iterator(samples):
                newsamples.append([dd.set_combined_counts(data, out_file)])
            return newsamples
        else:
            return samples
    files = [dd.get_count_file(data) for data in
            dd.sample_data_iterator(samples)
            if dd.get_count_file(data)]
    if not deduped:
        files = [os.path.splitext(x)[0] + "-dupes.mtx" for x in files]

    files = [fn for fn in files if file_exists(fn)]
    descriptions = [dd.get_sample_name(data) for data in
                    dd.sample_data_iterator(samples) if dd.get_count_file(data)]
    if not files:
        return samples
    counts = SparseMatrix()
    counts.read(filename=files.pop(), colprefix=descriptions.pop())
    for filename, description in zip(files, descriptions):
        newcounts = SparseMatrix()
        newcounts.read(filename=filename, colprefix=description)
        counts.cat(newcounts)
    counts.write(out_file)
    newsamples = []
    if deduped:
        for data in dd.sample_data_iterator(samples):
            newsamples.append([dd.set_combined_counts(data, out_file)])
        return newsamples
    return samples

def concatenate_cb_histograms(samples):
    work_dir = dd.get_in_samples(samples, dd.get_work_dir)
    umi_dir = os.path.join(work_dir, "umis")
    out_file = os.path.join(umi_dir, "cb-histogram.txt")

    files = [dd.get_histogram_counts(data) for data in
            dd.sample_data_iterator(samples)
            if dd.get_histogram_counts(data)]
    files = " ".join(files)
    cmd = "cat {files} > {out_file}"
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            message = "Concat cellular barcode histograms: %s." % files
            do.run(cmd.format(**locals()), message)
    newsamples = []
    for data in dd.sample_data_iterator(samples):
        newsamples.append([dd.set_combined_histogram(data, out_file)])
    return newsamples

def version(data):
    umis = _umis_cmd(data)
    version_cmd = "%s version" % umis
    try:
        output = subprocess.check_output(version_cmd, shell=True).decode().strip()
    except:
        output = None
    return output

def has_umi_matrix(data):
    umis_version = version(data)
    if not version:
        return False
    return LooseVersion(umis_version) >= "1.0.0"

def filter_barcode_histogram(filtered_out_file, out_file, cutoff):
    if file_exists(filtered_out_file):
        return filtered_out_file
    sample_name = os.path.basename(os.path.dirname(out_file))
    with file_transaction(filtered_out_file) as tx_out_file:
        with open(tx_out_file, "w") as outh:
            with open(out_file) as inh:
                for line in inh:
                    barcode, reads = line.strip().split()
                    if int(reads) > cutoff:
                        outh.write("%s-%s\t%s\n" % (sample_name, barcode, reads))
    return filtered_out_file

def is_transformed(fastq):
    """
    check the first 100 reads to see if a FASTQ file has already been transformed
    by umis
    """

    with open_fastq(fastq) as in_handle:
        for line in islice(in_handle, 400):
            if "UMI_" in line:
                return True
    return False

def use_installed_transcriptome(data):
    user_fa = dd.get_transcriptome_fasta(data)
    user_gtf = dd.get_transcriptome_gtf(data)
    if not user_fa and not user_gtf:
        return True
    else:
        return False
