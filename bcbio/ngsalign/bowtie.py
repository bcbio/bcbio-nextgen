"""Next gen sequence alignments with Bowtie (http://bowtie-bio.sourceforge.net).
"""
import os

from bcbio import utils
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.ngsalign import alignprep, novoalign, postalign
from bcbio.provenance import do

galaxy_location_file = "bowtie_indices.loc"

def _bowtie_args_from_config(data):
    """Configurable high level options for bowtie.
    """
    config = data['config']
    qual_format = config["algorithm"].get("quality_format", "")
    if qual_format.lower() == "illumina":
        qual_flags = ["--phred64-quals"]
    else:
        qual_flags = []
    multi_mappers = config["algorithm"].get("multiple_mappers", True)
    multi_flags = ["-M", 1] if multi_mappers else ["-m", 1]
    multi_flags = [] if data["analysis"].lower().startswith("smallrna-seq") else multi_flags
    cores = config.get("resources", {}).get("bowtie", {}).get("cores", None)
    num_cores = config["algorithm"].get("num_cores", 1)
    core_flags = ["-p", str(num_cores)] if num_cores > 1 else []
    return core_flags + qual_flags + multi_flags

def align(fastq_file, pair_file, ref_file, names, align_dir, data,
          extra_args=None):
    """Do standard or paired end alignment with bowtie.
    """
    num_hits = 1
    if data["analysis"].lower().startswith("smallrna-seq"):
        num_hits = 1000
    config = data['config']
    out_file = os.path.join(align_dir, "{0}-sort.bam".format(dd.get_sample_name(data)))
    if data.get("align_split"):
        final_file = out_file
        out_file, data = alignprep.setup_combine(final_file, data)
        fastq_file, pair_file = alignprep.split_namedpipe_cls(fastq_file, pair_file, data)
    else:
        final_file = None
        if fastq_file.endswith(".gz"):
            fastq_file = "<(gunzip -c %s)" % fastq_file
            if pair_file:
                pair_file = "<(gunzip -c %s)" % pair_file

    if not utils.file_exists(out_file) and (final_file is None or not utils.file_exists(final_file)):
        with postalign.tobam_cl(data, out_file, pair_file is not None) as (tobam_cl, tx_out_file):
            cl = [config_utils.get_program("bowtie", config)]
            cl += _bowtie_args_from_config(data)
            cl += extra_args if extra_args is not None else []
            cl += ["-q",
                   "-v", 2,
                   "-k", num_hits,
                   "-X", 2000, # default is too selective for most data
                   "--best",
                   "--strata",
                   "--sam",
                   ref_file]
            if pair_file:
                cl += ["-1", fastq_file, "-2", pair_file]
            else:
                cl += [fastq_file]
            cl = [str(i) for i in cl]
            fix_rg_cmd = r"samtools addreplacerg -r '%s' -" % novoalign.get_rg_info(names)
            cmd = " ".join(cl) + " | " + fix_rg_cmd + " | " + tobam_cl
            do.run(cmd, "Running Bowtie on %s and %s." % (fastq_file, pair_file), data)
    return out_file
