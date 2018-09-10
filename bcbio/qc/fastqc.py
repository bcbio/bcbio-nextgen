"""Run and parse output from FastQC.

http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
"""
import os
import shutil

import pandas as pd
try:
    from fadapa import Fadapa
except ImportError:
    Fadapa = None

from bcbio import bam, utils
from bcbio.distributed.transaction import tx_tmpdir
from bcbio.log import logger
from bcbio.provenance import do
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils

def run(bam_file, data, fastqc_out):
    """Run fastqc, generating report in specified directory and parsing metrics.

    Downsamples to 10 million reads to avoid excessive processing times with large
    files, unless we're running a Standard/smallRNA-seq/QC pipeline.

    Handles fastqc 0.11+, which use a single HTML file and older versions that use
    a directory of files + images. The goal is to eventually move to only 0.11+
    """
    sentry_file = os.path.join(fastqc_out, "fastqc_report.html")
    if not os.path.exists(sentry_file):
        work_dir = os.path.dirname(fastqc_out)
        utils.safe_makedir(work_dir)
        ds_file = (bam.downsample(bam_file, data, 1e7, work_dir=work_dir)
                   if data.get("analysis", "").lower() not in ["standard", "smallrna-seq"]
                   else None)
        if ds_file is not None:
            bam_file = ds_file
        frmt = "bam" if bam_file.endswith("bam") else "fastq"
        fastqc_name = utils.splitext_plus(os.path.basename(bam_file))[0]
        fastqc_clean_name = dd.get_sample_name(data)
        num_cores = data["config"]["algorithm"].get("num_cores", 1)
        with tx_tmpdir(data, work_dir) as tx_tmp_dir:
            with utils.chdir(tx_tmp_dir):
                cl = [config_utils.get_program("fastqc", data["config"]),
                      "-d", tx_tmp_dir,
                      "-t", str(num_cores), "--extract", "-o", tx_tmp_dir, "-f", frmt, bam_file]
                cl = "%s %s %s" % (utils.java_freetype_fix(),
                                   utils.local_path_export(), " ".join([str(x) for x in cl]))
                do.run(cl, "FastQC: %s" % dd.get_sample_name(data))
                tx_fastqc_out = os.path.join(tx_tmp_dir, "%s_fastqc" % fastqc_name)
                tx_combo_file = os.path.join(tx_tmp_dir, "%s_fastqc.html" % fastqc_name)
                if not os.path.exists(sentry_file) and os.path.exists(tx_combo_file):
                    utils.safe_makedir(fastqc_out)
                    # Use sample name for reports instead of bam file name
                    with open(os.path.join(tx_fastqc_out, "fastqc_data.txt"), 'r') as fastqc_bam_name, \
                            open(os.path.join(tx_fastqc_out, "_fastqc_data.txt"), 'w') as fastqc_sample_name:
                        for line in fastqc_bam_name:
                            fastqc_sample_name.write(line.replace(os.path.basename(bam_file), fastqc_clean_name))
                    shutil.move(os.path.join(tx_fastqc_out, "_fastqc_data.txt"), os.path.join(fastqc_out, 'fastqc_data.txt'))
                    shutil.move(tx_combo_file, sentry_file)
                    if os.path.exists("%s.zip" % tx_fastqc_out):
                        shutil.move("%s.zip" % tx_fastqc_out, os.path.join(fastqc_out, "%s.zip" % fastqc_clean_name))
                elif not os.path.exists(sentry_file):
                    raise ValueError("FastQC failed to produce output HTML file: %s" % os.listdir(tx_tmp_dir))
    logger.info("Produced HTML report %s" % sentry_file)
    parser = FastQCParser(fastqc_out, dd.get_sample_name(data))
    stats = parser.get_fastqc_summary()
    parser.save_sections_into_file()
    return stats

class FastQCParser:
    def __init__(self, base_dir, sample=None):
        self._dir = base_dir
        self.sample = sample

    def get_fastqc_summary(self):
        ignore = set(["Total Sequences", "Filtered Sequences",
                      "Filename", "File type", "Encoding"])
        stats = {}
        for stat_line in self._fastqc_data_section("Basic Statistics")[1:]:
            k, v = stat_line.split("\t")[:2]
            if k not in ignore:
                stats[k] = v
        return stats

    def _fastqc_data_section(self, section_name):
        out = []
        in_section = False
        data_file = os.path.join(self._dir, "fastqc_data.txt")
        if os.path.exists(data_file):
            with open(data_file) as in_handle:
                for line in in_handle:
                    if line.startswith(">>%s" % section_name):
                        in_section = True
                    elif in_section:
                        if line.startswith(">>END"):
                            break
                        out.append(line.rstrip("\r\n"))
        return out

    def save_sections_into_file(self):

        data_file = os.path.join(self._dir, "fastqc_data.txt")
        if os.path.exists(data_file) and Fadapa:
            parser = Fadapa(data_file)
            module = [m[1] for m in parser.summary()][2:9]
            for m in module:
                out_file = os.path.join(self._dir, m.replace(" ", "_") + ".tsv")
                dt = self._get_module(parser, m)
                dt.to_csv(out_file, sep="\t", index=False)

    def _get_module(self, parser, module):
        """
        Get module using fadapa package
        """
        dt = []
        lines = parser.clean_data(module)
        header = lines[0]
        for data in lines[1:]:
            if data[0].startswith("#"):  # some modules have two headers
                header = data
                continue
            if data[0].find("-") > -1:  # expand positions 1-3 to 1, 2, 3
                f, s = map(int, data[0].split("-"))
                for pos in range(f, s):
                    dt.append([str(pos)] + data[1:])
            else:
                dt.append(data)
        dt = pd.DataFrame(dt)
        dt.columns = [h.replace(" ", "_") for h in header]
        dt['sample'] = self.sample
        return dt
