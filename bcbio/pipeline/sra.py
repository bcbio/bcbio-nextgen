"""Deal with sra Id names as input"""
import os
import glob
import shutil
import subprocess
import json
import re

from bcbio.log import logger
from bcbio import utils
from bcbio.provenance import do
from bcbio.bam.fastq import is_fastq, combine_pairs
from bcbio.pipeline import fastq

def is_gsm(fn):
    p = re.compile("^GSM[0-9]+$")
    if p.match(fn) and not utils.file_exists(fn):
        return True

def _query_info(db, ids):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db={0}\&id={1}\&retmode=json".format(db, ids)
    cmd = "curl {0}".format(url)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out = process.stdout.read()
    data = json.loads(out)
    return data

def query_gsm(gsm, out_file, config = {}):
    gsm = gsm[0]
    out_dir = os.path.dirname(os.path.abspath(out_file))
    name = utils.splitext_plus(os.path.basename(out_file))[0]
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds\&term={0}\&retmode=json".format(gsm)
    cmd = "curl {0}".format(url)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out = process.stdout.read()
    data = json.loads(out)
    ids = data.get("esearchresult", {}).get("idlist", [])
    logger.debug("Get id sample for %s" % gsm)
    if ids:
        gsm_info = _query_info("gds", ids[-1])
        srxlist = gsm_info.get("result", {}).get(ids[-1], {}).get("extrelations", {})
        srxall = []
        for srxe in srxlist:
            if srxe.get("targetftplink", None):
                srxall.append(srxe["targetftplink"])
        logger.debug("Get FTP link for %s : %s" % (ids[-1], srxall))
        outs = []
        for srx in srxall:
            sra_dir = utils.safe_makedir(os.path.join(out_dir, name))
            srafiles = _download_srx(gsm, srx, sra_dir)
            logger.debug("Get SRA for %s: %s" % (gsm, " ".join(srafiles)))
            if srafiles:
                for sra in srafiles:
                    outs.extend(_convert_fastq(sra, out_dir))
            logger.debug("Get FASTQ for %s: %s" % (gsm, " ".join(outs)))
        if outs:
            files = combine_pairs(outs)
            out_file = fastq.merge(files, out_file, config)
            return out_file

def _create_link(sraid):
    sraprex = sraid[0:5]
    url = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/{sraprex}/{sraid}"
    return url.format(**locals())

def _download_srx(srxid, url, out_dir):
    cmd = "wget -N -r -nH -nd -np -nv {0}".format(url)
    out_dir = os.path.abspath(utils.safe_makedir(out_dir))
    with utils.chdir(out_dir):
        do.run(cmd, "Download %s" % url )
        # return [os.path.abspath(fn) for fn in glob.glob("*sra")]
    return [os.path.join(out_dir, fn) for fn in os.listdir(out_dir)]

def _download_sra(sraid, outdir):
    url = _create_link(sraid)
    cmd = "wget -O {out_file} {url}"
    out_file = os.path.join(outdir, "%s.sra" % sraid)
    if not utils.file_exists(out_file):
        do.run(cmd, "Download %s" % sraid)
    return out_file

def _convert_fastq(srafn, outdir, single=False):
    "convert sra to fastq"
    cmd = "fastq-dump --split-files --gzip {srafn}"
    sraid = os.path.basename(utils.splitext_plus(srafn)[0])
    if not single:
        out_file = [os.path.join(outdir, "%s_1.fastq.gz" % sraid),
                    os.path.join(outdir, "%s_2.fastq.gz" % sraid)]
        if not utils.file_exists(out_file[0]):
            with utils.chdir(outdir):
                do.run(cmd.format(**locals()), "Covert to fastq %s" % sraid)
        if not utils.file_exists(out_file[0]):
            raise IOError("SRA %s didn't convert, something happened." % srafn)
        return [out for out in out_file if utils.file_exists(out)]
    else:
        raise ValueError("Not supported single-end sra samples for now.")
