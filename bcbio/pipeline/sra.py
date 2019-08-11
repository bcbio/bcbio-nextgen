"""Deal with sra Id names as input"""
import os
import subprocess
import json
import re
import traceback

from bcbio.log import logger
from bcbio import utils
from bcbio.provenance import do
from bcbio.bam.fastq import combine_pairs
from bcbio.pipeline import fastq

def is_gsm(fn):
    p = re.compile(r"^GSM[0-9]+$")
    if p.match(fn) and not utils.file_exists(fn):
        return True

def is_srr(fn):
    p = re.compile(r"^SRR[0-9]+$")
    if p.match(fn) and not utils.file_exists(fn):
        return True

def _query_info(db, ids):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={0}\&id={1}".format(db, ids)
    cmd = "curl {0}".format(url)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out = process.stdout.read()
    data= []
    for line in out.split("RUN_SET")[1].split("RUN"):
        if line.find("accession=") > -1:
            srr = line.split("accession=")[1].split(" ")[0].replace("\"", "")
            data.append(srr)
    return data

def query_gsm(gsm, out_file, config = {}):
    gsm = gsm[0]
    out_dir = os.path.dirname(os.path.abspath(out_file))
    name = utils.splitext_plus(os.path.basename(out_file))[0]
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra\&term={0}\&retmode=json".format(gsm)
    cmd = "curl {0}".format(url)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out = process.stdout.read()
    data = json.loads(out)
    ids = data.get("esearchresult", {}).get("idlist", [])
    logger.debug("Get id sample for %s" % gsm)
    if ids:
        gsm_info = _query_info("sra", ids[-1])
        logger.debug("gsm_info:%s" % gsm_info)
        srrall = []
        for srr in gsm_info:
            srrall.append(_create_link(srr))
        logger.debug("Get FTP link for %s : %s" % (ids[-1], srrall))
        outs = []
        for srx in srrall:
            sra_dir = utils.safe_makedir(os.path.join(out_dir, name))
            srafiles = _download_srx(srx, sra_dir)
            if srafiles:
                logger.debug("Get SRA for %s: %s" % (gsm, " ".join(srafiles)))
                for sra in srafiles:
                    fastq_fn = _convert_fastq(sra, out_dir)
                    if fastq_fn:
                        outs.extend(fastq_fn)
            logger.debug("Get FASTQ for %s: %s" % (gsm, " ".join(outs)))
        if outs:
            files = combine_pairs(outs)
            out_file = fastq.merge(files, out_file, config)
            return out_file

def query_srr(sra, out_file, config = {}):
    sra = sra[0]
    outs = []
    out_dir = os.path.dirname(os.path.abspath(out_file))
    name = utils.splitext_plus(os.path.basename(out_file))[0]
    srrall = []
    for srr in sra:
        srrall.append(_create_link(srr))
    logger.debug("Get FTP link for %s : %s" % (name, srrall))
    for srx in srrall:
        sra_dir = utils.safe_makedir(os.path.join(out_dir, name))
        srafiles = _download_srx(srx, sra_dir)
        if srafiles:
            logger.debug("Get SRA for %s: %s" % (sra, " ".join(srafiles)))
            for sra in srafiles:
                fastq_fn = _convert_fastq(sra, out_dir)
                if fastq_fn:
                    outs.extend(fastq_fn)
        logger.debug("Get FASTQ for %s: %s" % (sra, " ".join(outs)))
    if outs:
        files = combine_pairs(outs)
        out_file = fastq.merge(files, out_file, config)
        return out_file

def _create_link(sraid):
    sraprex = sraid[0:6]
    url = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/{sraprex}/{sraid}"
    return url.format(**locals())

def _download_srx(url, out_dir):
    cmd = "wget -N -r -nH -nd -np -nv {0}".format(url)
    out_dir = os.path.abspath(utils.safe_makedir(out_dir))
    with utils.chdir(out_dir):
        try:
            do.run(cmd, "Download %s" % url )
        except:
            logger.warning("Sample path not found in database. Skipping.")
            traceback.print_exc()
            return None
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
    cmd = "%s %s" % (utils.local_path_export(), cmd)
    sraid = os.path.basename(utils.splitext_plus(srafn)[0])
    if not srafn:
        return None
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
