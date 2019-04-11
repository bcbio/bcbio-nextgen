#!/usr/bin/env python
"""Run LOH heterogeneity comparison amongst mutiple methods, focusing on HLA.

Includes LOHHLA with inputs from a bcbio OptiType hg38, PureCN, TitanCNA Cromwell run.

Requires:

    - pdftoppm (from poppler-utils) for plot generation from pdfs
"""
from __future__ import print_function
import collections
import csv
import glob
import StringIO as io
import os
import shutil
import subprocess
import sys

import yaml

from bcbio.pipeline import alignment
from bcbio import utils

HLA_GLOB="call-call_hla/shard-*/execution"
SV_GLOB="call-svcall/shard-*/wf-svcall.cwl/*/call-detect_sv/execution"
SVCALL_GLOB="structural/{sample}/{method}/*{ext}"
LOHHLA="../lohhla/lohhla"
ALIGNER="novoalign"
DEPTH_FILTER=5

# hg38 coordinates for HLA region https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC
hla_coords = ("chr6", 28510120, 33480577)

def run_sample(tumor, normal, work_dir, cromwell_dir, hla_fa):
    hla_fasta, hlas = prep_hla_ref(hla_fa, work_dir)
    hla_cromwell_dir = _get_cromwell_execution_dir(cromwell_dir, HLA_GLOB)
    sv_cromwell_dir = _get_cromwell_execution_dir(cromwell_dir, SV_GLOB)
    tumor_fastq, tumor_calls_orig = get_hla(tumor, hla_cromwell_dir, HLA_GLOB)
    tumor_bam = alignment.align_to_sort_bam(tumor_fastq, None, ALIGNER, get_data(tumor, hla_fasta, work_dir))["work_bam"]
    normal_fastq, normal_calls_orig = get_hla(normal, hla_cromwell_dir, HLA_GLOB)
    normal_bam = alignment.align_to_sort_bam(normal_fastq, None, ALIGNER, get_data(normal, hla_fasta, work_dir))["work_bam"]
    tumor_ploidy = prep_ploidy(work_dir, tumor, tumor_bam, sv_cromwell_dir, os.path.join(SV_GLOB, SVCALL_GLOB))
    tumor_calls = prep_hla(work_dir, tumor, normal_calls_orig, hlas, normal_bam, tumor_bam)
    normal_calls = prep_hla(work_dir, normal, normal_calls_orig, hlas, normal_bam, tumor_bam)

    bam_dir, normal_bam_ready = create_tumor_bamdir(tumor, tumor_bam, normal_bam, work_dir)
    out_dir = utils.safe_makedir(os.path.join(work_dir, tumor, "lohhla_out"))
    prep_bam_inputs(out_dir, tumor, tumor_calls, tumor_bam)
    prep_bam_inputs(out_dir, normal, normal_calls, normal_bam)
    lohhla_output = os.path.join(out_dir, "%s.%s.DNA.HLAlossPrediction_CI.xls" % (tumor, DEPTH_FILTER))
    cmd = ["Rscript", os.path.join(LOHHLA, "LOHHLAscript.R"),
           "--patientId", tumor, "--outputDir", out_dir,
           "--normalBAMfile", normal_bam_ready, "--BAMDir", bam_dir,
           "--hlaPath", normal_calls, "--HLAfastaLoc", hla_fasta,
           "--HLAexonLoc", os.path.join(LOHHLA, "data", "hla.dat"),
           "--CopyNumLoc", tumor_ploidy,
           "--mappingStep", "FALSE",
           "--minCoverageFilter", str(DEPTH_FILTER)]
    if not os.path.exists(lohhla_output):
        subprocess.check_call(cmd)
    compare_calls(tumor, lohhla_output, sv_cromwell_dir, os.path.join(SV_GLOB, SVCALL_GLOB))

def _create_plot(tumor, in_glob, out_ext, page=1):
    """Create an output plot for the given PDF in the images directory.
    """
    out_dir = utils.safe_makedir("images")
    out_name = os.path.join(out_dir, "%s-%s" % (tumor, out_ext))
    in_file = glob.glob(in_glob)[0]
    cmd = ["pdftoppm", in_file, out_name, "-png", "-f", page, "-singlefile"]
    if not os.path.exists(out_name + ".png"):
        subprocess.check_call([str(x) for x in cmd])
    return out_name + ".png"

def _get_loh_from_calls(calls):
    if calls["loh"]:
        return "mixed LOH" if calls["std"] else "LOH"
    else:
        return "no LOH"

def _compare_lohhla(lohhla_output):
    print("#### LOHHLA")
    print("```")
    seen = set([])
    calls = collections.defaultdict(int)
    with open(lohhla_output) as in_handle:
        header = in_handle.readline().strip().split("\t")
        for c in in_handle:
            vals = dict(zip(header, c.strip().split("\t")))
            key = (vals["HLA_A_type1"], vals["HLA_A_type2"])
            if key not in seen:
                print([vals[x] for x in ["PVal_unique", "HLA_A_type1", "HLA_A_type2", "HLA_type1copyNum_withBAFBin", "HLA_type2copyNum_withBAFBin"]])
                seen.add(key)
                if float(vals["PVal_unique"]) < 0.01:
                    calls["loh"] += 1
                else:
                    calls["std"] += 1
    print("```")
    return _get_loh_from_calls(calls)

def _compare_purecn(tumor, cromwell_dir, sv_glob):
    print("#### PureCN")
    calls = collections.defaultdict(int)
    pure_base_file = _get_cromwell_file(cromwell_dir, sv_glob, dict(sample=tumor, method="purecn", ext="-purecn.csv"))
    pure_cn_file = _get_cromwell_file(cromwell_dir, sv_glob, dict(sample=tumor, method="purecn", ext="loh.csv"))
    cov_plot = _create_plot(tumor, os.path.join(os.path.dirname(pure_cn_file), "%s*-purecn.pdf" % tumor), "purecn", 2)
    sun_plot = _create_plot(tumor, os.path.join(os.path.dirname(pure_cn_file), "%s*-purecn_local_optima.pdf" % tumor),
                            "purecn-sunrise")
    with open(pure_base_file) as in_handle:
        vals = dict(zip(in_handle.readline().strip().replace('"', '').split(","),
                        in_handle.readline().strip().split(",")))
        print()
        print("|  |  |")
        print("| ---  | --- |")
        print("| purity | %s |" % vals["Purity"])
        print("| ploidy | %s |" % vals["Ploidy"])
    print("```")
    with open(pure_cn_file) as in_handle:
        in_handle.readline()  # header
        for line in in_handle:
            _, chrom, start, end, _, cn, minor_cn = line.split(",")[:7]
            start = int(start)
            end = int(end)
            if chrom == hla_coords[0] and are_overlapping((start, end), hla_coords[1:]):
                print(line.strip().split(",")[1:])
                if int(minor_cn) == 0:
                    calls["loh"] += 1
                else:
                    calls["std"] += 1
    print("```")
    print("![%s PureCN coverage](%s)" % (tumor, cov_plot))
    print("![%s PureCN sunrise](%s)" % (tumor, sun_plot))
    return _get_loh_from_calls(calls)

def _compare_titancna(tumor, cromwell_dir, sv_glob):
    print("#### TitanCNA")
    calls = collections.defaultdict(int)
    titancna_file = _get_cromwell_file(cromwell_dir, sv_glob, dict(sample=tumor, method="titancna", ext="Clusters.txt"))
    with open(titancna_file) as in_handle:
        vals = dict(zip(in_handle.readline().strip().split("\t"), in_handle.readline().strip().split("\t")))
        path = vals["path"]
        init_dir, check_dir = os.path.dirname(titancna_file).split("/", 1)
        rel_path = init_dir + "/" + path[path.find(check_dir):]
        print()
        print("|  |  |")
        print("| ---  | --- |")
        print("| purity | %s |" % vals["purity"])
        print("| ploidy | %s |" % vals["ploidy"])
    cna_plot = _create_plot(tumor, os.path.join(rel_path, "%s*_CNA.pdf" % tumor), "titan-cna")
    loh_plot = _create_plot(tumor, os.path.join(rel_path, "%s*_LOH.pdf" % tumor), "titan-loh")
    seg_file = rel_path + ".segs.txt"
    out_keys = ["Chromosome", "Start_Position.bp.", "End_Position.bp.", "Copy_Number",
                "MinorCN", "MajorCN", "TITAN_call"]
    print("```")
    with open(seg_file) as in_handle:
        header = in_handle.readline().strip().split()
        for line in in_handle:
            val = dict(zip(header, line.strip().split()))
            start = int(val["Start_Position.bp."])
            end = int(val["End_Position.bp."])
            if val["Chromosome"] == hla_coords[0] and are_overlapping((start, end), hla_coords[1:]):
                print([val[k] for k in out_keys])
                if int(val["MinorCN"]) == 0:
                    calls["loh"] += 1
                else:
                    calls["std"] += 1
    print("```")
    print("![%s TitanCNA CNA](%s)" % (tumor, cna_plot))
    print("![%s TitanCNA LOH](%s)" % (tumor, loh_plot))
    return _get_loh_from_calls(calls)

def _compare_gatkcnv(tumor, cromwell_dir, sv_glob):
    print("#### GATK CNV")
    gatk_file = _get_cromwell_file(cromwell_dir, sv_glob, dict(sample=tumor, method="gatk-cnv", ext="-call.seg"))
    orig_model_file = gatk_file.replace("-call.seg", ".modeled.png")
    model_file = os.path.join("images", os.path.basename(orig_model_file))
    shutil.copy(orig_model_file, model_file)
    print("```")
    with open(gatk_file) as in_handle:
        for line in in_handle:
            if not line.startswith("@"):
                chrom, start, end = line.split()[:3]
                if chrom == hla_coords[0] and are_overlapping((int(start), int(end)), hla_coords[1:]):
                    print(line.strip())
    print("```")
    print("![%s GATK CNV](%s)" % (tumor, model_file))

def compare_calls(tumor, lohhla_output, cromwell_dir, sv_glob):
    summary = collections.OrderedDict()
    print("### %s" % tumor)
    orig_stdout = sys.stdout
    sys.stdout = io.StringIO()
    summary["LOHHLA"] = _compare_lohhla(lohhla_output)
    summary["PureCN"] = _compare_purecn(tumor, cromwell_dir, sv_glob)
    summary["TitanCNA"] = _compare_titancna(tumor, cromwell_dir, sv_glob)
    saved_stdout = sys.stdout
    sys.stdout = orig_stdout
    print()
    print("|  |  |")
    print("| ---  | --- |")
    for k, v in summary.items():
        print("| %s | %s |" % (k, v))
    sys.stdout.write(saved_stdout.getvalue())

    # print("#### CNVkit")
    # cnvkit_file = _get_cromwell_file(cromwell_dir, sv_glob, dict(sample=tumor, method="cnvkit", ext="-call.cns"))
    # out_keys = ["chromosome", "start", "end", "cn", "cn1", "cn2"]
    # print("```")
    # with open(cnvkit_file) as in_handle:
    #     header = in_handle.readline().strip().split()
    #     for line in in_handle:
    #         chrom, start, end = line.split()[:3]
    #         if chrom == hla_coords[0] and are_overlapping((int(start), int(end)), hla_coords[1:]):
    #             vals = dict(zip(header, line.strip().split()))
    #             print([vals[k] for k in out_keys])
    # print("```")

    print

def are_overlapping(r, s):
    """Test if two coordinates overlap.
    https://stackoverflow.com/a/27182551
    """
    return r[1] >= s[0] and s[1] >= r[0]

def _get_cromwell_file(cromwell_dir, file_glob, kwargs):
    fglob = os.path.join(cromwell_dir, file_glob.format(**kwargs))
    fs = glob.glob(fglob)
    assert len(fs) == 1, (fglob, fs)
    return fs[0]

def _get_cromwell_execution_dir(base_dir, target_glob):
    """Retrieve the baseline directory with cromwell output files.

    Handles Cromwell restarts where there are multiple work directories and
    we traverse symlinks back to the original.
    """
    cur_dir = glob.glob(os.path.join(base_dir, target_glob))[0]
    if os.path.exists(os.path.join(cur_dir, "cwl.output.json")):
        return base_dir
    else:
        symlink_dir = os.path.dirname(os.path.realpath(os.path.join(cur_dir, "script")))
        ref_base = os.path.dirname(base_dir)
        new_guid = symlink_dir[symlink_dir.find(ref_base) + len(ref_base) + 1:].split("/")[0]
        return _get_cromwell_execution_dir(os.path.join(ref_base, new_guid), target_glob)

def prep_bam_inputs(out_dir, sample, call_file, bam_file):
    """Prepare expected input BAM files from pre-aligned.
    """
    base = utils.splitext_plus(os.path.basename(bam_file))[0]
    with open(call_file) as in_handle:
        for cur_hla in (x.strip() for x in in_handle):
            out_file = os.path.join(utils.safe_makedir(os.path.join(out_dir, base)),
                                    "%s.type.%s.filtered.bam" % (base, cur_hla))
            if not os.path.exists(out_file):
                cmd = ["samtools", "view", "-b","-o", out_file, bam_file, cur_hla]
                subprocess.check_call(cmd)

def create_tumor_bamdir(tumor, tumor_bam, normal_bam, work_dir):
    """Create expected input directory with tumor/normal BAMs in one place.
    """
    bam_dir = utils.safe_makedir(os.path.join(work_dir, tumor, "in_bams"))
    normal_bam_ready = os.path.join(bam_dir, os.path.basename(normal_bam))
    utils.symlink_plus(normal_bam, normal_bam_ready)
    tumor_bam_ready = os.path.join(bam_dir, os.path.basename(tumor_bam))
    utils.symlink_plus(tumor_bam, tumor_bam_ready)
    return bam_dir, normal_bam_ready

def get_data(sample, hla_fasta, work_dir):
    return {"dirs": {"work": work_dir},
            "config": {"analysis": "variant", "algorithm": {"multiple_mappers": "All 9999"},
                       "resources": {"novoalign": {"options": ["-R", "0"]}}},
            "reference": {"bwa": {"indexes": hla_fasta + ".bwt"},
                          "novoalign": {"indexes": [hla_fasta + ".ndx"]}},
            "rgnames": {"sample": sample, "rg": sample, "pl": "illumina", "pu": sample, "lane": sample}}


def get_hla(sample, cromwell_dir, hla_glob):
    """Retrieve HLA calls and input fastqs for a sample.
    """
    hla_dir = glob.glob(os.path.join(cromwell_dir, hla_glob, "align", sample, "hla"))[0]
    fastq = os.path.join(hla_dir, "OptiType-HLA-A_B_C-input.fq")
    calls = os.path.join(hla_dir, "%s-optitype.csv" % sample)
    return fastq, calls

def name_to_absolute(x):
    """Convert standard hg38 HLA name into ABSOLUTE naming.
    """
    for c in ["-", "*", ":"]:
        x = x.replace(c, "_")
    x = x.lower()
    return x

def get_hla_choice(h, hlas, normal_bam, tumor_bam):
    """Retrieve matching HLA with best read support in both tumor and normal
    """
    def get_counts(bam_file):
        counts = {}
        for line in subprocess.check_output(["samtools", "idxstats", bam_file]).split("\n"):
            if line.startswith(h):
                name, _, count, _ = line.split()
                counts[name] = int(count)
        return counts
    tcounts = get_counts(tumor_bam)
    ncounts = get_counts(normal_bam)
    check_hlas = [x for x in hlas if x.startswith(h) and tcounts.get(x, 0) > 0 and ncounts.get(x, 0) > 0]
    cur_hlas = sorted(check_hlas, key=lambda x: ncounts[x], reverse=True)
    #print(cur_hlas[0], tcounts.get(cur_hlas[0]), ncounts.get(cur_hlas[0]))
    return cur_hlas[0]

def prep_hla(work_dir, sample, calls, hlas, normal_bam, tumor_bam):
    """Convert HLAs into ABSOLUTE format for use with LOHHLA.

    LOHHLA hard codes names to hla_a, hla_b, hla_c so need to move
    """
    work_dir = utils.safe_makedir(os.path.join(work_dir, sample, "inputs"))
    hla_file = os.path.join(work_dir, "%s-hlas.txt" % sample)
    with open(calls) as in_handle:
        with open(hla_file, "w") as out_handle:
            next(in_handle)  # header
            for line in in_handle:
                _, _, a, _, _ = line.strip().split(",")
                a1, a2 = a.split(";")
                out_handle.write(get_hla_choice(name_to_absolute(a1), hlas, normal_bam, tumor_bam) + "\n")
                out_handle.write(get_hla_choice(name_to_absolute(a2), hlas, normal_bam, tumor_bam) + "\n")
    return hla_file

def prep_ploidy(work_dir, sample, bam_file, cromwell_dir, sv_glob):
    """Create LOHHLA compatible input ploidy file from PureCN output.
    """
    purecn_file = _get_cromwell_file(cromwell_dir, sv_glob, dict(sample=sample, method="purecn", ext="purecn.csv"))
    work_dir = utils.safe_makedir(os.path.join(work_dir, sample, "inputs"))
    out_file = os.path.join(work_dir, "%s-solutions.txt" % sample)
    with open(purecn_file) as in_handle:
        reader = csv.reader(in_handle)
        purecn_stats = dict(zip(next(reader), next(reader)))
    with open(out_file, "w") as out_handle:
        out_handle.write("Ploidy\ttumorPurity\ttumorPloidy\n")
        lohhla_name = utils.splitext_plus(os.path.basename(bam_file))[0]
        out_handle.write("%s\t%s\t%s\t%s\n" % (lohhla_name, purecn_stats["Ploidy"],
                                               purecn_stats["Purity"], purecn_stats["Ploidy"]))
    return out_file

def prep_hla_ref(hla_fasta, work_dir):
    work_dir = utils.safe_makedir(os.path.join(work_dir, "hlaref"))
    out_file = os.path.join(work_dir, os.path.basename(hla_fasta))
    seen_names = set([])
    if not utils.file_uptodate(out_file, hla_fasta):
        with open(hla_fasta) as in_handle:
            with open(out_file, "w") as out_handle:
                for line in in_handle:
                    if line.startswith(">"):
                        cur_name = name_to_absolute(line.strip().split()[1])
                        if cur_name not in seen_names:
                            out_handle.write(">%s\n" % cur_name)
                            seen_names.add(cur_name)
                            write_seq = True
                        else:
                            write_seq = False
                    elif write_seq:
                        out_handle.write(line)
    if not os.path.exists(out_file + ".bwt"):
        subprocess.check_call(["bwa", "index", out_file])
    if not os.path.exists(out_file + ".ndx"):
        subprocess.check_call(["novoindex", out_file + ".ndx", out_file])
    hlas = []
    with open(out_file) as in_handle:
        for line in in_handle:
            if line.startswith(">"):
                hlas.append(line[1:].strip())
    return out_file, hlas

def samples_from_config(sample_yaml):
    with open(sample_yaml) as in_handle:
        config = yaml.safe_load(in_handle)
    by_batch = collections.defaultdict(dict)
    for s in config["details"]:
        by_batch[s["metadata"]["batch"]][s["metadata"]["phenotype"]] = s["description"]
    for bid in sorted(by_batch.keys()):
        yield by_batch[bid]["tumor"], by_batch[bid]["normal"]

if __name__ == "__main__":
    sample_config, hla_fa, cromwell_dir = sys.argv[1:]
    work_dir = utils.safe_makedir(os.path.join(os.getcwd(), "work_lohhla"))
    for t, n in sorted(samples_from_config(sample_config)):
        run_sample(t, n, work_dir, cromwell_dir, hla_fa)
