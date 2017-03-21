"""Quality control using Qualimap.

http://qualimap.bioinfo.cipf.es/
"""
import glob
import os
import shutil

import pandas as pd
import pybedtools
import toolz as tz
import toolz.dicttoolz as dtz

from bcbio.log import logger
from bcbio import bam, utils
from bcbio.ngsalign import postalign
from bcbio.provenance import do
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.distributed.transaction import file_transaction
from bcbio.rnaseq import gtf
from bcbio.variation import bedutils

# ## Standard Qualimap

def run(bam_file, data, out_dir):
    """Run qualimap to assess alignment quality metrics.
    """
    # Qualimap results should be saved to a directory named after sample.
    # MultiQC (for parsing additional data) picks the sample name after the dir as follows:
    #   <sample name>/raw_data_qualimapReport/insert_size_histogram.txt
    results_dir = os.path.join(out_dir, dd.get_sample_name(data))
    resources = config_utils.get_resources("qualimap", data["config"])
    options = " ".join(resources.get("options", ""))
    results_file = os.path.join(results_dir, "genome_results.txt")
    report_file = os.path.join(results_dir, "qualimapReport.html")
    utils.safe_makedir(results_dir)
    pdf_file = "qualimapReport.pdf"
    if not utils.file_exists(results_file) and not utils.file_exists(os.path.join(results_dir, pdf_file)):
        if "qualimap_full" in tz.get_in(("config", "algorithm", "tools_on"), data, []):
            logger.info("Full qualimap analysis for %s may be slow." % bam_file)
            ds_bam = bam_file
        else:
            ds_bam = bam.downsample(bam_file, data, 1e7, work_dir=out_dir)
            bam_file = ds_bam if ds_bam else bam_file
        if options.find("PDF") > -1:
            options = "%s -outfile %s" % (options, pdf_file)
        num_cores = data["config"]["algorithm"].get("num_cores", 1)
        qualimap = config_utils.get_program("qualimap", data["config"])
        max_mem = config_utils.adjust_memory(resources.get("memory", "1G"),
                                             num_cores)

        with file_transaction(data, results_dir) as tx_results_dir:
            utils.safe_makedir(tx_results_dir)

            export = utils.local_path_export()
            cmd = ("unset DISPLAY && {export} {qualimap} bamqc -bam {bam_file} -outdir {tx_results_dir} "
                   "--skip-duplicated --skip-dup-mode 0 "
                   "-nt {num_cores} --java-mem-size={max_mem} {options}")
            species = None
            if (tz.get_in(("genome_resources", "aliases", "human"), data, "")
                  or dd.get_genome_build(data).startswith(("hg", "GRCh"))):
                species = "HUMAN"
            elif dd.get_genome_build(data).startswith(("mm", "GRCm")):
                species = "MOUSE"
            if species in ["HUMAN", "MOUSE"]:
                cmd += " -gd {species}"
            regions = (dd.get_coverage(data) if dd.get_coverage(data) not in [None, False, "None"]
                       else dd.get_variant_regions_merged(data))
            if regions:
                regions = bedutils.merge_overlaps(bedutils.clean_file(regions, data), data)
                bed6_regions = _bed_to_bed6(regions, out_dir)
                cmd += " -gff {bed6_regions}"
            bcbio_env = utils.get_bcbio_env()
            do.run(cmd.format(**locals()), "Qualimap: %s" % dd.get_sample_name(data), env=bcbio_env)
            tx_results_file = os.path.join(tx_results_dir, "genome_results.txt")
            cmd = "sed -i 's/bam file = .*/bam file = %s.bam/' %s" % (dd.get_sample_name(data), tx_results_file)
            do.run(cmd, "Fix Name Qualimap for {}".format(dd.get_sample_name(data)))
    # Qualimap output folder (results_dir) needs to be named after the sample (see comments above). However, in order 
    # to keep its name after upload, we need to put  the base QC file (results_file) into the root directory (out_dir):
    base_results_file = os.path.join(out_dir, os.path.basename(results_file))
    shutil.copyfile(results_file, base_results_file)
    return {"base": base_results_file,
            "secondary": _find_qualimap_secondary_files(results_dir, base_results_file)}

def _parse_qualimap_metrics(report_file, data):
    """Extract useful metrics from the qualimap HTML report file.
    """
    if not utils.file_exists(report_file):
        return {}
    import lxml.html
    out = {}
    parsers = {"Globals": _parse_qualimap_globals,
               "Globals (inside of regions)": _parse_qualimap_globals_inregion,
               "Coverage": _parse_qualimap_coverage,
               "Coverage (inside of regions)": _parse_qualimap_coverage,
               "Insert size": _parse_qualimap_insertsize,
               "Insert size (inside of regions)": _parse_qualimap_insertsize}
    root = lxml.html.parse(report_file).getroot()
    for table in root.xpath("//div[@class='table-summary']"):
        header = table.xpath("h3")[0].text
        if header in parsers:
            out.update(parsers[header](table))
    new_names = []
    for metric in out:
        if "qualimap_full" not in tz.get_in(("config", "algorithm", "tools_on"), data, []):
            metric += "_qualimap_1e7reads_est"
        new_names.append(metric)
    out = dict(zip(new_names, out.values()))
    return out

def _parse_num_pct(k, v):
    num, pct = v.split(" / ")
    return {k: num.replace(",", "").strip(), "%s pct" % k: pct.strip()}

def _parse_qualimap_globals(table):
    """Retrieve metrics of interest from globals table.
    """
    out = {}
    want = {"Mapped reads": _parse_num_pct,
            "Duplication rate": lambda k, v: {k: v}}
    for row in table.xpath("table/tr"):
        col, val = [x.text for x in row.xpath("td")]
        if col in want:
            out.update(want[col](col, val))
    return out

def _parse_qualimap_globals_inregion(table):
    """Retrieve metrics from the global targeted region table.
    """
    out = {}
    for row in table.xpath("table/tr"):
        col, val = [x.text for x in row.xpath("td")]
        if col == "Mapped reads":
            out.update(_parse_num_pct("%s (in regions)" % col, val))
    return out

def _parse_qualimap_coverage(table):
    """Parse summary qualimap coverage metrics.
    """
    out = {}
    for row in table.xpath("table/tr"):
        col, val = [x.text for x in row.xpath("td")]
        if col == "Mean":
            out["Coverage (Mean)"] = val
    return out

def _parse_qualimap_insertsize(table):
    """Parse insert size metrics.
    """
    out = {}
    for row in table.xpath("table/tr"):
        col, val = [x.text for x in row.xpath("td")]
        if col == "Median":
            out["Insert size (Median)"] = val
    return out

def _bed_to_bed6(orig_file, out_dir):
    """Convert bed to required bed6 inputs.
    """
    bed6_file = os.path.join(out_dir, "%s-bed6%s" % os.path.splitext(os.path.basename(orig_file)))
    if not utils.file_exists(bed6_file):
        with open(bed6_file, "w") as out_handle:
            for i, region in enumerate(list(x) for x in pybedtools.BedTool(orig_file)):
                region = [x for x in list(region) if x]
                fillers = [str(i), "1.0", "+"]
                full = region + fillers[:6 - len(region)]
                out_handle.write("\t".join(full) + "\n")
    return bed6_file

# ## RNAseq Qualimap

def _parse_metrics(metrics):
    # skipped metrics can sometimes be in unicode, replace unicode with NA if it exists
    metrics = dtz.valmap(lambda x: 'nan' if isinstance(x, unicode) else x, metrics)

    # missing = set(["Genes Detected", "Transcripts Detected", "Mean Per Base Cov."])
    correct = set(["rRNA", "rRNA_rate"])
    percentages = set(["Intergenic pct", "Intronic pct", "Exonic pct"])
    to_change = dict({"5'-3' bias": 1,
                      "Intergenic pct": "Intergenic Rate",
                      "Intronic pct": "Intronic Rate",
                      "Exonic pct": "Exonic Rate",
                      "Duplication Rate of Mapped": 1,
                      "Average_insert_size": 1,
                      })
    total = ["Not aligned", "Aligned to genes", "No feature assigned"]

    out = {}
    total_reads = sum([int(metrics[name]) for name in total])
    out.update({key: val for key, val in metrics.items() if key in correct})
    [metrics.update({name: 1.0 * float(metrics[name]) / 100}) for name in
     percentages]

    for name in to_change:
        if not to_change[name]:
            continue
        try:
            if to_change[name] == 1:
                out.update({name: float(metrics[name])})
            else:
                out.update({to_change[name]: float(metrics[name])})
        # if we can't convert metrics[name] to float (?'s or other non-floats)
        except ValueError:
            continue
    return out

def _detect_duplicates(bam_file, out_dir, data):
    """
    count duplicate percentage
    """
    out_file = os.path.join(out_dir, "dup_metrics.txt")
    if not utils.file_exists(out_file):
        dup_align_bam = postalign.dedup_bam(bam_file, data)
        num_cores = dd.get_num_cores(data)
        with file_transaction(data, out_file) as tx_out_file:
            sambamba = config_utils.get_program("sambamba", data, default="sambamba")
            dup_count = ("{sambamba} view --nthreads {num_cores} --count "
                         "-F 'duplicate and not unmapped' "
                         "{dup_align_bam} >> {tx_out_file}")
            message = "Counting duplicates in {bam_file}.".format(bam_file=bam_file)
            do.run(dup_count.format(**locals()), message)
            tot_count = ("{sambamba} view --nthreads {num_cores} --count "
                         "-F 'not unmapped' "
                         "{dup_align_bam} >> {tx_out_file}")
            message = "Counting reads in {bam_file}.".format(bam_file=bam_file)
            do.run(tot_count.format(**locals()), message)
    with open(out_file) as in_handle:
        dupes = float(in_handle.next().strip())
        total = float(in_handle.next().strip())
    return {"Duplication Rate of Mapped": dupes / total}

def _transform_browser_coor(rRNA_interval, rRNA_coor):
    """
    transform interval format to browser coord: chr:start-end
    """
    with open(rRNA_coor, 'w') as out_handle:
        with open(rRNA_interval, 'r') as in_handle:
            for line in in_handle:
                c, bio, source, s, e = line.split("\t")[:5]
                if bio.startswith("rRNA"):
                    out_handle.write(("{0}:{1}-{2}\n").format(c, s, e))

def _detect_rRNA(data):
    sample = dd.get_sample_name(data)
    gtf_file = dd.get_gtf_file(data)
    tidy_file = dd.get_sailfish_tidy(data)
    rrna_features = gtf.get_rRNA(gtf_file)
    transcripts = set([x[1] for x in rrna_features if x])
    if not (transcripts and tidy_file):
        return {'rRNA': "NA", "rRNA_rate": "NA"}
    count_table = pd.read_csv(tidy_file, sep="\t", dtype={'sample': str})
    sample_table = count_table[count_table["sample"].isin([sample])]
    rrna_exp = map(float, sample_table[sample_table["id"].isin(transcripts)]["numreads"])
    total_exp = map(float, sample_table["numreads"])
    rrna = sum(rrna_exp)
    if sum(total_exp) == 0:
        return {'rRNA': str(rrna), 'rRNA_rate': "NA"}
    rrna_rate = float(rrna) / sum(total_exp)
    return {'rRNA': str(rrna), 'rRNA_rate': str(rrna_rate)}

def _parse_qualimap_rnaseq(table):
    """
    Retrieve metrics of interest from globals table.
    """
    out = {}
    for row in table.xpath("table/tr"):
        col, val = [x.text for x in row.xpath("td")]
        col = col.replace(":", "").strip()
        val = val.replace(",", "")
        m = {col: val}
        if val.find("/") > -1:
            m = _parse_num_pct(col, val.replace("%", ""))
        out.update(m)
    return out

def _parse_rnaseq_qualimap_metrics(report_file):
    """Extract useful metrics from the qualimap HTML report file.
    """
    import lxml.html
    out = {}
    parsers = ["Reads alignment", "Reads genomic origin", "Transcript coverage profile"]
    root = lxml.html.parse(report_file).getroot()
    for table in root.xpath("//div[@class='table-summary']"):
        header = table.xpath("h3")[0].text
        if header in parsers:
            out.update(_parse_qualimap_rnaseq(table))
    return out

def run_rnaseq(bam_file, data, out_dir):
    """
    Run qualimap for a rnaseq bam file and parse results
    """
    strandedness = {"firststrand": "strand-specific-reverse",
                    "secondstrand": "strand-specific-forward",
                    "unstranded": "non-strand-specific"}

    # Qualimap results should be saved to a directory named after sample.
    # MultiQC (for parsing additional data) picks the sample name after the dir as follows:
    #   <sample name>/raw_data_qualimapReport/insert_size_histogram.txt
    results_dir = os.path.join(out_dir, dd.get_sample_name(data))
    results_file = os.path.join(results_dir, "rnaseq_qc_results.txt")
    report_file = os.path.join(results_dir, "qualimapReport.html")
    config = data["config"]
    gtf_file = dd.get_gtf_file(data)
    single_end = not bam.is_paired(bam_file)
    library = strandedness[dd.get_strandedness(data)]
    if not utils.file_exists(results_file):
        with file_transaction(data, results_dir) as tx_results_dir:
            utils.safe_makedir(tx_results_dir)
            bam.index(bam_file, config)
            cmd = _rnaseq_qualimap_cmd(data, bam_file, tx_results_dir, gtf_file, single_end, library)
            do.run(cmd, "Qualimap for {}".format(dd.get_sample_name(data)))
            tx_results_file = os.path.join(tx_results_dir, "rnaseq_qc_results.txt")
            cmd = "sed -i 's/bam file = .*/bam file = %s.bam/' %s" % (dd.get_sample_name(data), tx_results_file)
            do.run(cmd, "Fix Name Qualimap for {}".format(dd.get_sample_name(data)))
    metrics = _parse_rnaseq_qualimap_metrics(report_file)
    metrics.update(_detect_duplicates(bam_file, results_dir, data))
    metrics.update(_detect_rRNA(data))
    metrics.update({"Average_insert_size": bam.estimate_fragment_size(bam_file)})
    metrics = _parse_metrics(metrics)
    # Qualimap output folder (results_dir) needs to be named after the sample (see comments above). However, in order
    # to keep its name after upload, we need to put  the base QC file (results_file) into the root directory (out_dir):
    base_results_file = os.path.join(out_dir, os.path.basename(results_file))
    shutil.copyfile(results_file, base_results_file)
    return {"base": base_results_file,
            "secondary": _find_qualimap_secondary_files(results_dir, base_results_file),
            "metrics": metrics}

def _rnaseq_qualimap_cmd(data, bam_file, out_dir, gtf_file=None, single_end=None, library="non-strand-specific"):
    """
    Create command lines for qualimap
    """
    config = data["config"]
    qualimap = config_utils.get_program("qualimap", config)
    resources = config_utils.get_resources("qualimap", config)
    num_cores = resources.get("cores", dd.get_num_cores(data))
    max_mem = config_utils.adjust_memory(resources.get("memory", "2G"),
                                         num_cores)
    export = utils.local_path_export()
    cmd = ("unset DISPLAY && {export} {qualimap} rnaseq -outdir {out_dir} "
           "-a proportional -bam {bam_file} -p {library} "
           "-gtf {gtf_file} --java-mem-size={max_mem}").format(**locals())
    return cmd

def _find_qualimap_secondary_files(results_dir, base_file):
    """Retrieve additional files, avoiding double uploading the base file.
    """
    def not_dup(x):
        is_dup = (os.path.basename(x) == os.path.basename(base_file) and
                  os.path.getsize(x) == os.path.getsize(base_file))
        return not is_dup
    return filter(not_dup,
                  glob.glob(os.path.join(results_dir, 'qualimapReport.html')) +
                  glob.glob(os.path.join(results_dir, '*.txt')) +
                  glob.glob(os.path.join(results_dir, "css", "*")) +
                  glob.glob(os.path.join(results_dir, "raw_data_qualimapReport", "*")) +
                  glob.glob(os.path.join(results_dir, "images_qualimapReport", "*")))
