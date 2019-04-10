"""Handle disambiguation of reads from a chimeric input, splitting by organism.

Given specification of mixed input samples, splits a sample into multiple
sub-samples for alignment to individual genomes, then runs third-party
disambiguation scripts to reconcile.

Uses disambiguation scripts contributed by AstraZeneca, incorporated into
bcbio-nextgen: https://github.com/mjafin/disambiguate
"""

from __future__ import print_function
import collections
import copy
import os

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline.disambiguate.run import main as disambiguate_main
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils, merge, run_info
from bcbio.provenance import do
from bcbio import bam


import six


def split(*items):
    """Split samples into all possible genomes for alignment.
    """
    out = []
    for data in [x[0] for x in items]:
        dis_orgs = data["config"]["algorithm"].get("disambiguate")
        if dis_orgs:
            if not data.get("disambiguate", None):
                data["disambiguate"] = {"genome_build": data["genome_build"],
                                        "base": True}
            out.append([data])
            # handle the instance where a single organism is disambiguated
            if isinstance(dis_orgs, six.string_types):
                dis_orgs = [dis_orgs]
            for dis_org in dis_orgs:
                dis_data = copy.deepcopy(data)
                dis_data["disambiguate"] = {"genome_build": dis_org}
                dis_data["genome_build"] = dis_org
                dis_data["config"]["algorithm"]["effects"] = False
                dis_data = run_info.add_reference_resources(dis_data)
                out.append([dis_data])
        else:
            out.append([data])
    return out

def resolve(items, run_parallel):
    """Combine aligned and split samples into final set of disambiguated reads.
    """
    out = []
    to_process = collections.defaultdict(list)
    for data in [x[0] for x in items]:
        if "disambiguate" in data:
            split_part = tuple([int(x) for x in data["align_split"].split("-")]) if data.get("combine") else None
            to_process[(dd.get_sample_name(data), split_part)].append(data)
        else:
            out.append([data])
    if len(to_process) > 0:
        dis1 = run_parallel("run_disambiguate",
                            [(xs, xs[0]["config"]) for xs in to_process.values()])
        disambigs_by_name = collections.defaultdict(list)
        print(len(dis1))
        for xs in dis1:
            assert len(xs) == 1
            data = xs[0]
            disambigs_by_name[dd.get_sample_name(data)].append(data)
        dis2 = run_parallel("disambiguate_merge_extras",
                            [(xs, xs[0]["config"]) for xs in disambigs_by_name.values()])
    else:
        dis2 = []
    return out + dis2

def merge_extras(items, config):
    """Merge extra disambiguated reads into a final BAM file.
    """
    final = {}
    for extra_name in items[0]["disambiguate"].keys():
        in_files = []
        for data in items:
            in_files.append(data["disambiguate"][extra_name])
        out_file = "%s-allmerged%s" % os.path.splitext(in_files[0])
        if in_files[0].endswith(".bam"):
            merged_file = merge.merge_bam_files(in_files, os.path.dirname(out_file), items[0],
                                                out_file=out_file)
        else:
            assert extra_name == "summary", extra_name
            merged_file = _merge_summary(in_files, out_file, items[0])
        final[extra_name] = merged_file
    out = []
    for data in items:
        data["disambiguate"] = final
        out.append([data])
    return out

def _merge_summary(in_files, out_file, data):
    """Create one big summary file for disambiguation from multiple splits.
    """
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                for i, in_file in enumerate(in_files):
                    with open(in_file) as in_handle:
                        for j, line in enumerate(in_handle):
                            if j == 0:
                                if i == 0:
                                    out_handle.write(line)
                            else:
                                out_handle.write(line)
    return out_file

def run(items, config):
    """Run third party disambiguation script, resolving into single set of calls.
    """
    assert len(items) == 2, "Can only resolve two organism disambiguation"
    # check aligner, handling tophat/tophat2 distinctions
    aligner = config["algorithm"].get("aligner")
    aligner = "tophat" if aligner.startswith("tophat") else aligner
    assert aligner in ["bwa", "hisat2", "tophat", "star"], "Disambiguation only supported for bwa, hisat2, star and tophat alignments."
    if items[0]["disambiguate"].get("base"):
        data_a, data_b = items
    else:
        data_b, data_a = items
    work_bam_a = bam.sort(data_a["work_bam"], config, "queryname")
    work_bam_b = bam.sort(data_b["work_bam"], config, "queryname")
    if data_a.get("align_split"):
        base_dir = utils.safe_makedir(os.path.normpath(os.path.join(os.path.dirname(work_bam_a),
                                                                    os.pardir, os.pardir,
                                                                    "disambiguate_%s" % aligner)))
        out_dir = os.path.join(base_dir, "_".join([str(x) for x in data_a["align_split"].split("-")]))
    else:
        out_dir = os.path.normpath(os.path.join(os.path.dirname(work_bam_a),
                                                os.pardir, "disambiguate_%s" % aligner))
    base_name = os.path.join(out_dir, os.path.splitext(os.path.basename(work_bam_a))[0])
    summary_file = "%s_summary.txt" % base_name
    if not utils.file_exists(summary_file):
        with file_transaction(items[0], out_dir) as tx_out_dir:
            _run_cplusplus(work_bam_a, work_bam_b, tx_out_dir, aligner, os.path.basename(base_name), items)
    data_a["disambiguate"] = \
      {data_b["genome_build"]: bam.sort("%s.disambiguatedSpeciesB.bam" % base_name, config),
       "%s-ambiguous" % data_a["genome_build"]: bam.sort("%s.ambiguousSpeciesA.bam" % base_name, config),
       "%s-ambiguous" % data_b["genome_build"]: bam.sort("%s.ambiguousSpeciesB.bam" % base_name, config),
       "summary": summary_file}
    data_a["work_bam"] = bam.sort("%s.disambiguatedSpeciesA.bam" % base_name, config)
    bam.index(dd.get_work_bam(data_a), data_a["config"])
    return [[data_a]]

def _run_python(work_bam_a, work_bam_b, out_dir, aligner, prefix, items):
    """Run python version of disambiguation
    """
    Args = collections.namedtuple("Args", "A B output_dir intermediate_dir "
                                    "no_sort prefix aligner")
    args = Args(work_bam_a, work_bam_b, out_dir, out_dir, True, "", aligner)
    disambiguate_main(args)

def _run_cplusplus(work_bam_a, work_bam_b, out_dir, aligner, prefix, items):
    """Run third party disambiguation script, resolving into single set of calls.
    """
    ngs_disambiguate = config_utils.get_program("ngs_disambiguate", items[0]["config"])
    cmd = [ngs_disambiguate, "--no-sort", "--prefix", prefix, "--aligner", aligner,
           "--output-dir", out_dir, work_bam_a, work_bam_b]
    do.run(cmd, "Disambiguation", items[0])
