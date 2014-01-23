"""Handle disambiguation of reads from a chimeric input, splitting by organism.

Given specification of mixed input samples, splits a sample into multiple
sub-samples for alignment to individual genomes, then runs third-party disambiguation
scripts to reconcile.

Uses disambiguation scripts contributed by AstraZeneca, incorporated into bcbio-nextgen:
https://github.com/mjafin/disambiguate
"""
import collections
import copy
import os

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline.disambiguate.run import main as disambiguate_main
from bcbio.pipeline import run_info
from bcbio.provenance import do
from bcbio import bam


def split(items):
    """Split samples into all possible genomes for alignment.
    """
    out = []
    for data in [x[0] for x in items]:
        dis_orgs = data["config"]["algorithm"].get("disambiguate")
        if dis_orgs:
            data["disambiguate"] = {"genome_build": data["genome_build"],
                                    "base": True}
            out.append([data])
            for dis_org in dis_orgs:
                dis_data = copy.deepcopy(data)
                dis_data["disambiguate"] = {"genome_build": dis_org}
                dis_data["genome_build"] = dis_org
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
            to_process[data["name"][-1]].append(data)
        else:
            out.append([data])
    return out + run_parallel("run_disambiguate",
                              [(xs, xs[0]["config"]) for xs in to_process.itervalues()])

def run(items, config):
    """Run third party disambiguation script, resolving into single set of calls.
    """
    assert len(items) == 2, "Can only resolve two organism disambiguation"
    # check aligner, handling tophat/tophat2 distinctions
    aligner = config["algorithm"].get("aligner")
    aligner = "tophat" if aligner.startswith("tophat") else aligner
    assert aligner in ["bwa", "tophat"], "Disambiguation only supported for bwa and tophat alignments."
    if items[0]["disambiguate"].get("base"):
        data_a, data_b = items
    else:
        data_b, data_a = items
    work_bam_a = bam.sort(data_a["work_bam"], config, "queryname")
    work_bam_b = bam.sort(data_b["work_bam"], config, "queryname")
    out_dir = os.path.normpath(os.path.join(os.path.dirname(work_bam_a),
                                            os.pardir, os.pardir, "disambiguate"))
    base_name = os.path.join(out_dir, os.path.splitext(os.path.basename(work_bam_a))[0])
    summary_file = "%s_summary.txt" % base_name
    if not utils.file_exists(summary_file):
        with file_transaction(out_dir) as tx_out_dir:
            Args = collections.namedtuple("Args", "A B output_dir intermediate_dir "
                                          "no_sort prefix aligner")
            args = Args(work_bam_a, work_bam_b, tx_out_dir, tx_out_dir,
                        True, "", aligner)
            disambiguate_main(args)
    data_a["disambiguate"] = \
      {data_b["genome_build"]: "%s.disambiguatedSpeciesB.bam" % base_name,
       "%s-ambiguous" % data_a["genome_build"]: "%s.ambiguousSpeciesA.bam" % base_name,
       "%s-ambiguous" % data_b["genome_build"]: "%s.ambiguousSpeciesB.bam" % base_name,
       "summary": summary_file}
    data_a["work_bam"] = bam.sort("%s.disambiguatedSpeciesA.bam" % base_name, config)
    return [[data_a]]

def run_cplusplus(items, config):
    """Run third party disambiguation script, resolving into single set of calls.
    """
    assert len(items) == 2, "Can only resolve two organism disambiguation"
    # check aligner, handling tophat/tophat2 distinctions
    aligner = config["algorithm"].get("aligner")
    aligner = "tophat" if aligner.startswith("tophat") else aligner
    assert aligner in ["bwa", "tophat"], "Disambiguation only supported for bwa and tophat alignments."
    if items[0]["disambiguate"].get("base"):
        data_a, data_b = items
    else:
        data_b, data_a = items
    work_bam_a = bam.sort(data_a["work_bam"], config, "queryname")
    work_bam_b = bam.sort(data_b["work_bam"], config, "queryname")
    out_dir = os.path.normpath(os.path.join(os.path.dirname(work_bam_a),
                                            os.pardir, os.pardir, "disambiguate"))
    base_name = os.path.join(out_dir, os.path.splitext(os.path.basename(work_bam_a))[0])
    summary_file = "%s_summary.txt" % base_name
    if not utils.file_exists(summary_file):
        with file_transaction(out_dir) as tx_out_dir:
            raise NotImplementedError("Still need to test and support C++ version")
            cmd = ""
            do.run(cmd.format(**locals()), "Disambiguation", data_a)
    data_a["disambiguate"] = \
      {data_b["genome_build"]: "%s.disambiguatedSpeciesB.bam" % base_name,
       "%s-ambiguous" % data_a["genome_build"]: "%s.ambiguousSpeciesA.bam" % base_name,
       "%s-ambiguous" % data_b["genome_build"]: "%s.ambiguousSpeciesB.bam" % base_name,
       "summary": summary_file}
    data_a["work_bam"] = bam.sort("%s.disambiguatedSpeciesA.bam" % base_name, config)
    return [[data_a]]
