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
import shutil

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline.disambiguate.run import main as disambiguate_main
from bcbio.pipeline.disambiguate.pdxfilter import PDXFilter
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils, merge, run_info
from bcbio.provenance import do
from bcbio.pipeline import merge, run_info
from bcbio import bam

from bcbio.log import logger


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
            if isinstance(dis_orgs, basestring):
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

        disamb_config = [xs[0]["config"] for xs in to_process.itervalues()][0]
        logger.debug("Disambiguation resources {}".format(str(disamb_config)))

        dis1 = run_parallel("run_disambiguate",
                            [(xs, xs[0]["config"])
                             for xs in to_process.itervalues()])
        disambigs = []
        for xs in dis1:
            assert len(xs) == 1
            disambigs.append(xs[0])

        final = collections.defaultdict(dict)

        # Iterate over disambiguation types (explant, human, ambigious, summary)
        disambiguation_config = disambigs[0]["config"]
        disambiguation_types  = disambigs[0]["disambiguate"].keys()

        # Iterate over disambiguated sample splits and extract sample names
        items_by_sample = collections.defaultdict(list)
        for data in disambigs:
            items_by_sample[dd.get_sample_name(data)].append(data)

        to_merge  = []
        to_delete = set()

        # Split by disambiguation type
        for disambiguation_type in disambiguation_types:

            for sample_name, type_items in items_by_sample.items():

                in_files = []
                for data in type_items:
                    in_files.append(data["disambiguate"][disambiguation_type])

                merged_dir    = os.path.dirname(os.path.dirname(in_files[0]))
                merged_base   = '_'.join(os.path.basename(in_files[0]).split('_')[:2])
                merged_suffix = os.path.splitext(in_files[0])[1]
                merged_type   = disambiguation_type.replace(' ', '_')
                merged_out    = os.path.join(merged_dir, merged_base + '.' + merged_type +
                                             '.allmerged' + merged_suffix)

                delete_dir    = os.path.dirname(in_files[0])
                to_delete.add(delete_dir)

                if in_files[0].endswith(".bam"):
                    to_merge.append([in_files, merged_out, disambigs[0]])

                elif disambiguation_type == "summary":
                    _merge_summary(in_files, merged_out, type_items[0])

                else:
                    continue

                final[sample_name][disambiguation_type] = merged_out

        # Merge BAM files
        run_parallel("disambiguate_merge_extras", to_merge)

        # TODO: delete inputs to merge (the whole directory)
        # Warning: requires filtering by existing merged_out at the beginning of this function
        #for delete_dir in to_delete:
        #    os.rmtree(delete_dir)

        for data in disambigs:
            data["disambiguate"] = final[dd.get_sample_name(data)]
            out.append([data])

    return out


def merge_extras(in_files, out_file, config):
    """Merge extra disambiguated reads into a final BAM file.
    """

    merged_file = merge.merge_bam_files(in_files,
                                        os.path.dirname(out_file),
                                        config,
                                        out_file=out_file)
    return merged_file


def _merge_summary(in_files, out_file, data):
    """Create one big summary file for disambiguation from multiple splits.
    Drop the headers on all files except the first one.
    """
    if utils.file_exists(out_file):
        return out_file

    with file_transaction(data, out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            for in_file in in_files:
                with open(in_file) as in_handle:
                    if in_file != in_files[0]:
                        next(in_handle)  # skip the header row
                    for line in in_handle:
                        out_handle.write(line)
    return out_file


def run(items, config):
    """Run third party disambiguation script, resolving into single set of calls.
    """
    assert len(items) == 2, "Can only resolve two organism disambiguation"
    # check aligner, handling tophat/tophat2 distinctions
    aligner = config["algorithm"].get("aligner")
    if items[0]["disambiguate"].get("base"):
        data_a, data_b = items
    else:
        data_b, data_a = items

    # Construct name of sorted input files
    work_bam_a_nsorted = os.path.splitext(data_a["work_bam"])[0] + '.nsorted.bam'
    work_bam_b_nsorted = os.path.splitext(data_b["work_bam"])[0] + '.nsorted.bam'

    # logger.info('Disambiguate prep of input BAM {} and {}'.format(work_bam_a_nsorted, work_bam_b_nsorted))
    if data_a.get("align_split"):
        base_dir = utils.safe_makedir(os.path.normpath(
            os.path.join(os.path.dirname(work_bam_a_nsorted), os.pardir, os.pardir,
                         "disambiguate_%s" % aligner)))
        logger.info('Disambiguate prep of prepped work bam BAM {} with base dir {}'.format(work_bam_a_nsorted, base_dir))
        split_name = "_".join([str(x) for x in data_a["align_split"].split("-")])
        out_dir = os.path.join(base_dir, split_name)
        logger.info('Disambiguate prep of prepped work bam BAM {} with out dir {}'.format(work_bam_a_nsorted, out_dir))
    else:
        out_dir = os.path.normpath(os.path.join(os.path.dirname(work_bam_a_nsorted),
                                                os.pardir,
                                                "disambiguate_%s" % aligner))

    base_name = os.path.join(out_dir,
                             os.path.splitext(os.path.basename(work_bam_a_nsorted))[0])
    logger.info('Disambiguate prep of prepped work bam BAM {} with base name {}'.format(work_bam_a_nsorted, base_name))

    summary_file = "%s_summary.txt" % base_name
    explant_bam = "%s.explant.sorted.bam" % base_name
    ambiguous_bam = "%s.ambiguous.sorted.bam" % base_name
    work_bam = "%s.human.sorted.bam" % base_name

    logger.info('Disambiguate prep with work bam {}'.format(work_bam))

    logger.info('Deciding if disambiguation is required. Checking for existence of {}, {}, {} and {}'.format(summary_file, explant_bam, ambiguous_bam, work_bam))

    if not utils.file_exists(summary_file) or not utils.file_exists(explant_bam) or not utils.file_exists(ambiguous_bam) or not utils.file_exists(work_bam):
        logger.info('Disambiguating work bam a {} since outputs are not already existing'.format(work_bam_a_nsorted))
        work_bam_a = bam.sort(data_a["work_bam"], config, "queryname")
        work_bam_b = bam.sort(data_b["work_bam"], config, "queryname")
        logger.info('Disambiguate run with work bam a {}'.format(work_bam_a))
        logger.info('Disambiguate run with work bam b {}'.format(work_bam_b))
        with file_transaction(items[0], out_dir) as tx_out_dir:
            logger.info('Disambiguate run with sorted prep work bam a {} and tx out dir {}'.format(work_bam_a_nsorted, tx_out_dir))
            tmp_base_name = os.path.join(tx_out_dir, os.path.basename(base_name))
            logger.info('Disambiguate run with sorted prep work bam a {} and tmp_base_name {}'.format(work_bam_a_nsorted, tmp_base_name))
            pdx_filter = PDXFilter(work_bam_a, work_bam_b,
                                   "%s.human.bam" % tmp_base_name,
                                   # Must be bam else it will not be merged
                                   "%s.explant.bam" % tmp_base_name,
                                   # Must be bam else it will not be merged
                                   "%s.ambiguous.bam" % tmp_base_name,
                                   # Must be bam else it will not be merged
                                   "%s_summary.txt" % tmp_base_name,
                                   hard_filter=True,
                                   debug=True)
            pdx_filter.run()

        # Perhaps this can be removed since it has been fixed in bcbio
        if data_a.get("align_split"):
            split_dir = os.path.join(out_dir, split_name)
            logger.info('Disambiguate post-run with sorted prep work bam a {} and split dir {}'.format(work_bam_a_nsorted, split_dir))
            if os.path.isdir(split_dir):
                for tmp_file in os.listdir(split_dir):
                    logger.info('Disambiguate post-run with sorted prep work bam a {} aiming to move file {}'.format(work_bam_a_nsorted, tmp_file))
                    src = os.path.join(split_dir, tmp_file)
                    if os.path.isfile(src):
                        dest = os.path.join(out_dir, tmp_file)
                        logger.info('Disambiguate post-run with sorted prep work bam a {} moving file {} from {} to {}'.format(work_bam_a_nsorted, tmp_file, src, dest))
                        shutil.move(src, dest)
                shutil.rmtree(split_dir)

        try:
            if work_bam_a != data_a["work_bam"]:
                os.remove(work_bam_a)
        except:
            pass
        try:
            if work_bam_b != data_b["work_bam"]:
                os.remove(work_bam_b)
        except:
            pass

    else:
        logger.info('Skipping disambiguation for work bam a {} since outputs are already existing'.format(work_bam_a_nsorted))

    explant_bam = os.path.isfile(explant_bam) and explant_bam or bam.sort("%s.explant.bam" % base_name, config)
    ambiguous_bam = os.path.isfile(ambiguous_bam) and ambiguous_bam or bam.sort("%s.ambiguous.bam" % base_name, config)
    work_bam = os.path.isfile(work_bam) and work_bam or bam.sort("%s.human.bam" % base_name, config)
    # logger.info('Disambiguate run with post work_bam {}'.format(work_bam))

    data_a["disambiguate"] = {data_b["genome_build"]: explant_bam,
                              "%s-ambiguous" % data_a["genome_build"]: ambiguous_bam,
                              "summary": summary_file}
    data_a["work_bam"] = work_bam
    try:
        os.remove("%s.explant.bam" % base_name)
    except:
        pass
    try:
        os.remove("%s.human.bam" % base_name)
    except:
        pass
    try:
        os.remove("%s.ambiguous.bam" % base_name)
    except:
        pass

    return [[data_a]]
