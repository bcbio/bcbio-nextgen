"""Runs Preseq, a tool that estimates the complexity of a library.
http://smithlabresearch.org/software/preseq/

Executes `lc_extrap` or `c_curve` commands that predict the yield
for future experiments (for higher or lower input correspondingly).

"""
import os
import math
import pybedtools
import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.provenance import do
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.bam import ref, readstats
from bcbio.variation import coverage as cov
from bcbio.qc import samtools


def run(bam_file, data, out_dir):
    out = {}
    preseq_cmd = tz.get_in(["config", "algorithm", "preseq"], data)
    if not preseq_cmd:
        return out

    samtools_stats_dir = os.path.join(out_dir, os.path.pardir, "samtools")
    samtools_stats = samtools.run(bam_file, data, samtools_stats_dir)["metrics"]

    stats_file = os.path.join(out_dir, "%s.txt" % dd.get_sample_name(data))
    if not utils.file_exists(stats_file):
        utils.safe_makedir(out_dir)
        preseq = config_utils.get_program("preseq", data["config"])
        params = _get_preseq_params(data, preseq_cmd, int(samtools_stats["Total_reads"]))
        param_line = "{options} -step {step} -seg_len {seg_len} "
        if preseq_cmd == "lc_extrap":
            param_line += "-extrap {extrap} "
        param_line = param_line.format(**params)
        with file_transaction(data, stats_file) as tx_out_file:
            cmd = "{preseq} {preseq_cmd} -bam -pe {bam_file} -o {tx_out_file} {param_line}".format(**locals())
            do.run(cmd.format(**locals()), "preseq " + preseq_cmd, data)

    out = _prep_real_counts(bam_file, data, samtools_stats)

    return {"base": stats_file,
            "metrics": out}

def _get_preseq_params(data, preseq_cmd, read_count):
    """ Get parameters through resources.
        If "step" or "extrap" limit are not provided, then calculate optimal values based on read count.
    """
    defaults = {
        'seg_len': 100000,        # maximum segment length when merging paired end bam reads
        'steps': 300,             # number of points on the plot
        'extrap_fraction': 3,     # extrapolate up to X times read_count
        'extrap': None,           # extrapolate up to X reads
        'step': None,             # step size (number of reads between points on the plot)
        'options': '',
    }
    params = {}

    main_opts = [("-e", "-extrap"), ("-l", "-seg_len"), ("-s", "-step")]
    other_opts = config_utils.get_resources("preseq", data["config"]).get("options", [])
    if isinstance(other_opts, str):
        other_opts = [other_opts]
    for sht, lng in main_opts:
        if sht in other_opts:
            i = other_opts.index(sht)
        elif lng in other_opts:
            i = other_opts.index(lng)
        else:
            i = None
        if i is not None:
            params[lng[1:]] = other_opts[i + 1]
            other_opts = other_opts[:i] + other_opts[i + 2:]
    params['options'] = ' '.join(other_opts)
    for k, v in config_utils.get_resources("preseq", data["config"]).items():
        if k != 'options':
            params[k] = v

    params['steps'] = params.get('steps', defaults['steps'])

    if preseq_cmd == 'c_curve':
        params['extrap_fraction'] = 1

    else:
        if params.get('step') is None:
            if params.get('extrap') is None:
                unrounded__extrap = read_count * params.get('extrap_fraction', defaults['extrap_fraction'])
                unrounded__step = unrounded__extrap // params['steps']
                if params.get('extrap_fraction') is not None:  # extrap_fraction explicitly provided
                    params['extrap'] = unrounded__extrap
                    params['step'] = unrounded__step
                else:
                    power_of_10 = 10 ** math.floor(math.log(unrounded__step, 10))
                    rounded__step = int(math.floor(unrounded__step // power_of_10) * power_of_10)
                    rounded__extrap = int(rounded__step) * params['steps']
                    params['step'] = rounded__step
                    params['extrap'] = rounded__extrap
            else:
                params['step'] = params['extrap'] // params['steps']

        elif params.get('extrap') is None:
            params['extrap'] = params['step'] * params['steps']

    params['step'] = params.get('step', defaults['step'])
    params['extrap'] = params.get('extrap', defaults['extrap'])
    params['seg_len'] = params.get('seg_len', defaults['seg_len'])

    logger.info("Preseq: running {steps} steps of size {step}, extap limit {extrap}".format(**params))
    return params

def _prep_real_counts(bam_file, data, samtools_stats):
    out = {}

    if dd.get_coverage(data) and dd.get_coverage(data) not in ["None"]:
        bed = dd.get_coverage_merged(data)
        target_name = "coverage"
    elif dd.get_coverage_interval(data) != "genome":
        bed = dd.get_variant_regions_merged(data) or dd.get_sample_callable(data)
        target_name = "variant_regions"
    else:
        bed = None
        target_name = "genome"

    dedupped = utils.get_in(data, ("config", "algorithm", "mark_duplicates"), True)

    if bed:
        out["Preseq_genome_size"] = pybedtools.BedTool(bed).total_coverage()
        out["Preseq_read_count"] = readstats.number_of_mapped_reads(
            data, bam_file, keep_dups=True, bed_file=bed, target_name=target_name)
        ontrg_unique_depth = cov.get_average_coverage(target_name, bed, data, bam_file)
        if dedupped:
            out["Preseq_unique_count"] = readstats.number_of_mapped_reads(
                data, bam_file, keep_dups=False, bed_file=bed, target_name=target_name)

        # Counting average on-target alignment length, based on the equation:
        #    avg depth ~~ num (unique) on-target alignments * avg on-target aln length / target size
        total_alignments = out.get("Preseq_unique_count") or out["Preseq_read_count"]
        out["Preseq_read_length"] = ontrg_unique_depth * out["Preseq_genome_size"] // total_alignments

    else:  # WGS
        out["Preseq_genome_size"] = sum([c.size for c in ref.file_contigs(dd.get_ref_file(data), data["config"])])
        out["Preseq_read_count"] = int(samtools_stats["Total_reads"])
        out["Preseq_read_length"] = int(samtools_stats["Average_read_length"])
        if dedupped:
            out["Preseq_unique_count"] = out["Preseq_read_count"] - int(samtools_stats["Duplicates"])

    return out
