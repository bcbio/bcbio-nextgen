"""Subset the genome into standard sets of regions surrounding transcripts.

Provides a central place to bin the genome into smaller transcript-based regions
for structural variant calling and prioritization.
"""
import itertools
import os

import pybedtools
import toolz as tz

from bcbio import utils
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import bedutils

def get_sv_bed(data, method=None, out_dir=None):
    """Retrieve a BED file of regions for SV and heterogeneity calling using the provided method.

    method choices:
      - exons: Raw BED file of exon regions
      - transcripts: Full collapsed regions with the min and max of each transcript.
      - transcriptsXXXX: Collapsed regions around transcripts with a window size of
        XXXX.
      - A custom BED file of regions
    """
    if method is None:
        method = tz.get_in(["config", "algorithm", "sv_regions"], data)
    gene_file = dd.get_gene_bed(data)
    if not gene_file or not method:
        return None
    elif os.path.isfile(method):
        return method
    elif method == "exons":
        return gene_file
    elif method.startswith("transcripts"):
        window = method.split("transcripts")[-1]
        window = int(float(window)) if window else 0
        return _collapse_transcripts(gene_file, window, data, out_dir)
    else:
        raise ValueError("Unexpected transcript retrieval method: %s" % method)

def _collapse_transcripts(in_file, window, data, out_dir):
    """Collapse transcripts into min/max coordinates and optionally add windows.
    """
    if out_dir is None:
        out_dir = os.path.dirname(in_file)
    out_file = os.path.join(out_dir,
                            "%s-transcripts_w%s.bed" % (os.path.splitext(os.path.basename(in_file))[0],
                                                        window))
    chrom_sizes = {}
    for contig in ref.file_contigs(dd.get_ref_file(data), data["config"]):
        chrom_sizes[contig.name] = contig.size
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(data, out_file) as tx_out_file:
            prep_file = "%s-sortprep%s" % os.path.splitext(tx_out_file)
            sort_cmd = bedutils.get_sort_cmd()
            cmd = "{sort_cmd} -k4,4 -k1,1 {in_file} > {prep_file}"
            do.run(cmd.format(**locals()), "Sort BED file by transcript name")
            with open(tx_out_file, "w") as out_handle:
                # Work around for segmentation fault issue with groupby
                # https://github.com/daler/pybedtools/issues/131#issuecomment-89832476
                x = pybedtools.BedTool(prep_file)
                def gen():
                    for r in x:
                        yield r
                for name, rs in itertools.groupby(gen(), lambda r: (r.name, r.chrom)):
                    rs = list(rs)
                    r = rs[0]
                    for gcoords in _group_coords(rs):
                        min_pos = max(min(gcoords) - window, 0)
                        max_pos = min(max(gcoords) + window, chrom_sizes[r.chrom])
                        out_handle.write("%s\t%s\t%s\t%s\n" % (r.chrom, min_pos, max_pos, r.name))
    return bedutils.sort_merge(out_file, data)

def _group_coords(rs):
    """Organize coordinate regions into groups for each transcript.

    Avoids collapsing very large introns or repetitive genes spread across
    the chromosome by limiting the intron size to 100kb for creating a single transcript
    """
    max_intron_size = 1e5
    coords = []
    for r in rs:
        coords.append(r.start)
        coords.append(r.end)
    coord_groups = []
    cur_group = []
    for coord in sorted(coords):
        if not cur_group or coord - cur_group[-1] < max_intron_size:
            cur_group.append(coord)
        else:
            coord_groups.append(cur_group)
            cur_group = [coord]
    if cur_group:
        coord_groups.append(cur_group)
    return coord_groups
