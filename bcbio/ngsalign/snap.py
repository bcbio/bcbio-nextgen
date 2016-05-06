"""Alignment with SNAP: http://snap.cs.berkeley.edu/
"""
import os

from bcbio import utils
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.ngsalign import alignprep, novoalign, postalign
from bcbio.provenance import do

def align(fastq_file, pair_file, index_dir, names, align_dir, data):
    """Perform piped alignment of fastq input files, generating sorted, deduplicated BAM.

    Pipes in input, handling paired and split inputs, using interleaving magic
    from: https://biowize.wordpress.com/2015/03/26/the-fastest-darn-fastq-decoupling-procedure-i-ever-done-seen/

    Then converts a tab delimited set of outputs into interleaved fastq.

    awk changes spaces to underscores since SNAP only takes the initial name.
    SNAP requires /1 and /2 at the end of read names. If these are not present
    in the initial fastq may need to expand awk code to do this.
    """
    out_file = os.path.join(align_dir, "{0}-sort.bam".format(dd.get_sample_name(data)))
    num_cores = data["config"]["algorithm"].get("num_cores", 1)
    resources = config_utils.get_resources("snap", data["config"])
    rg_info = novoalign.get_rg_info(names)

    if data.get("align_split"):
        final_file = out_file
        out_file, data = alignprep.setup_combine(final_file, data)
        fastq_file, pair_file = alignprep.split_namedpipe_cls(fastq_file, pair_file, data)
        fastq_file = fastq_file[2:-1]
        if pair_file:
            pair_file = pair_file[2:-1]
            stream_input = (r"paste <({fastq_file} | paste - - - -) "
                            r"<({pair_file} | paste - - - -) | "
                            r"""awk 'BEGIN {{FS="\t"; OFS="\n"}} """
                            r"""{{ """
                            r"""split($1, P1, " "); split($5, P5, " "); """
                            r"""if ($1 !~ /\/1$/) $1 = P1[1]"/1"; if ($5 !~ /\/2$/) $5 = P5[1]"/2"; """
                            r"""gsub(" ", "_", $1); gsub(" ", "_", $5); """
                            r"""print $1, $2, "+", $4, $5, $6, "+", $8}}' """)
        else:
            stream_input = fastq_file[2:-1]
    else:
        final_file = None
        assert fastq_file.endswith(".gz")
        if pair_file:
            stream_input = (r"paste <(zcat {fastq_file} | paste - - - -) "
                            r"<(zcat {pair_file} | paste - - - -) | "
                            r"""awk 'BEGIN {{FS="\t"; OFS="\n"}} """
                            r"""{{ """
                            r"""split($1, P1, " "); split($5, P5, " "); """
                            r"""if ($1 !~ /\/1$/) $1 = P1[1]"/1"; if ($5 !~ /\/2$/) $5 = P5[1]"/2"; """
                            r"""gsub(" ", "_", $1); gsub(" ", "_", $5); """
                            r"""print $1, $2, "+", $4, $5, $6, "+", $8}}' """)
        else:
            stream_input = "zcat {fastq_file}"

    pair_file = pair_file if pair_file else ""
    if not utils.file_exists(out_file) and (final_file is None or not utils.file_exists(final_file)):
        with postalign.tobam_cl(data, out_file, pair_file is not None) as (tobam_cl, tx_out_file):
            if pair_file:
                sub_cmd = "paired"
                input_cmd = "-pairedInterleavedFastq -"
            else:
                sub_cmd = "single"
                input_cmd = "-fastq -"
            stream_input = stream_input.format(**locals())
            tmp_dir = os.path.dirname(tx_out_file)
            cmd = ("export TMPDIR={tmp_dir} && {stream_input} | snap-aligner {sub_cmd} {index_dir} {input_cmd} "
                   "-R '{rg_info}' -t {num_cores} -M -o -sam - | ")
            do.run(cmd.format(**locals()) + tobam_cl, "SNAP alignment: %s" % names["sample"])
    data["work_bam"] = out_file
    return data

# Optional galaxy location file. Falls back on remap_index_fn if not found
galaxy_location_file = "snap_indices.loc"

def remap_index_fn(ref_file):
    """Map sequence references to snap reference directory, using standard layout.
    """
    snap_dir = os.path.join(os.path.dirname(ref_file), os.pardir, "snap")
    assert os.path.exists(snap_dir) and os.path.isdir(snap_dir), snap_dir
    return snap_dir
