"""
Create log files to be parsed by multiqc
"""

import os
import pandas as pd

from bcbio import utils
from bcbio.provenance.programs import get_version_manifest
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd

def run(bam_file, data, out_dir):
    """Create several log files"""
    _mirbase_stats(data, out_dir)
    _seqcluster_stats(data, out_dir)

def _mirbase_stats(data, out_dir):
    """Create stats from miraligner"""
    utils.safe_makedir(out_dir)
    mirbase_fn = data.get("seqbuster", None)
    if mirbase_fn:
        out_file = os.path.join(out_dir, "%s_bcbio_mirbase.txt" % dd.get_sample_name(data))
        _get_stats_from_miraligner(mirbase_fn, out_file, "seqbuster")
    mirdeep_fn = data.get("seqbuster_novel", None)
    if mirdeep_fn:
        out_file = os.path.join(out_dir, "%s_bcbio_mirdeeep2.txt" % dd.get_sample_name(data))
        _get_stats_from_miraligner(mirdeep_fn, out_file, "mirdeep2")

def _get_stats_from_miraligner(fn, out_file, name):
    df = pd.read_csv(fn, sep="\t", dtype={"mism": "string",
                                          "add": "string",
                                          "t5": "string",
                                          "t3": "string"},
                     na_values=["."])
    dfmirs = df[['mir', 'freq']].groupby(['mir']).count()
    df5 = df.loc[df.t5 != "0", ['mir', 't5']].groupby(['mir']).count()
    df3 = df.loc[df.t3 != "0", ['mir', 't3']].groupby(['mir']).count()
    dfadd = df.loc[df["add"] != "0", ['mir', 'add']].groupby(['mir']).count()
    dfmut = df.loc[df.mism != "0", ['mir', 'mism']].groupby(['mir']).count()
    if not utils.file_exists(out_file):
        version = get_version_manifest("seqbuster")
        with file_transaction(out_file) as tx_out:
            with open(tx_out, "w") as out_handle:
                print >>out_handle, "# stats {name}, version: {version}".format(**locals())
                print >>out_handle, ("mirs\t{mirs}\nisomirs\t{isomirs}").format(
                        mirs=len(dfmirs.index), isomirs=len(df.index))
                print >>out_handle, ("mirs_mutations\t{muts}\nmirs_additions\t{add}").format(
                        muts=len(dfmut.index), add=len(dfadd.index))
                print >>out_handle, ("mirs_5-trimming\t{t5}\nmirs_3-trimming\t{t3}").format(
                        t5=len(df5.index), t3=len(df3.index))
                print >>out_handle, ("iso_mutations\t{muts}\niso_additions\t{add}").format(
                        muts=sum(dfmut.mism), add=sum(dfadd["add"]))
                print >>out_handle, ("iso_5-trimming\t{t5}\niso_3-trimming\t{t3}").format(
                        t5=sum(df5.t5), t3=sum(df3.t3))
    return out_file

def _seqcluster_stats(data, out_dir):
    """Parse seqcluster output"""
    name = dd.get_sample_name(data)
    fn = data.get("seqcluster", {}).get("stat_file", None)
    if not fn:
        return None
    out_file = os.path.join(out_dir, "%s.txt" % name)
    df = pd.read_csv(fn, sep="\t", names = ["reads", "sample", "type"])
    df_sample = df[df["sample"] == name]
    df_sample.to_csv(out_file, sep="\t")

