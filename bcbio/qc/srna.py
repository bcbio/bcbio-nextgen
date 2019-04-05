"""
Create log files to be parsed by multiqc
"""
from __future__ import print_function

import os
import pandas as pd

from bcbio import utils
from bcbio.provenance.programs import get_version_manifest
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd

def run(bam_file, data, out_dir):
    """Create several log files"""
    m = {"base": None, "secondary": []}
    m.update(_mirbase_stats(data, out_dir))
    m["secondary"].append(_seqcluster_stats(data, out_dir))

def _mirbase_stats(data, out_dir):
    """Create stats from miraligner"""
    utils.safe_makedir(out_dir)
    out_file = os.path.join(out_dir, "%s_bcbio_mirbase.txt" % dd.get_sample_name(data))
    out_file_novel = os.path.join(out_dir, "%s_bcbio_mirdeeep2.txt" % dd.get_sample_name(data))
    mirbase_fn = data.get("seqbuster", None)
    if mirbase_fn:
        _get_stats_from_miraligner(mirbase_fn, out_file, "seqbuster")
    mirdeep_fn = data.get("seqbuster_novel", None)
    if mirdeep_fn:
        _get_stats_from_miraligner(mirdeep_fn, out_file_novel, "mirdeep2")
    return {"base": out_file, "secondary": [out_file_novel]}

def _get_stats_from_miraligner(fn, out_file, name):
    df = pd.read_csv(fn, sep="\t", dtype={"mism": "str",
                                          "add": "str",
                                          "t5": "str",
                                          "t3": "str"},
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
                print(("# stats {name}, version: {version}").format(**locals()), file=out_handle)
                print(("mirs\t{mirs}\nisomirs\t{isomirs}").format(
                    mirs=len(dfmirs.index), isomirs=len(df.index)), file=out_handle)
                print(("mirs_mutations\t{muts}\nmirs_additions\t{add}").format(
                    muts=len(dfmut.index), add=len(dfadd.index)), file=out_handle)
                print(("mirs_5-trimming\t{t5}\nmirs_3-trimming\t{t3}").format(
                    t5=len(df5.index), t3=len(df3.index)), file=out_handle)
                print(("iso_mutations\t{muts}\niso_additions\t{add}").format(
                    muts=sum(dfmut.mism), add=sum(dfadd["add"])), file=out_handle)
                print(("iso_5-trimming\t{t5}\niso_3-trimming\t{t3}").format(
                    t5=sum(df5.t5), t3=sum(df3.t3)), file=out_handle)
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
    return out_file

