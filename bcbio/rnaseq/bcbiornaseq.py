import os
import shutil
import toolz as tz
from string import Template
from bcbio.utils import file_exists, Rscript_cmd, safe_makedir, chdir
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.provenance import do
from bcbio.pipeline import datadict as dd

def make_bcbiornaseq_object(data):
    if "bcbiornaseq" not in dd.get_tools_on(data):
        return data
    upload_dir = tz.get_in(("upload", "dir"), data)
    report_dir = os.path.join(upload_dir, "bcbioRNASeq")
    safe_makedir(report_dir)
    organism = dd.get_bcbiornaseq(data).get("organism", None)
    groups = dd.get_bcbiornaseq(data).get("interesting_groups", None)
    loadstring = create_load_string(upload_dir, groups, organism)
    r_file = os.path.join(report_dir, "load_bcbioRNAseq.R")
    with file_transaction(r_file) as tmp_file:
        write_load_bcbiornaseq_file(loadstring, tmp_file)
    rcmd = Rscript_cmd()
    with chdir(report_dir):
        do.run([rcmd, r_file], "Loading bcbioRNASeq object.")
    return data

def create_load_string(upload_dir, groups=None, organism=None):
    libraryline = 'library(bcbioRNASeq)'
    load_template = Template(
        ('bcb <- bcbioRNASeq(uploadDir="$upload_dir",'
         'interestingGroups=$groups,'
         'organism="$organism")'))
    load_noorganism_template = Template(
        ('bcb <- bcbioRNASeq(uploadDir="$upload_dir",'
         'interestingGroups=$groups,'
         'organism=NULL)'))
    flatline = 'flat <- flatFiles(bcb)'
    saveline = 'saveData(bcb, flat, dir="data")'
    if groups:
        groups = _list2Rlist(groups)
    else:
        groups = _quotestring("sampleName")
    if organism:
        load_bcbio = load_template.substitute(upload_dir=upload_dir,
                                              groups=groups,
                                              organism=organism)
    else:
        load_bcbio = load_noorganism_template.substitute(upload_dir=upload_dir,
                                                groups=groups)
    return ";\n".join([libraryline, load_bcbio, flatline, saveline])

def write_load_bcbiornaseq_file(string, filename):
    if file_exists(filename):
        return filename
    with file_transaction(filename) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            out_handle.write(string)
    return filename

def _quotestring(string):
    """ escape quote a string """
    return "\"" + string + "\""

def _list2Rlist(xs):
    """ convert a python list to an R list """
    if isinstance(xs, basestring):
        xs = [xs]
    rlist = ",".join([_quotestring(x) for x in xs])
    return "c(" + rlist + ")"
