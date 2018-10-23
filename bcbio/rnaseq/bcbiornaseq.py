import os
import shutil
import toolz as tz
from string import Template
from bcbio.utils import file_exists, Rscript_cmd, safe_makedir, chdir
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.provenance import do
from bcbio.pipeline import datadict as dd

def make_bcbiornaseq_object(data):
    """
    load the initial bcb.rda object using bcbioRNASeq
    """
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
        memoize_write_file(loadstring, tmp_file)
    rcmd = Rscript_cmd()
    with chdir(report_dir):
        do.run([rcmd, "--no-environ", r_file], "Loading bcbioRNASeq object.")
    make_quality_report(data)
    return data

def make_quality_report(data):
    """
    create and render the bcbioRNASeq quality report
    """
    if "bcbiornaseq" not in dd.get_tools_on(data):
        return data
    upload_dir = tz.get_in(("upload", "dir"), data)
    report_dir = os.path.join(upload_dir, "bcbioRNASeq")
    safe_makedir(report_dir)
    quality_rmd = os.path.join(report_dir, "quality_control.Rmd")
    quality_html = os.path.join(report_dir, "quality_control.html")
    quality_rmd = rmarkdown_draft(quality_rmd, "quality_control", "bcbioRNASeq")

def rmarkdown_draft(filename, template, package):
    """
    create a draft rmarkdown file from an installed template
    """
    draft_template = Template(
        'rmarkdown::draft("$filename", template="$template", package="$package", edit=FALSE)'
    )
    draft_string = draft_template.substitute(
        filename=filename, template=template, package=package)
    report_dir = os.path.dirname(filename)
    rcmd = Rscript_cmd()
    with chdir(report_dir):
        do.run([rcmd, "--no-environ", "-e", draft_string], "Creating bcbioRNASeq quality control template.")
        do.run(["sed", "-i", "s/YYYY-MM-DD\///g", filename], "Editing bcbioRNAseq quality control template.")
    render_rmarkdown_file(filename)
    return filename

def render_rmarkdown_file(filename):
    """
    render a rmarkdown file using the rmarkdown library
    """
    render_template = Template(
        'rmarkdown::render("$filename")'
    )
    render_string = render_template.substitute(
        filename=filename)
    report_dir = os.path.dirname(filename)
    rcmd = Rscript_cmd()
    with chdir(report_dir):
        do.run([rcmd, "--no-environ", "-e", render_string], "Rendering bcbioRNASeq quality control report.")
    return filename

def create_load_string(upload_dir, groups=None, organism=None):
    """
    create the code necessary to load the bcbioRNAseq object
    """
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
        load_bcbio = load_template.substitute(
            upload_dir=upload_dir, groups=groups, organism=organism)
    else:
        load_bcbio = load_noorganism_template.substitute(upload_dir=upload_dir,
                                                groups=groups)
    return ";\n".join([libraryline, load_bcbio, flatline, saveline])

def memoize_write_file(string, filename):
    if file_exists(filename):
        return filename
    with file_transaction(filename) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            out_handle.write(string)
    return filename

def _quotestring(string, double=True):
    """ escape quote a string """
    if double:
        return "\"" + string + "\""
    else:
        return "\'" + string + "\'"

def _list2Rlist(xs):
    """ convert a python list to an R list """
    if isinstance(xs, basestring):
        xs = [xs]
    rlist = ",".join([_quotestring(x) for x in xs])
    return "c(" + rlist + ")"
