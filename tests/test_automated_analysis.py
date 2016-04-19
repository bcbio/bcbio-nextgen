"""This directory is setup with configurations to run the main functional test.

It exercises a full analysis pipeline on a smaller subset of data.
"""
import os
import subprocess
import unittest
import shutil
import contextlib
import collections
import functools

from nose import SkipTest
from nose.plugins.attrib import attr
import yaml

from bcbio.pipeline.config_utils import load_system_config

@contextlib.contextmanager
def make_workdir():
    remove_old_dir = True
    #remove_old_dir = False
    dirname = os.path.join(os.path.dirname(__file__), "test_automated_output")
    if remove_old_dir:
        if os.path.exists(dirname):
            shutil.rmtree(dirname)
        os.makedirs(dirname)
    orig_dir = os.getcwd()
    try:
        os.chdir(dirname)
        yield dirname
    finally:
        os.chdir(orig_dir)

def expected_failure(test):
    """Small decorator to mark tests as expected failure.
    Useful for tests that are work-in-progress.
    """
    @functools.wraps(test)
    def inner(*args, **kwargs):
        try:
            test(*args, **kwargs)
        except Exception:
            raise SkipTest
        else:
            raise AssertionError('Failure expected')
    return inner

def get_post_process_yaml(data_dir, workdir):
    """Prepare a bcbio_system YAML file pointing to test data.
    """
    try:
        from bcbiovm.docker.defaults import get_datadir
        datadir = get_datadir()
        system = os.path.join(datadir, "galaxy", "bcbio_system.yaml") if datadir else None
    except ImportError:
        system = None
    if system is None or not os.path.exists(system):
        try:
            _, system = load_system_config("bcbio_system.yaml")
        except ValueError:
            system = None
    if system is None or not os.path.exists(system):
        system = os.path.join(data_dir, "post_process-sample.yaml")
    # create local config pointing to reduced genomes
    test_system = os.path.join(workdir, "bcbio_system.yaml")
    with open(system) as in_handle:
        config = yaml.load(in_handle)
        config["galaxy_config"] = os.path.join(data_dir, "universe_wsgi.ini")
        with open(test_system, "w") as out_handle:
            yaml.dump(config, out_handle)
    return test_system

class AutomatedAnalysisTest(unittest.TestCase):
    """Setup a full automated analysis and run the pipeline.
    """
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "data", "automated")

    def _install_test_files(self, data_dir):
        """Download required sequence and reference files.
        """
        DlInfo = collections.namedtuple("DlInfo", "fname dirname version")
        download_data = [DlInfo("110106_FC70BUKAAXX.tar.gz", None, None),
                         DlInfo("genomes_automated_test.tar.gz", "genomes", 29),
                         DlInfo("110907_ERP000591.tar.gz", None, None),
                         DlInfo("100326_FC6107FAAXX.tar.gz", None, 10),
                         DlInfo("tcga_benchmark.tar.gz", None, 3),
                         DlInfo("singlecell-rnaseq-test-data.tar.gz", "Harvard-inDrop", 1)]
        for dl in download_data:
            url = "http://chapmanb.s3.amazonaws.com/{fname}".format(fname=dl.fname)
            dirname = os.path.join(data_dir, os.pardir,
                                   dl.fname.replace(".tar.gz", "") if dl.dirname is None
                                   else dl.dirname)
            if os.path.exists(dirname) and dl.version is not None:
                version_file = os.path.join(dirname, "VERSION")
                is_old = True
                if os.path.exists(version_file):
                    with open(version_file) as in_handle:
                        version = int(in_handle.read())
                    is_old = version < dl.version
                if is_old:
                    shutil.rmtree(dirname)
            if not os.path.exists(dirname):
                self._download_to_dir(url, dirname)

    def _download_to_dir(self, url, dirname):
        print dirname
        cl = ["wget", url]
        subprocess.check_call(cl)
        cl = ["tar", "-xzvpf", os.path.basename(url)]
        subprocess.check_call(cl)
        shutil.move(os.path.basename(dirname), dirname)
        os.remove(os.path.basename(url))

    @attr(speed=3)
    def IGNOREtest_3_full_pipeline(self):
        """Run full automated analysis pipeline with multiplexing.

        XXX Multiplexing not supporting in latest versions.
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, os.pardir, "110106_FC70BUKAAXX"),
                  os.path.join(self.data_dir, "run_info.yaml")]
            subprocess.check_call(cl)

    @attr(speed=3)
    def IGNOREtest_4_empty_fastq(self):
        """Handle analysis of empty fastq inputs from failed runs.

        XXX Multiplexing not supporting in latest versions.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, os.pardir, "110221_empty_FC12345AAXX"),
                  os.path.join(self.data_dir, "run_info-empty.yaml")]
            subprocess.check_call(cl)

    @attr(stranded=True)
    @attr(rnaseq=True)
    def test_2_stranded(self):
        """Run an RNA-seq analysis with TopHat and generate gene-level counts.
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, os.pardir, "test_stranded"),
                  os.path.join(self.data_dir, "run_info-stranded.yaml")]
            subprocess.check_call(cl)

    @attr(rnaseq=True)
    @attr(tophat=True)
    def test_2_rnaseq(self):
        """Run an RNA-seq analysis with TopHat and generate gene-level counts.
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, os.pardir, "110907_ERP000591"),
                  os.path.join(self.data_dir, "run_info-rnaseq.yaml")]
            subprocess.check_call(cl)

    @attr(fusion=True)
    def test_2_fusion(self):
        """Run an RNA-seq analysis and test fusion genes
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, os.pardir, "test_fusion"),
                  os.path.join(self.data_dir, "run_info-fusion.yaml")]
            subprocess.check_call(cl)

    @attr(rnaseq=True)
    @attr(rnaseq_standard=True)
    @attr(star=True)
    def test_2_star(self):
        """Run an RNA-seq analysis with STAR and generate gene-level counts.
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, os.pardir, "110907_ERP000591"),
                  os.path.join(self.data_dir, "run_info-star.yaml")]
            subprocess.check_call(cl)

    @attr(rnaseq=True)
    @attr(fastrnaseq=True)
    def test_2_fastrnaseq(self):
        """Run a fast RNA-seq analysis
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, os.pardir, "110907_ERP000591"),
                  os.path.join(self.data_dir, "run_info-fastrnaseq.yaml")]
            subprocess.check_call(cl)

    # XXX Turned off until umis library installed via conda
    @expected_failure
    @attr(rnaseq=True)
    @attr(scrnaseq=True)
    def test_2_scrnaseq(self):
        """Run a single-cell RNA-seq analysis
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, os.pardir, "Harvard-inDrop"),
                  os.path.join(self.data_dir, "run_info-scrnaseq.yaml")]
            subprocess.check_call(cl)

    @attr(rnaseq=True)
    @attr(rnaseq_standard=True)
    @attr(hisat2=True)
    def test_2_hisat2(self):
        """Run an RNA-seq analysis with hisat2 and generate gene-level counts.
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, os.pardir, "110907_ERP000591"),
                  os.path.join(self.data_dir, "run_info-hisat2.yaml")]
            subprocess.check_call(cl)

    @attr(explant=True)
    @attr(singleend=True)
    @attr(rnaseq=True)
    def test_explant(self):
        """
        Run an explant RNA-seq analysis with TopHat and generate gene-level counts.
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, os.pardir, "1_explant"),
                  os.path.join(self.data_dir, "run_info-explant.yaml")]
            subprocess.check_call(cl)

    @attr(srnaseq=True)
    @attr(srnaseq_star=True)
    def test_srnaseq_star(self):
        """Run an sRNA-seq analysis.
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, "run_info-srnaseq_star.yaml")]
            subprocess.check_call(cl)

    @attr(srnaseq=True)
    @attr(srnaseq_bowtie=True)
    def test_srnaseq_bowtie(self):
        """Run an sRNA-seq analysis.
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, "run_info-srnaseq_bowtie.yaml")]
            subprocess.check_call(cl)

    @attr(chipseq=True)
    def test_chipseq(self):
        """
        Run a chip-seq alignment with Bowtie2
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, os.pardir, "test_chipseq"),
                  os.path.join(self.data_dir, "run_info-chipseq.yaml")]
            subprocess.check_call(cl)

    @attr(speed=1)
    @attr(ensemble=True)
    def test_1_variantcall(self):
        """Test variant calling with GATK pipeline.
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, os.pardir, "100326_FC6107FAAXX"),
                  os.path.join(self.data_dir, "run_info-variantcall.yaml")]
            subprocess.check_call(cl)

    @attr(speed=1)
    @attr(devel=True)
    def test_5_bam(self):
        """Allow BAM files as input to pipeline.
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, "run_info-bam.yaml")]
            subprocess.check_call(cl)

    @attr(speed=2)
    def test_6_bamclean(self):
        """Clean problem BAM input files that do not require alignment.
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, os.pardir, "100326_FC6107FAAXX"),
                  os.path.join(self.data_dir, "run_info-bamclean.yaml")]
            subprocess.check_call(cl)

    @attr(speed=2)
    @attr(cancer=True)
    @attr(cancermulti=True)
    def test_7_cancer(self):
        """Test paired tumor-normal calling using multiple calling approaches: MuTect, VarScan, FreeBayes.
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, "run_info-cancer.yaml")]
            subprocess.check_call(cl)

    @attr(cancer=True)
    @attr(cancerpanel=True)
    def test_7_cancer_nonormal(self):
        """Test cancer calling without normal samples or with normal VCF panels.
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, "run_info-cancer2.yaml")]
            subprocess.check_call(cl)

    @attr(speed=1)
    @attr(template=True)
    def test_8_template(self):
        """Create a project template from input files and metadata configuration.
        """
        self._install_test_files(self.data_dir)
        fc_dir = os.path.join(self.data_dir, os.pardir, "100326_FC6107FAAXX")
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py", "-w", "template", "--only-metadata",
                  "freebayes-variant",
                  os.path.join(fc_dir, "100326.csv"),
                  os.path.join(fc_dir, "7_100326_FC6107FAAXX_1_fastq.txt"),
                  os.path.join(fc_dir, "7_100326_FC6107FAAXX_2_fastq.txt"),
                  os.path.join(fc_dir, "8_100326_FC6107FAAXX.bam")]
            subprocess.check_call(cl)

    @attr(joint=True)
    def test_9_joint(self):
        """Perform joint calling/backfilling/squaring off following variant calling.
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_nextgen.py",
                  get_post_process_yaml(self.data_dir, workdir),
                  os.path.join(self.data_dir, "run_info-joint.yaml")]
            subprocess.check_call(cl)

    @attr(docker=True)
    def test_docker(self):
        """Run an analysis with code and tools inside a docker container.

        Requires https://github.com/chapmanb/bcbio-nextgen-vm
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_vm.py",
                  "--datadir=%s" % self.data_dir,
                  "run",
                  "--systemconfig=%s" % get_post_process_yaml(self.data_dir, workdir),
                  "--fcdir=%s" % os.path.join(self.data_dir, os.pardir, "100326_FC6107FAAXX"),
                  os.path.join(self.data_dir, "run_info-bam.yaml")]
            subprocess.check_call(cl)

    @attr(docker_ipython=True)
    def test_docker_ipython(self):
        """Run an analysis with code and tools inside a docker container, driven via IPython.

        Requires https://github.com/chapmanb/bcbio-nextgen-vm
        """
        self._install_test_files(self.data_dir)
        with make_workdir() as workdir:
            cl = ["bcbio_vm.py",
                  "--datadir=%s" % self.data_dir,
                  "ipython",
                  "--systemconfig=%s" % get_post_process_yaml(self.data_dir, workdir),
                  "--fcdir=%s" % os.path.join(self.data_dir, os.pardir, "100326_FC6107FAAXX"),
                  os.path.join(self.data_dir, "run_info-bam.yaml"),
                  "lsf", "localrun"]
            subprocess.check_call(cl)

class CWLTest(unittest.TestCase):
    """ Run simple CWL workflows.

    Requires https://github.com/chapmanb/bcbio-nextgen-vm
    """
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "data", "automated")

    @attr(speed=2)
    @attr(cwl=True)
    @attr(cwl_local=True)
    def test_1_cwl_local(self):
        """Create a common workflow language description and run on local installation.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_vm.py", "cwl", "../data/automated/run_info-cwl.yaml",
                  "--systemconfig", get_post_process_yaml(self.data_dir, workdir)]
            subprocess.check_call(cl)
            out_base = "run_info-cwl-workflow/main-run_info-cwl"
            cl = ["cwltool", "--verbose", "--preserve-environment", "PATH", "HOME", "--no-container",
                  out_base + ".cwl", out_base + "-samples.json"]
            subprocess.check_call(cl)
            print
            print "To run with a CWL tool, cd test_automated_output and:"
            print " ".join(cl)

    @attr(speed=2)
    @attr(cwl=True)
    @attr(cwl_docker=True)
    def test_2_cwl_docker(self):
        """Create a common workflow language description and run on a Docker installation.
        """
        with make_workdir() as workdir:
            cl = ["bcbio_vm.py", "cwl", "../data/automated/run_info-cwl.yaml",
                  "--systemconfig", get_post_process_yaml(self.data_dir, workdir)]
            subprocess.check_call(cl)
            out_base = "run_info-cwl-workflow/main-run_info-cwl"
            cl = ["cwltool", "--verbose", out_base + ".cwl", out_base + "-samples.json"]
            subprocess.check_call(cl)
            print
            print "To run with a CWL tool, cd test_automated_output and:"
            print " ".join(cl)
