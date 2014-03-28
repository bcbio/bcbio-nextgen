import unittest
from nose.plugins.attrib import attr
from bcbio.pipeline.main import *
from bcbio.pipeline.main import _get_pipeline, _pair_lanes_with_pipelines, _add_provenance
from bcbio.pipeline.run_info import organize, _run_info_from_yaml, add_reference_resources
from bcbio.pipeline.config_utils import update_w_custom
from bcbio.provenance.programs import _get_versions
from bcbio.provenance import versioncheck


class TestPipelineSetup(unittest.TestCase):
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data/")
        self.config_samples_dir = os.path.join(self.data_dir, "config_samples/")
        self.galaxy_dir = os.path.join(self.data_dir, "galaxy/")
        self.bin_dir = os.path.join(self.data_dir, "bin/")

    # main.run_main(config, ...) --> main._run_toplevel(config, ...), which calls:
    #  1) run_info = organize(config)
    #  2) pipelines = _pair_lanes_with_pipelines(run_info)
    #  3) pipeline_items = _add_provenance(pipeline_items, ...)
    #  4) versioncheck.testall(pipeline_items)
    #  5) pipeline.run(config, ...)

    # 1) organize(config) calls:
    #  1) run_info = _run_info_from_yaml(config)
    #  2) run_info["config"] = update_w_custom(config, run_info)
    #  3) add_reference_resources(run_info)
    def test_organize(self):
        dirs = {
            "config": None,
            "fastq": None,
            "flowcell": None,
            "galaxy": self.galaxy_dir,
            "work": None,
        }
        config = {}
        run_info_yaml = os.path.join(self.config_samples_dir, "pipeline.yaml")
        run_info = organize(dirs, config, run_info_yaml)

    def test_run_info_from_yaml(self):
        fc_dir = None
        run_info_yaml = os.path.join(self.config_samples_dir, "pipeline.yaml")
        config = {}
        run_info = _run_info_from_yaml(fc_dir, run_info_yaml, config)

    def test_update_w_custom(self):
        config = {}
        lane_info = {}
        config = update_w_custom(config, lane_info)

    def test_add_reference_resources(self):
        run_info = {
            "config": {
                "algorithm": {}
            },
            "genome_build": "GRCh37",
            "dirs": {
                "galaxy": self.galaxy_dir
            }
        }
        run_info = add_reference_resources(run_info)

        reference = {
            'snpeff': {},
            'fasta': {
                'base': os.path.join(self.galaxy_dir, 'tool-data/human_g1k_v37'),
                'indexes': [
                    os.path.join(self.galaxy_dir, 'tool-data/human_g1k_v37.dict.gz'),
                    os.path.join(self.galaxy_dir, 'tool-data/human_g1k_v37.fasta.fai.gz'),
                    os.path.join(self.galaxy_dir, 'tool-data/human_g1k_v37.fasta.gz')
                ]
            }
        }
        sam_ref = os.path.join(self.galaxy_dir, 'tool-data/human_g1k_v37')
        genome_resources = {'version': 8}

        assert run_info.get("reference") == reference
        assert run_info.get("sam_ref") == sam_ref
        assert run_info.get("genome_resources") == genome_resources


    # 2) main._pair_lanes_with_pipelines(run_info) calls:
    #  1) main._get_pipeline(run_info)
    SUPPORTED_PIPELINES = {x.name.lower(): x for x in
                           utils.itersubclasses(AbstractPipeline)}

    def test_pair_lanes_with_pipelines(self):
        run_info = [
            {"analysis": "variant2"},
            {"analysis": "chip-seq"},
        ]
        pipelines = _pair_lanes_with_pipelines(run_info)

    def test_supported_pipelines(self):
        print(self.SUPPORTED_PIPELINES.keys())
        assert len(self.SUPPORTED_PIPELINES) > 0

    def test_get_pipeline(self):
        item = {"analysis": "variant2"}
        assert _get_pipeline(item) == self.SUPPORTED_PIPELINES["variant2"]

    # 3) pipeline_items = _add_provenance(pipeline_items, ...)
    #  1) programs.write_versions(...)
    #  2) system.write_info(...)
    def test_add_provenance(self):
        items = [{"description": "", "config": {"resources": {}}}]
        dirs = {
            "config": None,
            "fastq": None,
            "flowcell": None,
            "galaxy": None,
            "work": None,
        }
        parallel = {"type": "local"}
        config = {
            "resources": {
                "picard": {
                    "dir": self.bin_dir
                },
                "gatk": {
                    "dir": self.bin_dir
                },
                "mutect": {
                    "dir": self.bin_dir
                },
            }
        }
        pipeline_items = _add_provenance(items, dirs, parallel, config)

    # comes down to the _alt_progs using _broad_versioner
    def test_get_versions(self):
        config = {
            "resources": {
                "picard": {
                    "dir": self.bin_dir
                },
                "gatk": {
                    "dir": self.bin_dir
                },
                "mutect": {
                    "dir": self.bin_dir
                },
            }
        }
        versions = _get_versions(config)

    # 4) versioncheck.testall(pipeline_items)
    @attr('current')
    def test_versioncheck(self):
        config = {
            "resources": {
                "samtools": {
                    "cmd": os.path.join(self.bin_dir, "samtools")
                },
                "picard": {
                    "dir": self.bin_dir
                },
                "gatk": {
                    "dir": self.bin_dir
                },
                "mutect": {
                    "dir": self.bin_dir
                },
            }
        }
        pipeline_items = [[{"config": config}]]
        versioncheck.testall(pipeline_items)

