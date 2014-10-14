"""Prepare a workflow for running on AWS using STORMSeq as a front end.

http://www.stormseq.org/
"""
import argparse
import json
import os

import yaml

from bcbio import utils
from bcbio.upload import s3
from bcbio.workflow import xprize

def parse_args(args):
    parser = xprize.HelpArgParser(description="Run STORMSeq processing on AWS")
    parser.add_argument("config_file", help="JSON configuration file with form parameters")
    parser.add_argument("base_dir", help="Base directory to process in")
    parser.add_argument("bcbio_config_file", help="bcbio system YAML config")
    args = parser.parse_args(args)
    return args

def _get_s3_files(local_dir, file_info, params):
    """Retrieve s3 files to local directory, handling STORMSeq inputs.
    """
    assert len(file_info) == 1
    files = file_info.values()[0]
    fnames = []
    for k in ["1", "2"]:
        if files[k] not in fnames:
            fnames.append(files[k])
    out = []
    for fname in fnames:
        bucket, key = fname.replace("s3://", "").split("/", 1)
        if params["access_key_id"] == "TEST":
            out.append(os.path.join(local_dir, os.path.basename(key)))
        else:
            out.append(s3.get_file(local_dir, bucket, key, params))
    return out

def setup(args):
    configdir = utils.safe_makedir(os.path.join(args.base_dir, "config"))
    inputdir = utils.safe_makedir(os.path.join(args.base_dir, "inputs"))
    workdir = utils.safe_makedir(os.path.join(args.base_dir, "work"))
    finaldir = utils.safe_makedir(os.path.join(args.base_dir, "ready"))
    out_config_file = os.path.join(configdir, "%s.yaml" %
                                   os.path.splitext(os.path.basename(args.config_file))[0])
    with open(args.config_file) as in_handle:
        ss_config = json.load(in_handle)
        ss_params = ss_config["parameters"]
    out = {"fc_date": xprize.get_fc_date(out_config_file),
           "fc_name": ss_config["sample"],
           "upload": {"dir": finaldir,
                      "method": "s3",
                      "bucket": ss_params["s3_bucket"],
                      "access_key_id": ss_params["access_key_id"],
                      "secret_access_key": ss_params["secret_access_key"]},
            "details": [{
                "files": _get_s3_files(inputdir, ss_config["files"], ss_params),
                "lane": 1,
                "description": ss_params["sample"],
                "analysis": "variant",
                "genome_build": ss_params["genome_version"],
                "algorithm": {
                    "aligner": ss_params["alignment_pipeline"],
                    "variantcaller": ss_params["calling_pipeline"],
                    "quality_format": "Standard",
                    "coverage_interval": "genome" if ss_params["data_type"] == "data_wgs" else "exome",
                    }}]}
    with open(out_config_file, "w") as out_handle:
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)

    return workdir, {"config_file": args.bcbio_config_file,
                     "run_info_yaml": out_config_file}
