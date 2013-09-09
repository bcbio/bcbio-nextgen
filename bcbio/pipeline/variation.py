"""Next-gen variant detection and evaluation with GATK and SnpEff.
"""
import os
import json

from bcbio.log import logger
from bcbio.structural import hydra
from bcbio.variation.genotype import variant_filtration, get_variantcaller
from bcbio.variation import effects

# ## Genotyping

def postprocess_variants(data):
    """Provide post-processing of variant calls: filtering and effects annotation.
    """
    cur_name = "%s, %s" % (data["name"][-1], get_variantcaller(data))
    logger.info("Finalizing variant calls: %s" % cur_name)
    if data["work_bam"] and data.get("vrn_file"):
        data["vrn_file"] = variant_filtration(data["vrn_file"], data["sam_ref"],
                                              data["genome_resources"]["variation"],
                                              data["config"])
        logger.info("Calculating variation effects for %s" % cur_name)
        ann_vrn_file = effects.snpeff_effects(data)
        if ann_vrn_file:
            data["vrn_file"] = ann_vrn_file
    return [[data]]

# ## Structural variation

def detect_sv(data):
    """Detect structural variation for input sample.
    """
    sv_todo = data["config"]["algorithm"].get("sv_detection", None)
    if sv_todo is not None and data.get("fastq2"):
        if sv_todo == "hydra":
            sv_calls = hydra.detect_sv(data["work_bam"], data["genome_build"],
                                       data["dirs"], data["config"])
        else:
            raise ValueError("Unexpected structural variation method:{}".format(sv_todo))
    return [[data]]
