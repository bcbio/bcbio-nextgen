"""Next-gen variant detection and evaluation with GATK and SnpEff.
"""
from bcbio.log import logger
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
                                              data)
        logger.info("Calculating variation effects for %s" % cur_name)
        ann_vrn_file = effects.snpeff_effects(data)
        if ann_vrn_file:
            data["vrn_file"] = ann_vrn_file
    return [[data]]
