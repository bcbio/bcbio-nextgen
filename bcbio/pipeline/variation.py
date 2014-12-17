"""Next-gen variant detection and evaluation with GATK and SnpEff.
"""
import toolz as tz

from bcbio.log import logger
from bcbio.variation.genotype import variant_filtration, get_variantcaller
from bcbio.variation import effects, prioritize

# ## Genotyping

def postprocess_variants(data):
    """Provide post-processing of variant calls: filtering and effects annotation.
    """
    cur_name = "%s, %s" % (data["name"][-1], get_variantcaller(data))
    logger.info("Finalizing variant calls: %s" % cur_name)
    if data.get("align_bam") and data.get("vrn_file"):
        logger.info("Calculating variation effects for %s" % cur_name)
        effect_todo = effects.get_type(data)
        if effect_todo:
            if effect_todo == "snpeff":
                ann_vrn_file = effects.snpeff_effects(data)
            elif effect_todo == "vep":
                ann_vrn_file = effects.run_vep(data)
            else:
                raise ValueError("Unexpected variant effects configuration: %s" % effect_todo)
            if ann_vrn_file:
                data["vrn_file"] = ann_vrn_file
        logger.info("Filtering for %s" % cur_name)
        data["vrn_file"] = variant_filtration(data["vrn_file"], data["sam_ref"],
                                              tz.get_in(("genome_resources", "variation"), data, {}),
                                              data)
        logger.info("Prioritization for %s" % cur_name)
        data["vrn_file"] = prioritize.handle_vcf_calls(data["vrn_file"], data)
    return [[data]]
