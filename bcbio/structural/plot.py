"""Provide plots of structural variations to manually validate results.

Uses existing plots from CNVkit along with custom plotting of coverage to
provide the ability to quickly validate and explore predicted structural
variants.
"""
from bcbio.pipeline import datadict as dd
from bcbio.variation import vcfutils

def _sort_by_type(x):
    """Simple prioritization to identify 'lead' items within a batch.
    """
    if vcfutils.get_paired_phenotype(x) == "tumor":
        priority = 0
    else:
        priority = 1
    return [priority, dd.get_sample_name(x)]

def by_regions(items):
    """Plot for a union set of combined ensemble regions across all of the data items.
    """
    import pybedtools
    items = sorted(items, key=_sort_by_type)
    calls = []
    for data in items:
        for sv in data["sv"]:
            if sv["variantcaller"] == "sv-ensemble":
                calls.append(sv)
    if len(calls) > 0:
        # Merge SV calls into a union set
        # Make summary level plots
        # Add plots to SV information for lead item
        pass
    print [x["description"] for x in items]
    return items
