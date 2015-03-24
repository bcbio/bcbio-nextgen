"""Annotate structural variant calls with associated genes.
"""

def with_genes(svcalls):
    out = []
    for sv in svcalls:
        if sv["variantcaller"] == "sv-ensemble":
            # Do annotations
            pass
        out.append(sv)
    return out
