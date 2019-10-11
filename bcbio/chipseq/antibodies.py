#from dataclasses import dataclass
from collections import namedtuple

VALID_PEAKTYPES = ["narrow", "broad"]

# @dataclass
# class Antibody:
#     """
#     ChIP-seq antibody
#     """
#     name: str
#     # call narrow or broad peaks
#     peaktype: str
#     # remove duplicates?
#     rmdup: bool = True

#     def __post_init__(self):
#         if self.peaktype not in VALID_PEAKTYPES:
#             raise TypeError(f"peaktype {self.peatktype} is not one of {VALID_PEAKTYPES}")

Antibody = namedtuple('Antibody', 'name peaktype rmdup')

_ANTIBODIES = [
    Antibody("h3f3a", "broad", True),
    Antibody("h3k27me3", "broad", True),
    Antibody("h3k36me3", "broad", True),
    Antibody("h3k4me1", "broad", True),
    Antibody("h3k79me2", "broad", True),
    Antibody("h3k79me3", "broad", True),
    Antibody("h3k9me1", "broad", True),
    Antibody("h3k9me2", "broad", True),
    Antibody("h4k20me1", "broad", True),
    Antibody("h2afz", "narrow", True),
    Antibody("h3ac", "narrow", True),
    Antibody("h3k27ac", "narrow", True),
    Antibody("h3k4me2", "narrow", True),
    Antibody("h3k4me3", "narrow", True),
    Antibody("h3k9ac", "narrow", True),
    Antibody("h3k9me3", "broad", False)
]

ANTIBODIES = {x.name: x for x in _ANTIBODIES}
SUPPORTED_ANTIBODIES = {x.name for x in _ANTIBODIES}
