from dataclasses import dataclass

VALID_PEAKTYPES = ["narrow", "broad"]

@dataclass
class Antibody:
    """
    ChIP-seq antibody
    """
    name: str
    # call narrow or broad peaks
    peaktype: str
    # remove duplicates?
    rmdup: bool = True

    def __post_init__(self):
        if self.peaktype not in VALID_PEAKTYPES:
            raise TypeError(f"peaktype {self.peatktype} is not one of {VALID_PEAKTYPES}")

_ANTIBODIES = [
    Antibody("h3f3a", "broad"),
    Antibody("h3k27me3", "broad"),
    Antibody("h3k36me3", "broad"),
    Antibody("h3k4me1", "broad"),
    Antibody("h3k79me2", "broad"),
    Antibody("h3k79me3", "broad"),
    Antibody("h3k9me1", "broad"),
    Antibody("h3k9me2", "broad"),
    Antibody("h4k20me1", "broad"),
    Antibody("h2afz", "narrow"),
    Antibody("h3ac", "narrow"),
    Antibody("h3k27ac", "narrow"),
    Antibody("h3k4me2", "narrow"),
    Antibody("h3k4me3", "narrow"),
    Antibody("h3k9ac", "narrow"),
    Antibody("h3k9me3", "broad", False)
]

ANTIBODIES = {x.name: x for x in _ANTIBODIES}
SUPPORTED_ANTIBODIES = {x.name for x in _ANTIBODIES}
