from collections import namedtuple

Kit = namedtuple('Kit', 'name clip_r1_5 clip_r1_3 clip_r2_5 clip_r2_3')

_KITS = [
    Kit("truseq", 8, 8, 8, 8),
    Kit("accelngs", 0, 15, 15, 0)
]

KITS = {x.name: x for x in _KITS}
SUPPORTED_KITS = {x.name for x in _KITS}
