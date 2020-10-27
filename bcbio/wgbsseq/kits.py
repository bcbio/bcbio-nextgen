from collections import namedtuple

Kit = namedtuple('Kit', 'name clip_r1_5 clip_r1_3 clip_r2_5 clip_r2_3 is_directional')

_KITS = [
    Kit("truseq", 8, 8, 8, 8, False),
    Kit("accelngs", 0, 19, 19, 0, True),
    Kit("nebemseq", 0, 0, 0, 0, True)
]

KITS = {x.name: x for x in _KITS}
SUPPORTED_KITS = {x.name for x in _KITS}
