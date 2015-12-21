"""Inspector that examines changes in the world object and output directory.

Work toward preparing bcbio to enumerate all inputs and outputs for CWL.
"""
import os

import toolz as tz
import yaml

from bcbio import utils
from bcbio.pipeline import datadict as dd


class WorldWatcher:
    """Watch changes in the world and output directory and report.

    Used to create input files we can feed into CWL creation about
    the changed state of the world.
    """
    def __init__(self, work_dir, is_on=True):
        self._work_dir = work_dir
        self._is_on = is_on
        if not self._is_on:
            return
        self._out_dir = utils.safe_makedir(os.path.join(work_dir, "world2cwl"))
        self._lworld = {}
        self._lfiles = set([])

    def _find_files(self):
        out = []
        for (dir, _, files) in os.walk(self._work_dir):
            out += [os.path.join(dir, f).replace(self._work_dir + "/", "") for f in files]
        return set(out)

    def _items_to_world(self, items):
        world = {}
        for item in items:
            assert len(item) == 1
            world[dd.get_sample_name(item[0])] = item[0]
        return world

    def _compare_dicts(self, orig, new, ns):
        out = {}
        for key, val in new.items():
            nskey = ns + [key]
            orig_val = tz.get_in([key], orig)
            if isinstance(val, dict) and isinstance(orig_val, dict):
                for nkey, nval in self._compare_dicts(orig_val or {}, val or {}, nskey).items():
                    out = tz.update_in(out, [nkey], lambda x: nval)
            elif val != orig_val:
                out = tz.update_in(out, nskey, lambda x: val)
        return out

    def initialize(self, world):
        if not self._is_on:
            return
        self._lfiles = self._find_files()
        self._lworld = self._items_to_world(world)

    def report(self, step, world):
        if not self._is_on:
            return
        new_files = self._find_files()
        file_changes = new_files - self._lfiles
        self._lfiles = new_files
        world_changes = self._compare_dicts(self._lworld, self._items_to_world(world), [])
        self._lworld = world
        out_file = os.path.join(self._out_dir, "%s.yaml" % step)
        with open(out_file, "w") as out_handle:
            out = {"files": file_changes,
                   "world": world_changes}
            yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)