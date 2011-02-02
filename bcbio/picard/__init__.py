"""Ease working with Broad's Picard BAM manipulation and analysis library.
"""
import os
import subprocess

import run as run_fns

class PicardRunner:
    """Simplify running Picard commands.
    """
    def __init__(self, picard_dir, max_memory="6g"):
        self._memory_args = []
        if max_memory:
            self._memory_args.append("-Xmx%s" % max_memory)
        self._picard_dir = picard_dir

    def run_fn(self, name, *args, **kwds):
        fn = getattr(run_fns, name)
        return fn(self, *args, **kwds)

    def run(self, command, options):
        options = ["%s=%s" % (x, y) for x, y in options]
        options.append("VALIDATION_STRINGENCY=SILENT")
        dist_file = self._get_jar(command)
        cl = ["java"] + self._memory_args +["-jar", dist_file] + options
        subprocess.check_call(cl)

    def run_gatk(self, params, tmp_dir=None):
        gatk_jar = self._get_jar("GenomeAnalysisTK")
        local_args = []
        if tmp_dir:
            local_args.append("-Djava.io.tmpdir=%s" % tmp_dir)
        cl = ["java"] + self._memory_args + local_args + \
                ["-jar", gatk_jar] + [str(x) for x in params]
        #print " ".join(cl)
        subprocess.check_call(cl)

    def _get_jar(self, command):
        """Retrieve the jar for running the specified command.
        """
        dirs = [self._picard_dir,
                os.path.join(self._picard_dir, "dist"),
                os.path.join(self._picard_dir, "Picard-private", "dist"),
                os.path.join(self._picard_dir, "GATK"),
                os.path.join(self._picard_dir, "GATK", "dist")]
        for dir_check in dirs:
            check_file = os.path.join(dir_check, "%s.jar" % command)
            if os.path.exists(check_file):
                return check_file
        raise ValueError("Could not find jar %s in %s" % (command, self._picard_dir))
