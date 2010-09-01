"""Ease working with Broad's Picard BAM manipulation and analysis library.
"""
import os
import subprocess

class PicardRunner:
    """Simplify running Picard commands.
    """
    def __init__(self, picard_dir, max_memory="4g"):
        self._memory_args = []
        if max_memory:
            self._memory_args.append("-Xmx%s" % max_memory)
        self._picard_dir = picard_dir

    def run(self, command, options):
        options = ["%s=%s" % (x, y) for x, y in options]
        options.append("VALIDATION_STRINGENCY=SILENT")
        dist_file = os.path.join(self._picard_dir, "dist", "%s.jar" % command)
        private_dist_file = os.path.join(self._picard_dir, "Picard-private",
                "dist", "%s.jar" % command)
        if not os.path.exists(dist_file):
            if os.path.exists(private_dist_file):
                dist_file = private_dist_file
            else:
                raise ValueError("Could not find jar %s in %s" % (command,
                    self._picard_dir))
        #cl = ["java", "-XX:-UseGCOverheadLimit", "-Xmx1800m", "-Xms1800m", "-jar", dist_file] + options
        cl = ["java"] + self._memory_args +["-jar", dist_file] + options
        subprocess.check_call(cl)

    def run_gatk(self, params, tmp_dir=None):
        gatk_jar = os.path.join(self._picard_dir, "GATK", "dist",
                "GenomeAnalysisTK.jar")
         #       "GATK-Picard.jar")
        local_args = []
        if tmp_dir:
            local_args.append("-Djava.io.tmpdir=%s" % tmp_dir)
        cl = ["java"] + self._memory_args + local_args + \
                ["-jar", gatk_jar] + [str(x) for x in params]
        #print " ".join(cl)
        subprocess.check_call(cl)
