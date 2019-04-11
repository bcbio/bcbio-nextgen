"""Get log file from trimming step if file inside data["log_trimming"]"""

from bcbio import utils

def run(bam_file, data, dir_out):
    m = {}
    if "log_trimming" in data and utils.file_exists(data["log_trimming"]):
        with open(data["log_trimming"]) as inh:
            for line in inh:
                cols = line.strip().split(":")
                if line.find("Total reads processed") > -1:
                    m["reads_before_trimming"] = _get_value(cols[1])
                elif line.find("Reads with adapter") > -1:
                    m["read_with_adapter"] = _get_value(cols[1])
                elif line.find("Reads written") > -1:
                    m["read_pass_filter"] = _get_value(cols[1])
        return {"metrics": m, "base": data["log_trimming"]}

def _get_value(information):
    return int(information.strip().split()[0].strip().replace(",", ""))

