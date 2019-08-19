"""Utilities to manage processing flowcells and retrieving Galaxy stored info.
"""
import os
import glob
from six.moves import urllib, http_cookiejar
import json

def parse_dirname(fc_dir):
    """Parse the flow cell ID and date from a flow cell directory.
    """
    (_, fc_dir) = os.path.split(fc_dir)
    parts = fc_dir.split("_")
    name = None
    date = None
    for p in parts:
        if p.endswith(("XX", "xx", "XY", "X2", "X3")):
            name = p
        elif len(p) == 6:
            try:
                int(p)
                date = p
            except ValueError:
                pass
    if name is None or date is None:
        raise ValueError("Did not find flowcell name: %s" % fc_dir)
    return name, date

def get_qseq_dir(fc_dir):
    """Retrieve the qseq directory within Solexa flowcell output.
    """
    machine_bc = os.path.join(fc_dir, "Data", "Intensities", "BaseCalls")
    if os.path.exists(machine_bc):
        return machine_bc
    # otherwise assume we are in the qseq directory
    # XXX What other cases can we end up with here?
    else:
        return fc_dir

def get_fastq_dir(fc_dir):
    """Retrieve the fastq directory within Solexa flowcell output.
    """
    full_goat_bc = glob.glob(os.path.join(fc_dir, "Data", "*Firecrest*", "Bustard*"))
    bustard_bc = glob.glob(os.path.join(fc_dir, "Data", "Intensities", "*Bustard*"))
    machine_bc = os.path.join(fc_dir, "Data", "Intensities", "BaseCalls")
    if os.path.exists(machine_bc):
        return os.path.join(machine_bc, "fastq")
    elif len(full_goat_bc) > 0:
        return os.path.join(full_goat_bc[0], "fastq")
    elif len(bustard_bc) > 0:
        return os.path.join(bustard_bc[0], "fastq")
    # otherwise assume we are in the fastq directory
    # XXX What other cases can we end up with here?
    else:
        return fc_dir

class GalaxySqnLimsApi:
    """Manage talking with the Galaxy REST api for sequencing information.
    """
    def __init__(self, base_url, user, passwd):
        self._base_url = base_url
        # build cookies so we keep track of being logged in
        cj = http_cookiejar.LWPCookieJar()
        opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cj))
        urllib.request.install_opener(opener)
        login = dict(email=user, password=passwd, login_button='Login')
        req = urllib.request.Request("%s/user/login" % self._base_url,
                urllib.parse.urlencode(login))
        response = urllib.request.urlopen(req)

    def run_details(self, run):
        """Retrieve sequencing run details as a dictionary.
        """
        run_data = dict(run=run)
        req = urllib.request.Request("%s/nglims/api_run_details" % self._base_url,
                urllib.parse.urlencode(run_data))
        response = urllib.request.urlopen(req)
        info = json.loads(response.read())
        if "error" in info:
            raise ValueError("Problem retrieving info: %s" % info["error"])
        else:
            return info["details"]
