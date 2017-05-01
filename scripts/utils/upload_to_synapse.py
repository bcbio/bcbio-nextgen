#!/usr/bin/env python
"""Upload bcbio reference materials and inputs to a Synapse project.

Uses the API to recursively upload a directory of files.

Usage:
   upload_to_synapse.py wf_name to_upload synapse_parent_id

"""
import os
import sys

import synapseclient
import synapseutils

def main(wf_name, to_upload, synapse_parent_id):
    syn = synapseclient.Synapse()
    syn.login(os.environ["SYNAPSE_LOGIN"], apiKey=os.environ["SYNAPSE_API_KEY"])

    s_base_folder, remotes = _accumulate_remotes(synapse_parent_id, syn)

    for dirpath, _, filenames in os.walk(to_upload):
        remote_dirpath = os.path.join(s_base_folder.name, dirpath)
        if filenames:
            s_folder, remotes = _remote_folder(remote_dirpath, remotes, syn)
            for filename in filenames:
                remote_filename = os.path.join(remote_dirpath, filename)
                if remote_filename not in remotes:
                    filename = os.path.join(dirpath, filename)
                    if os.path.getsize(filename) > 0:
                        print("Uploading %s" % filename)
                        f = synapseclient.File(filename, parent=s_folder)
                        f.workflow = wf_name
                        f.workflowOption = "all"
                        s_filename = syn.store(f)
                        remotes[remote_filename] = s_filename.id

def _accumulate_remotes(synapse_parent_id, syn):
    """Retrieve references to all remote directories and files.
    """
    remotes = {}
    s_base_folder = syn.get(synapse_parent_id)
    for (s_dirpath, s_dirpath_id), _, s_filenames in synapseutils.walk(syn, synapse_parent_id):
        remotes[s_dirpath] = s_dirpath_id
        if s_filenames:
            for s_filename, s_filename_id in s_filenames:
                remotes[os.path.join(s_dirpath, s_filename)] = s_filename_id
    return s_base_folder, remotes

def _remote_folder(dirpath, remotes, syn):
    """Retrieve the remote folder for files, creating if necessary.
    """
    if dirpath in remotes:
        return remotes[dirpath], remotes
    else:
        parent_dir, cur_dir = os.path.split(dirpath)
        parent_folder, remotes = _remote_folder(parent_dir, remotes, syn)
        s_cur_dir = syn.store(synapseclient.Folder(cur_dir, parent=parent_folder))
        remotes[dirpath] = s_cur_dir.id
        return s_cur_dir.id, remotes

if __name__ == "__main__":
    main(*sys.argv[1:])
