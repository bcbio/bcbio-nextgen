"""
handle data from TCGA; from a list of filenames create a set of properly
paired up analyses

handles pairing primary and metastasized tumors with blood or solid normals
with the following rules:
1) use the blood sample as the matched normal if it exists
2) if more than one normal exists, choose the highest priority (blood > tissue)
normal and run the rest of the normals as tumors in their own batch
3 if multiple tumor sample types exist, combine them all together and run
   against the priority normal normal

So for example if there is:

1) 1 tumor (01) and 1 normal (10) run 01/10 as a batch
2) 1 tumor (01) and multiple normals (10, 11) run 01/10 as a batch
   and run 11 (the tissue tumor) as tumor alone
3) multiple tumors (01, 06) with 1 normal (10), run
   01+06/10 as a batch
4) multiple tumors (01, 06) with multiple normals (10, 11) run
   01+06/10 as a batch and 11 alone as a batch

normals and tumors are prioritized by the PRIOTIZED_TUMOR_CODES
and PRIORITIZED_NORMAL_CODES dictionaries
"""
import re
import os
import collections
import toolz as tz
import itertools
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from argparse import ArgumentParser


TCGA_RE = ("\w+.TCGA-"
           "(?P<tissue_source>\w{2})-"
           "(?P<participant>\w{4})-"
           "(?P<sample_type>\w{2})(?P<vial>\w{1})-"
           "(?P<portion>\w{2})(?P<analyte>\w{1})-"
           "(?P<plate>\w{4})-"
           "(?P<center>\w{2})."
           "(?P<version>\d*)")

# ordered to prefer the primary tumor
PRIORITIZED_TUMOR_CODES = {"01": 1, "02": 2, "03": 3, "04": 4, "05": 5,
                           "06": 6, "07": 7, "08": 8, "09": 9}
# ordered to prefer normals from blood
PRIORITIZED_NORMAL_CODES = {"10": 1, "11": 2, "12": 3, "13": 4, "14": 5}
VALID_TCGA_FIELDS = ["analyte", "batch", "center", "participant", "plate",
                     "portion", "sample_type", "tissue_source", "vial", "version"]
TCGA_HEADER = (["samplename", "description", "batch", "phenotype"] +
               VALID_TCGA_FIELDS)

def get_phenotype(sample):
    if "phenotype" in sample:
        return sample["phenotype"]
    sample_type = sample["sample_type"]
    if sample_type in PRIORITIZED_TUMOR_CODES.keys():
        return "tumor"
    elif sample_type in PRIORITIZED_NORMAL_CODES.keys():
        return "normal"
    else:
        return "unknown"

def is_tcga(fn):
    """
    detect if a filename could be from a TCGA dataset
    """
    return re.match(TCGA_RE, fn) is not None

def sample_to_bcbio(sample):
    samplename = os.path.splitext(sample["fn"])[0]
    phenotype = get_phenotype(sample)
    description = sample["batch"] + "-" + phenotype
    batch = sample["batch"]
    required =  ('{samplename},{description},{batch},'
                 '{phenotype}'.format(**locals()))
    return required + "," + ",".join([sample[x] for x in VALID_TCGA_FIELDS])

def prioritize_normals(metadata):
    normals = sorted(filter(lambda x: x["sample_type"] in
                            PRIORITIZED_NORMAL_CODES.keys(), metadata),
                     key=lambda x: PRIORITIZED_NORMAL_CODES[x["sample_type"]])
    if len(normals) == 0:
        return [], []
    normal_rest = [] if len(normals) < 2 else normals[1:]
    return normals[0], normal_rest

def rebatch_metadata_by_experiment(metadata):
    normal, normal_rest = prioritize_normals(metadata)
    batch = metadata[0]["participant"]
    tumor_batch = [tz.assoc(x, "batch", batch) for x in metadata
                   if x["sample_type"] in PRIORITIZED_TUMOR_CODES.keys()]
    normal = [tz.assoc(normal, "batch", batch)] if normal else []
    # run each non priority normal as its own tumor sample with no control
    normal_rest = [tz.assoc(x, "batch", batch + "-" + x["sample_type"]) for x
                   in normal_rest]
    normal_rest = [tz.assoc(x, "phenotype", "tumor") for x in normal_rest]
    all_batches = normal + normal_rest + tumor_batch
    return all_batches

def batch_tcga_metadata_by_participant(fns):
    metadata = [re.match(TCGA_RE, fn).groupdict() for fn in fns]
    metadata = [tz.assoc(d, "fn", fn) for d, fn in zip(metadata, fns)]
    participant = tz.groupby(lambda x: x["participant"], metadata).values()
    return participant

def keep_latest_sample(batch):
    groups = tz.groupby(lambda x: (x["participant"], x["sample_type"]),
                        batch).values()
    keep = [sorted(group, key=lambda x: int(x["version"]), reverse=True)
            for group in groups]
    return [x[0] for x in keep]

def tcga_to_bcbio_csv(fns, out_file):
    if file_exists(out_file):
        return out_file
    fns = set([x for x in fns if is_tcga(x)])
    batches = batch_tcga_metadata_by_participant(fns)
    batches = [keep_latest_sample(batch) for batch in batches]
    batches = [rebatch_metadata_by_experiment(x) for x in batches]
    batches = itertools.chain.from_iterable(batches)
    batches = set([sample_to_bcbio(x) for x in batches])
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            out_handle.write(",".join(TCGA_HEADER) + "\n")
            for batch in batches:
                out_handle.write(batch + "\n")
    return out_file

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("tcga_manifest", help="File of TCGA files to parse.")
    parser.add_argument("out_file", help="Name of bcbio_nextgen CSV file to write")
    args = parser.parse_args()

    with open(args.tcga_manifest) as in_handle:
        samples = [x.strip() for x in in_handle]
    tcga_to_bcbio_csv(samples, args.out_file)
