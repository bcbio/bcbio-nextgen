from bcbio.chipseq import atac

def run(_, data, out_dir):
    """Standard QC metrics for ATAC-seq
    """
    out = {}
    atac.run_ataqv(data)
    return out
