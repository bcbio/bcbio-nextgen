from bcbio.chipseq import atac

def run(_, data, out_dir):
    """Standard QC metrics for ATAC-seq
    """
    return atac.run_ataqv(data)
