"""Run Broad's RNA-SeqQC tool and handle reporting of useful summary metrics.
"""
# soft imports
try:
    import pandas as pd
    import statsmodels.formula.api as sm
except ImportError:
    pd, sm = None, None

from bcbio import bam
import bcbio.pipeline.datadict as dd

def starts_by_depth(bam_file, data, sample_size=10000000):
    """
    Return a set of x, y points where x is the number of reads sequenced and
    y is the number of unique start sites identified
    If sample size < total reads in a file the file will be downsampled.
    """
    binsize = (sample_size / 100) + 1
    seen_starts = set()
    counted = 0
    num_reads = []
    starts = []
    buffer = []
    downsampled = bam.downsample(bam_file, data, sample_size)
    with bam.open_samfile(downsampled) as samfile:
        for read in samfile:
            if read.is_unmapped:
                continue
            counted += 1
            buffer.append(str(read.tid) + ":" + str(read.pos))
            if counted % binsize == 0:
                seen_starts.update(buffer)
                buffer = []
                num_reads.append(counted)
                starts.append(len(seen_starts))
        seen_starts.update(buffer)
        num_reads.append(counted)
        starts.append(len(seen_starts))
    return pd.DataFrame({"reads": num_reads, "starts": starts})


def estimate_library_complexity(df, algorithm="RNA-seq"):
    """
    estimate library complexity from the number of reads vs.
    number of unique start sites. returns "NA" if there are
    not enough data points to fit the line
    """
    DEFAULT_CUTOFFS = {"RNA-seq": (0.25, 0.40)}
    cutoffs = DEFAULT_CUTOFFS[algorithm]
    if len(df) < 5:
        return {"unique_starts_per_read": 'nan',
                "complexity": "NA"}
    model = sm.ols(formula="starts ~ reads", data=df)
    fitted = model.fit()
    slope = fitted.params["reads"]
    if slope <= cutoffs[0]:
        complexity = "LOW"
    elif slope <= cutoffs[1]:
        complexity = "MEDIUM"
    else:
        complexity = "HIGH"

    # for now don't return the complexity flag
    return {"Unique Starts Per Read": float(slope)}
    # return {"unique_start_per_read": float(slope),
    #         "complexity": complexity}
