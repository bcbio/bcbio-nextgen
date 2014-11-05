#!/usr/bin/env python
"""Convert time stamps from bcbio logs into hourly timings per step.

Usage:
  /path/to/bcbio/anaconda/python summarize_timing.py your_work_dir/log/bcbio-nextgen.log

This parses out the "Timing" steps in the log file to provide timing
per analysis block. If your analysis had a lot of stops and restarts,
you'll need to manually reconstruct a reasonable estimate of timings
around the stops/failures. The best approach is to get a file of Timinig
only data with `grep Timing log/bcbio-nextgen.log` then manually add in
failure/stop points with timestamps and "unexpected error" for the
description.
"""
import sys

import arrow

def main(log_file):
    prev_time = None
    total_time = None
    cur_process = None
    with open(log_file) as in_handle:
        for line in in_handle:
            if line.strip() and line.find("Timing:") > 0 or line.lower().find("unexpected error") > 0:
                cur_time = arrow.get(line.split("]")[0])
                if prev_time:
                    tdiff = cur_time - prev_time
                    if cur_process.lower().find("unexpected error") == -1:
                        if total_time:
                            total_time += tdiff
                        else:
                            total_time = tdiff
                    if cur_process.lower().find("unexpected error") == -1:
                        print cur_process, "  ", tdiff
                cur_process = line.split(":")[-1].strip()
                prev_time = cur_time
    print "Total  ", total_time

if __name__ == "__main__":
    main(sys.argv[1])
