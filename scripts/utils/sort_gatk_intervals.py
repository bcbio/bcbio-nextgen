#!/usr/bin/env python
"""Sort GATK interval lists based on a sequence dictionary.

Usage:
    sort_gatk_intervals.py <interval file> [<sequence dictionary>]

The sequence dictionary is needed if the original interval file does not contain
it. The output file is a Picard style file with a sequence dictionary and
semi-SAM formatted lines.
"""
import sys
import os

def main(interval_file, seqdict_file=None):
    out_file = "%s-sort.interval_list" % (
            os.path.splitext(interval_file)[0].replace(".", "-"))
    with open(out_file, "w") as out_handle:
        with open(interval_file) as in_handle:
            if seqdict_file is None:
                chr_indexes, seqdict = read_dict(in_handle)
            else:
                with open(seqdict_file) as seqdict_handle:
                    chr_indexes, seqdict = read_dict(seqdict_handle)
            out_handle.write(seqdict)
        all_parts = []
        with open(interval_file) as in_handle:
            for parts in read_intervals(in_handle):
                try:
                    all_parts.append(((chr_indexes[parts[0]], int(parts[1]),
                        int(parts[2])), parts))
                except KeyError:
                    print parts[0]
        all_parts.sort()
        for (_, parts) in all_parts:
            out_handle.write("\t".join(parts) + "\n")

def read_intervals(in_handle):
    for i, line in enumerate(l for l in in_handle if not l.startswith("@")):
        parts = line.rstrip("\r\n").split()
        if len(parts) == 1:
            chr_name, loc = parts[0].split(":")
            start, end = loc.split("-")
            yield (chr_name, start, end, "+", "interval_%s" % i)
        elif len(parts) == 6:
            chr_name, start, end, strand, name, _ = parts
            #chr_name, start, end, name, _, strand = parts
            yield (chr_name, start, end, strand, name)
        elif len(parts) == 5:
            yield parts
        else:
            raise NotImplementedError(parts)

def read_dict(in_handle):
    parts = []
    chr_indexes = dict()
    cur_index = 0
    while 1:
        line = in_handle.readline()
        if not line.startswith("@"):
            break
        parts.append(line)
        if line.startswith("@SQ"):
            sn_part = [p for p in line.split("\t") if p.startswith("SN:")][0]
            (_, chr_name) = sn_part.split(":")
            chr_indexes[chr_name] = cur_index
            cur_index += 1
    return chr_indexes, "".join(parts)

if __name__ == "__main__":
    main(*sys.argv[1:])
