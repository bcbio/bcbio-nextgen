"""Useful utilities for handling CWL inputs and outputs.

This is shared functionality abstracted across multiple approaches, currently
mostly handling CWL records. This needs some generalization to apply across
non-variant calling workflows.
"""
import collections
import pprint

import toolz as tz

def _get_all_cwlkeys(items):
    """Retrieve cwlkeys from inputs, handling defaults which can be null.

    When inputs are null in some and present in others, this creates unequal
    keys in each sample, confusing decision making about which are primary and extras.
    """
    default_keys = set(["metadata__batch", "config__algorithm__validate",
                        "config__algorithm__validate_regions",
                        "config__algorithm__validate_regions_merged",
                        "validate__summary",
                        "validate__tp", "validate__fp", "validate__fn",
                        "config__algorithm__coverage", "config__algorithm__coverage_merged"])
    all_keys = set([])
    for data in items:
        all_keys.update(set(data["cwl_keys"]))
    all_keys.update(default_keys)
    return all_keys

def split_data_cwl_items(items):
    """Split a set of CWL output dictionaries into data samples and CWL items.

    Handles cases where we're arrayed on multiple things, like a set of regional
    VCF calls and data objects.
    """
    key_lens = set([])
    for data in items:
        key_lens.add(len(_get_all_cwlkeys([data])))
    extra_key_len = min(list(key_lens)) if len(key_lens) > 1 else None
    data_out = []
    extra_out = []
    for data in items:
        if extra_key_len and len(_get_all_cwlkeys([data])) == extra_key_len:
            extra_out.append(data)
        else:
            data_out.append(data)
    if len(extra_out) == 0:
        return data_out, {}
    else:
        cwl_keys = extra_out[0]["cwl_keys"]
        for extra in extra_out[1:]:
            cur_cwl_keys = extra["cwl_keys"]
            assert cur_cwl_keys == cwl_keys, pprint.pformat(extra_out)
        cwl_extras = collections.defaultdict(list)
        for data in items:
            for key in cwl_keys:
                cwl_extras[key].append(data[key])
        data_final = []
        for data in data_out:
            for key in cwl_keys:
                data.pop(key)
            data_final.append(data)
        return data_final, dict(cwl_extras)

def samples_to_records(samples):
    """Convert samples into output CWL records.
    """
    from bcbio.pipeline import run_info
    RECORD_CONVERT_TO_LIST = set(["config__algorithm__tools_on", "config__algorithm__tools_off",
                                  "reference__genome_context"])
    all_keys = _get_all_cwlkeys(samples)
    out = []
    for data in samples:
        for raw_key in sorted(list(all_keys)):
            key = raw_key.split("__")
            if tz.get_in(key, data) is None:
                data = tz.update_in(data, key, lambda x: None)
                data["cwl_keys"].append(raw_key)
            if raw_key in RECORD_CONVERT_TO_LIST:
                val = tz.get_in(key, data)
                if not val: val = []
                elif not isinstance(val, (list, tuple)): val = [val]
                data = tz.update_in(data, key, lambda x: val)
        data["metadata"] = run_info.add_metadata_defaults(data.get("metadata", {}))
        out.append(data)
    return out
