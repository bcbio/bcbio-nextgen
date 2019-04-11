"""Summarize amplification and loss of heterozygosity (LOH) from heterogeneity callers.

Provides high level summaries of calls in regions of interest.
"""
import csv
import collections
import os
import decimal
import uuid

import six
from six import StringIO
import toolz as tz
import yaml

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd

# Standard sets of coordinates we always include
_COORDS = {"LOH":
           {"hg38": {"HLA": ("chr6", 28510120, 33480577),
                     "B2M": ("chr15", 44711487, 44718877)},
            "hg19": {"HLA": ("chr6", 29640000, 33120000),
                     "B2M": ("chr15", 45003675, 45011075)},
            "GRCh37": {"HLA": ("6", 29640000, 33120000),
                       "B2M": ("15", 45003675, 45011075)}}}

def get_coords(data):
    """Retrieve coordinates of genes of interest for prioritization.

    Can read from CIViC input data or a supplied BED file of chrom, start, end
    and gene information.
    """
    for category, vtypes in [("LOH", {"LOSS", "HETEROZYGOSITY"}),
                             ("amplification", {"AMPLIFICATION"})]:
        out = tz.get_in([category, dd.get_genome_build(data)], _COORDS, {})
        priority_file = dd.get_svprioritize(data)
        if priority_file:
            if os.path.basename(priority_file).find("civic") >= 0:
                for chrom, start, end, gene in _civic_regions(priority_file, vtypes, dd.get_disease(data)):
                    out[gene] = (chrom, start, end)
            elif os.path.basename(priority_file).find(".bed") >= 0:
                for line in utils.open_gzipsafe(priority_file):
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        chrom, start, end, gene = parts[:4]
                        out[gene] = (chrom, int(start), int(end))
        yield category, out

def _matches(tocheck, target):
    for t in target:
        t = t.lower()
        for c in tocheck:
            if c.lower().find(t) >= 0:
                return True

def _civic_regions(civic_file, variant_types=None, diseases=None, drugs=None):
    """Retrieve gene regions and names filtered by variant_types and diseases.
    """
    if isinstance(diseases, six.string_types):
        diseases = [diseases]
    with utils.open_gzipsafe(civic_file) as in_handle:
        reader = csv.reader(in_handle, delimiter="\t")
        for chrom, start, end, info_str in reader:
            info = edn_loads(info_str)
            if not variant_types or _matches(info["support"]["variants"], variant_types):
                if not diseases or _matches(info["support"]["diseases"], diseases):
                    if not drugs or _matches(info["support"]["drugs"], drugs):
                        yield (chrom, int(start), int(end), list(info["name"])[0])

def summary_status(call, data):
    """Retrieve status in regions of interest, along with heterogeneity metrics.

    Provides output with overall purity and ploidy, along with region
    specific calls.
    """
    out_file = None
    if call.get("vrn_file") and os.path.exists(call.get("vrn_file")):
        out_file = os.path.join(os.path.dirname(call["vrn_file"]),
                                "%s-%s-lohsummary.yaml" % (dd.get_sample_name(data), call["variantcaller"]))
        if not utils.file_uptodate(out_file, call["vrn_file"]):
            out = {}
            if call["variantcaller"] == "titancna":
                out.update(_titancna_summary(call, data))
                pass
            elif call["variantcaller"] == "purecn":
                out.update(_purecn_summary(call, data))
            if out:
                out["description"] = dd.get_sample_name(data)
                out["variantcaller"] = call["variantcaller"]
                with file_transaction(data, out_file) as tx_out_file:
                    with open(tx_out_file, "w") as out_handle:
                        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file if out_file and os.path.exists(out_file) else None

def _check_copy_number_changes(svtype, cn, minor_cn, data):
    """Check if copy number changes match the expected svtype.
    """
    if svtype == "LOH" and minor_cn == 0:
        return svtype
    elif svtype == "amplification" and cn > dd.get_ploidy(data):
        return svtype
    else:
        return "std"

def _to_cn(v):
    return int(round(float(v)))

def _titancna_summary(call, data):
    """Summarize purity, ploidy and LOH for TitanCNA.
    """
    out = {}
    for svtype, coords in get_coords(data):
        cur_calls = {k: collections.defaultdict(int) for k in coords.keys()}
        with open(call["subclones"]) as in_handle:
            header = in_handle.readline().strip().split()
            for line in in_handle:
                val = dict(zip(header, line.strip().split()))
                start = int(val["Start_Position.bp."])
                end = int(val["End_Position.bp."])
                for region, cur_coords in coords.items():
                    if val["Chromosome"] == cur_coords[0] and are_overlapping((start, end), cur_coords[1:]):
                        cur_calls[region][_check_copy_number_changes(svtype, _to_cn(val["Copy_Number"]),
                                                                     _to_cn(val["MinorCN"]), data)] += 1
        out[svtype] = {r: _merge_cn_calls(c, svtype) for r, c in cur_calls.items()}

    with open(call["hetsummary"]) as in_handle:
        vals = dict(zip(in_handle.readline().strip().split("\t"), in_handle.readline().strip().split("\t")))
    out["purity"] = vals["purity"]
    out["ploidy"] = vals["ploidy"]
    return out

def _purecn_summary(call, data):
    """Summarize purity, ploidy and LOH for PureCN.
    """
    out = {}
    for svtype, coords in get_coords(data):
        cur_calls = {k: collections.defaultdict(int) for k in coords.keys()}
        with open(call["loh"]) as in_handle:
            in_handle.readline()  # header
            for line in in_handle:
                _, chrom, start, end, _, cn, minor_cn = line.split(",")[:7]
                start = int(start)
                end = int(end)
                for region, cur_coords in coords.items():
                    if chrom == cur_coords[0] and are_overlapping((start, end), cur_coords[1:]):
                        cur_calls[region][_check_copy_number_changes(svtype, _to_cn(cn), _to_cn(minor_cn), data)] += 1
        out[svtype] = {r: _merge_cn_calls(c, svtype) for r, c in cur_calls.items()}
    with open(call["hetsummary"]) as in_handle:
        vals = dict(zip(in_handle.readline().strip().replace('"', '').split(","),
                        in_handle.readline().strip().split(",")))
    out["purity"] = vals["Purity"]
    out["ploidy"] = vals["Ploidy"]
    return out

def _merge_cn_calls(calls, svtype):
    if calls[svtype]:
        return "mixed" if calls["std"] else svtype
    else:
        return "no"

def are_overlapping(r, s):
    """Test if two coordinates overlap.
    https://stackoverflow.com/a/27182551
    """
    return r[1] >= s[0] and s[1] >= r[0]

# ## EDN parser
# Thanks to https://github.com/sunng87/pyclj
# Slightly adapter to avoid external dependencies

def edn_load(fp):
    decoder = CljDecoder(fp)
    return decoder.decode()

def edn_loads(s):
    buf = StringIO(s)
    result = edn_load(buf)
    buf.close()
    return result

def _number(v):
    if v.endswith('M'):
        out = decimal.Decimal(v[:-1])
    else:
        try:
            out = int(v)
        except ValueError as e:
            out = float(v)
    return out


_STOP_CHARS = [" ", ",", "\n", "\r", "\t"]
_COLL_OPEN_CHARS = ["#", "[", "{", "("]
_COLL_CLOSE_CHARS = ["]", "}", ")"]
_EXTRA_NUM_CHARS = ["-", "+", ".", "e", "E", "M"]

class CljDecoder(object):
    def __init__(self, fd):
        self.fd = fd
        self.cur_line = 1
        self.cur_pos = 1
        self.value_stack = []
        self.terminator = None ## for collection type

    def decode(self):
        while True:
            v = self.__read_token()
            if len(self.value_stack) == 0:
                return v

    def __seek_back(self, size):
        self.fd.seek(self.fd.tell()-size, 0)

    def __read_and_back(self, size):
        s = self.fd.read(size)
        self.__seek_back(size)
        return s

    def __get_type_from_char(self, c):
        """return a tuple of type information
        * type name
        * a flag to indicate if it's a collection
        """
        if c.isdigit() or c =='-':
            return ("number", False, None)
        elif c == 't' or c == 'f': ## true/false
            return ("boolean", False, None)
        elif c == 'n': ## nil
            return ("nil", False, None)
        elif c == '\\' :
            return ("char", False, None)
        elif c == ':':
            return ("keyword", False, None)
        elif c == '"':
            return ("string", False, None)
        elif c == '#':
            if self.__read_and_back(1) == '{':
                return ("set", True, "}")
            if self.__read_and_back(1) == ':':
                return ("namespaced_dict", True, "}")
            if self.__read_and_back(4) == 'inst':
                return ("datetime", False, None)
            if self.__read_and_back(4) == 'uuid':
                return ("uuid", False, None)
        elif c == '{':
            return ("dict", True, "}")
        elif c == '(':
            return ("list", True, ")")
        elif c == '[':
            return ('list', True, "]")

        return (None, False, None)

    def __read_fd(self, size):
        if size == 1:
            c = self.fd.read(size)
            if  c == '\n':
                self.cur_pos = 0
                self.cur_line = self.cur_line + 1
            return c
        else:
            self.cur_pos = self.cur_pos + size
            cs = self.fd.read(size)
            return cs

    def __read_token(self):
        c = self.__read_fd(1)

        ## skip all stop chars if necessary
        while c in _STOP_CHARS:
            c = self.__read_fd(1)

        ## raise exception when unexpected EOF found
        if c == '':
            raise ValueError("Unexpected EOF")

        t, coll, term = self.__get_type_from_char(c)
        if coll:
            ## move cursor
            if t == "set":
                ## skip {
                self.__read_fd(1)
            namespace = None
            if t == "namespaced_dict":
                ## skip :
                self.__read_fd(1)
                ## get namespace
                buf = []
                while c != '{':
                    c = self.__read_fd(1)
                    buf.append(c)
                namespace = ''.join(buf[:-1])

            self.terminator = term

            self.value_stack.append(([], self.terminator, t, namespace))
            return None
        else:
            v = None ## token value
            e = None ## end char
            r = True ## the token contains data or not

            if t == "boolean":
                if c == 't':
                    chars = self.__read_fd(4)
                    if chars[:3] != 'rue':
                        raise ValueError('Expect true, got t%s at line %d, col %d' % (chars[:3], self.cur_line, self.cur_pos))
                    e = chars[-1]
                    v = True
                else:
                    chars = self.__read_fd(5)
                    if chars[:4] != 'alse':
                        raise ValueError('Expect true, got t%s at line %d, col %d' % (chars[:3], self.cur_line, self.cur_pos))
                    e = chars[-1]
                    v = False

            elif t == "char":
                buf = []
                while c is not self.terminator and c is not "" and c not in _STOP_CHARS:
                    c = self.__read_fd(1)
                    buf.append(c)

                e = c
                v = ''.join(buf[:-1])

            elif t == "nil":
                chars = self.__read_fd(3)
                if chars[:2] != 'il':
                    raise ValueError('Expect nil, got n%s at line %d, col %d' % (chars[:2], self.cur_line, self.cur_pos))
                e = chars[-1]
                v = None

            elif t == "number":
                buf = []
                while c.isdigit() or (c in _EXTRA_NUM_CHARS):
                    buf.append(c)
                    c = self.__read_fd(1)
                e = c
                numstr = ''.join(buf)
                v = _number(numstr)

                ## special case for
                ## [23[12]]
                ## this is a valid clojure form
                if e in _COLL_OPEN_CHARS:
                    self.__seek_back(1)

            elif t == "keyword":
                buf = []    ##skip the leading ":"
                while c is not self.terminator and c is not "" and c not in _STOP_CHARS:
                    c = self.__read_fd(1)
                    buf.append(c)

                e = c
                v = ''.join(buf[:-1])

            elif t == "string":
                buf = []
                cp = c = self.__read_fd(1) ## to check escaping character \

                while not(c == '"' and cp != '\\'):
                    buf.append(c)
                    cp = c
                    c = self.__read_fd(1)
                e = c
                v = unicode(''.join(buf).decode('unicode-escape'))

            elif t == "datetime":
                ## skip "inst"
                self.__read_fd(4)

                ## read next value as string
                s = self.__read_token()
                if not isinstance(s, six.string_types):
                    raise ValueError('Str expected, but got %s' % str(s))

                ## remove read string from the value_stack
                if len(self.value_stack) > 0:
                    self.value_stack[-1][0].pop()
                e = '"'
                v = pyrfc3339.parse(s)

            elif t == "uuid":
                ## skip "uuid"
                self.__read_fd(4)

                ## read next value as string
                s = self.__read_token()
                if not isinstance(s, six.string_types):
                    raise ValueError('Str expected, but got %s' % str(s))

                ## remove read string from the value_stack
                if len(self.value_stack) > 0:
                    self.value_stack[-1][0].pop()
                e = '"'
                v = uuid.UUID(s)

            else:
                if c not in _COLL_CLOSE_CHARS:
                    raise ValueError('Unexpected char: "%s" at line %d, col %d' % (c, self.cur_line, self.cur_pos))
                r = False
                e = c

            if e == self.terminator:
                current_scope, _, container, namespace = self.value_stack.pop()

                if r:
                    current_scope.append(v)

                if container == "set":
                    try:
                        v = set(current_scope)
                    except TypeError:
                        v = tuple(current_scope)
                elif container == "list":
                    v = current_scope
                elif container in ["dict", "namespaced_dict"]:
                    v = {}
                    for i in range(0, len(current_scope), 2):
                        key = '%s/%s' % (namespace, current_scope[i]) if namespace else current_scope[i]
                        v[key] = current_scope[i+1]
                r = True

            if r and len(self.value_stack) > 0:
                self.value_stack[-1][0].append(v)
                self.terminator = self.value_stack[-1][1]

            return v
