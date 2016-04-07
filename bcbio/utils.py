"""Helpful utilities for building analysis pipelines.
"""
import gzip
import os
import tempfile
import time
import shutil
import contextlib
import itertools
import functools
import random
import ConfigParser
import collections
import fnmatch
import subprocess
import sys
import types

import toolz as tz
import yaml
try:
    from concurrent import futures
except ImportError:
    try:
        import futures
    except ImportError:
        futures = None

@contextlib.contextmanager
def cpmap(cores=1):
    """Configurable parallel map context manager.

    Returns appropriate map compatible function based on configuration:
    - Local single core (the default)
    - Multiple local cores
    """
    if int(cores) == 1:
        yield itertools.imap
    else:
        if futures is None:
            raise ImportError("concurrent.futures not available")
        pool = futures.ProcessPoolExecutor(cores)
        yield pool.map
        pool.shutdown()

def map_wrap(f):
    """Wrap standard function to easily pass into 'map' processing.
    """
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return apply(f, *args, **kwargs)
    return wrapper

def transform_to(ext):
    """
    Decorator to create an output filename from an output filename with
    the specified extension. Changes the extension, in_file is transformed
    to a new type.

    Takes functions like this to decorate:
    f(in_file, out_dir=None, out_file=None) or,
    f(in_file=in_file, out_dir=None, out_file=None)

    examples:
    @transform(".bam")
    f("the/input/path/file.sam") ->
        f("the/input/path/file.sam", out_file="the/input/path/file.bam")

    @transform(".bam")
    f("the/input/path/file.sam", out_dir="results") ->
        f("the/input/path/file.sam", out_file="results/file.bam")

    """

    def decor(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            out_file = kwargs.get("out_file", None)
            if not out_file:
                in_path = kwargs.get("in_file", args[0])
                out_dir = kwargs.get("out_dir", os.path.dirname(in_path))
                safe_makedir(out_dir)
                out_name = replace_suffix(os.path.basename(in_path), ext)
                out_file = os.path.join(out_dir, out_name)
            kwargs["out_file"] = out_file
            if not file_exists(out_file):
                out_file = f(*args, **kwargs)
            return out_file
        return wrapper
    return decor


def filter_to(word):
    """
    Decorator to create an output filename from an input filename by
    adding a word onto the stem. in_file is filtered by the function
    and the results are written to out_file. You would want to use
    this over transform_to if you don't know the extension of the file
    going in. This also memoizes the output file.

    Takes functions like this to decorate:
    f(in_file, out_dir=None, out_file=None) or,
    f(in_file=in_file, out_dir=None, out_file=None)

    examples:
    @filter_to(".foo")
    f("the/input/path/file.sam") ->
        f("the/input/path/file.sam", out_file="the/input/path/file.foo.bam")

    @filter_to(".foo")
    f("the/input/path/file.sam", out_dir="results") ->
        f("the/input/path/file.sam", out_file="results/file.foo.bam")

    """

    def decor(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            out_file = kwargs.get("out_file", None)
            if not out_file:
                in_path = kwargs.get("in_file", args[0])
                out_dir = kwargs.get("out_dir", os.path.dirname(in_path))
                safe_makedir(out_dir)
                out_name = append_stem(os.path.basename(in_path), word)
                out_file = os.path.join(out_dir, out_name)
            kwargs["out_file"] = out_file
            if not file_exists(out_file):
                out_file = f(*args, **kwargs)
            return out_file
        return wrapper
    return decor


def memoize_outfile(ext=None, stem=None):
    """
    Memoization decorator.

    See docstring for transform_to and filter_to for details.
    """
    if ext:
        return transform_to(ext)
    if stem:
        return filter_to(stem)

def to_single_data(input):
    """Convert an input to a single bcbio data/world object.

    Handles both single sample cases (CWL) and all sample cases (standard bcbio).
    """
    if (isinstance(input, (list, tuple)) and len(input) == 1):
        return input[0]
    else:
        assert isinstance(input, dict), input
        return input

def unpack_worlds(items):
    """Handle all the ways we can pass multiple samples for back-compatibility.
    """
    # Unpack nested lists of samples grouped together (old IPython style)
    if isinstance(items[0], (list, tuple)) and len(items[0]) == 1:
        out = []
        for d in items:
            assert len(d) == 1 and isinstance(d[0], dict), len(d)
            out.append(d[0])
    # Unpack a single argument with multiple samples (CWL style)
    elif isinstance(items, (list, tuple)) and len(items) == 1 and isinstance(items[0], (list, tuple)):
        out = items[0]
    else:
        out = items
    return out

def safe_makedir(dname):
    """Make a directory if it doesn't exist, handling concurrent race conditions.
    """
    if not dname:
        return dname
    num_tries = 0
    max_tries = 5
    while not os.path.exists(dname):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dname)
        except OSError:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dname

@contextlib.contextmanager
def chdir(new_dir):
    """Context manager to temporarily change to a new directory.

    http://lucentbeing.com/blog/context-managers-and-the-with-statement-in-python/
    """
    cur_dir = os.getcwd()
    safe_makedir(new_dir)
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(cur_dir)

@contextlib.contextmanager
def tmpfile(*args, **kwargs):
    """Make a tempfile, safely cleaning up file descriptors on completion.
    """
    (fd, fname) = tempfile.mkstemp(*args, **kwargs)
    try:
        yield fname
    finally:
        os.close(fd)
        if os.path.exists(fname):
            os.remove(fname)

def file_exists(fname):
    """Check if a file exists and is non-empty.
    """
    try:
        return fname and os.path.exists(fname) and os.path.getsize(fname) > 0
    except OSError:
        return False

def file_uptodate(fname, cmp_fname):
    """Check if a file exists, is non-empty and is more recent than cmp_fname.
    """
    try:
        return (file_exists(fname) and file_exists(cmp_fname) and
                os.path.getmtime(fname) >= os.path.getmtime(cmp_fname))
    except OSError:
        return False

def create_dirs(config, names=None):
    if names is None:
        names = config["dir"].keys()
    for dname in names:
        d = config["dir"][dname]
        safe_makedir(d)

def save_diskspace(fname, reason, config):
    """Overwrite a file in place with a short message to save disk.

    This keeps files as a sanity check on processes working, but saves
    disk by replacing them with a short message.
    """
    if config["algorithm"].get("save_diskspace", False):
        with open(fname, "w") as out_handle:
            out_handle.write("File removed to save disk space: %s" % reason)

def read_galaxy_amqp_config(galaxy_config, base_dir):
    """Read connection information on the RabbitMQ server from Galaxy config.
    """
    galaxy_config = add_full_path(galaxy_config, base_dir)
    config = ConfigParser.ConfigParser()
    config.read(galaxy_config)
    amqp_config = {}
    for option in config.options("galaxy_amqp"):
        amqp_config[option] = config.get("galaxy_amqp", option)
    return amqp_config

def add_full_path(dirname, basedir=None):
    if basedir is None:
        basedir = os.getcwd()
    if not dirname.startswith("/"):
        dirname = os.path.join(basedir, dirname)
    return dirname

def splitext_plus(f):
    """Split on file extensions, allowing for zipped extensions.
    """
    base, ext = os.path.splitext(f)
    if ext in [".gz", ".bz2", ".zip"]:
        base, ext2 = os.path.splitext(base)
        ext = ext2 + ext
    return base, ext

def remove_safe(f):
    try:
        if os.path.isdir(f):
            shutil.rmtree(f)
        else:
            os.remove(f)
    except OSError:
        pass

def move_safe(origin, target):
    """
    Move file, skip if exists
    """
    if origin == target:
        return origin
    if file_exists(target):
        return target
    shutil.move(origin, target)
    return target

def file_plus_index(fname):
    """Convert a file name into the file plus required indexes.
    """
    exts = {".vcf": ".idx", ".bam": ".bai", ".vcf.gz": ".tbi", ".bed.gz": ".tbi",
            ".fq.gz": ".gbi"}
    ext = splitext_plus(fname)[-1]
    if ext in exts:
        return [fname, fname + exts[ext]]
    else:
        return [fname]

def copy_plus(orig, new):
    """Copy a fils, including biological index files.
    """
    for ext in ["", ".idx", ".gbi", ".tbi", ".bai"]:
        if os.path.exists(orig + ext) and (not os.path.lexists(new + ext) or not os.path.exists(new + ext)):
            shutil.copyfile(orig + ext, new + ext)

def symlink_plus(orig, new):
    """Create relative symlinks and handle associated biological index files.
    """
    for ext in ["", ".idx", ".gbi", ".tbi", ".bai"]:
        if os.path.exists(orig + ext) and (not os.path.lexists(new + ext) or not os.path.exists(new + ext)):
            with chdir(os.path.dirname(new)):
                remove_safe(new + ext)
               # Work around symlink issues on some filesystems. Randomly
               # fail to symlink.
                try:
                    os.symlink(os.path.relpath(orig + ext), os.path.basename(new + ext))
                except OSError:
                    if not os.path.exists(new + ext) or not os.path.lexists(new + ext):
                        remove_safe(new + ext)
                        shutil.copyfile(orig + ext, new + ext)
    orig_noext = splitext_plus(orig)[0]
    new_noext = splitext_plus(new)[0]
    for sub_ext in [".bai"]:
        if os.path.exists(orig_noext + sub_ext) and not os.path.lexists(new_noext + sub_ext):
            with chdir(os.path.dirname(new_noext)):
                os.symlink(os.path.relpath(orig_noext + sub_ext), os.path.basename(new_noext + sub_ext))

def open_gzipsafe(f):
    return gzip.open(f) if f.endswith(".gz") else open(f)

def append_stem(to_transform, word):
    """
    renames a filename or list of filenames with 'word' appended to the stem
    of each one:
    example: append_stem("/path/to/test.sam", "_filtered") ->
    "/path/to/test_filtered.sam"

    """
    if is_sequence(to_transform):
        return [append_stem(f, word) for f in to_transform]
    elif is_string(to_transform):
        (base, ext) = splitext_plus(to_transform)
        return "".join([base, word, ext])
    else:
        raise ValueError("append_stem takes a single filename as a string or "
                         "a list of filenames to transform.")

def replace_suffix(to_transform, suffix):
    """
    replaces the suffix on a filename or list of filenames
    example: replace_suffix("/path/to/test.sam", ".bam") ->
    "/path/to/test.bam"

    """
    if is_sequence(to_transform):
        transformed = []
        for f in to_transform:
            (base, _) = os.path.splitext(f)
            transformed.append(base + suffix)
        return transformed
    elif is_string(to_transform):
        (base, _) = os.path.splitext(to_transform)
        return base + suffix
    else:
        raise ValueError("replace_suffix takes a single filename as a string or "
                         "a list of filenames to transform.")

# ## Functional programming

def partition_all(n, iterable):
    """Partition a list into equally sized pieces, including last smaller parts
    http://stackoverflow.com/questions/5129102/python-equivalent-to-clojures-partition-all
    """
    it = iter(iterable)
    while True:
        chunk = list(itertools.islice(it, n))
        if not chunk:
            break
        yield chunk

def robust_partition_all(n, iterable):
    """
    replaces partition_all with a more robust version.
    Workaround for a segfault in pybedtools when using a BedTool as an iterator:
    https://github.com/daler/pybedtools/issues/88 for the discussion
    """
    it = iter(iterable)
    while True:
        x = []
        for _ in range(n):
            try:
                x.append(it.next())
            except StopIteration:
                yield x
                # Omitting this StopIteration results in a segfault!
                raise StopIteration
        yield x

def partition(pred, iterable, tolist=False):
    'Use a predicate to partition entries into false entries and true entries'
    # partition(is_odd, range(10)) --> 0 2 4 6 8   and  1 3 5 7 9
    t1, t2 = itertools.tee(iterable)
    ifalse = itertools.ifilterfalse(pred, t1)
    itrue = itertools.ifilter(pred, t2)
    if tolist:
        return list(ifalse), list(itrue)
    else:
        return ifalse, itrue

# ## Dealing with configuration files

def merge_config_files(fnames):
    """Merge configuration files, preferring definitions in latter files.
    """
    def _load_yaml(fname):
        with open(fname) as in_handle:
            config = yaml.load(in_handle)
        return config
    out = _load_yaml(fnames[0])
    for fname in fnames[1:]:
        cur = _load_yaml(fname)
        for k, v in cur.iteritems():
            if k in out and isinstance(out[k], dict):
                out[k].update(v)
            else:
                out[k] = v
    return out

def deepish_copy(org):
    """Improved speed deep copy for dictionaries of simple python types.

    Thanks to Gregg Lind:
    http://writeonly.wordpress.com/2009/05/07/deepcopy-is-a-pig-for-simple-data/
    """
    out = dict().fromkeys(org)
    for k, v in org.iteritems():
        if isinstance(v, dict):
            out[k] = deepish_copy(v)
        else:
            try:
                out[k] = v.copy()   # dicts, sets
            except AttributeError:
                try:
                    out[k] = v[:]   # lists, tuples, strings, unicode
                except TypeError:
                    out[k] = v      # ints
    return out

def get_in(d, t, default=None):
    """
    look up if you can get a tuple of values from a nested dictionary,
    each item in the tuple a deeper layer

    example: get_in({1: {2: 3}}, (1, 2)) -> 3
    example: get_in({1: {2: 3}}, (2, 3)) -> {}
    """
    return tz.get_in(t, d, default)

def flatten(l):
    """
    flatten an irregular list of lists
    example: flatten([[[1, 2, 3], [4, 5]], 6]) -> [1, 2, 3, 4, 5, 6]
    lifted from: http://stackoverflow.com/questions/2158395/

    """
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el,
                                                                   basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el


def is_sequence(arg):
    """
    check if 'arg' is a sequence

    example: arg([]) -> True
    example: arg("lol") -> False

    """
    return (not hasattr(arg, "strip") and
            hasattr(arg, "__getitem__") or
            hasattr(arg, "__iter__"))


def is_pair(arg):
    """
    check if 'arg' is a two-item sequence

    """
    return is_sequence(arg) and len(arg) == 2

def is_string(arg):
    return isinstance(arg, basestring)


def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.'''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


def itersubclasses(cls, _seen=None):
    """
    snagged from:  http://code.activestate.com/recipes/576949/
    itersubclasses(cls)

    Generator over all subclasses of a given class, in depth first order.

    >>> list(itersubclasses(int)) == [bool]
    True
    >>> class A(object): pass
    >>> class B(A): pass
    >>> class C(A): pass
    >>> class D(B,C): pass
    >>> class E(D): pass
    >>>
    >>> for cls in itersubclasses(A):
    ...     print(cls.__name__)
    B
    D
    E
    C
    >>> # get ALL (new-style) classes currently defined
    >>> [cls.__name__ for cls in itersubclasses(object)] #doctest: +ELLIPSIS
    ['type', ...'tuple', ...]
    """

    if not isinstance(cls, type):
        raise TypeError('itersubclasses must be called with '
                        'new-style classes, not %.100r' % cls)
    if _seen is None:
        _seen = set()
    try:
        subs = cls.__subclasses__()
    except TypeError:  # fails only when cls is type
        subs = cls.__subclasses__(cls)
    for sub in subs:
        if sub not in _seen:
            _seen.add(sub)
            yield sub
            for sub in itersubclasses(sub, _seen):
                yield sub

def replace_directory(out_files, dest_dir):
    """
    change the output directory to dest_dir
    can take a string (single file) or a list of files

    """
    if is_sequence(out_files):
        filenames = map(os.path.basename, out_files)
        return [os.path.join(dest_dir, x) for x in filenames]
    elif is_string(out_files):
        return os.path.join(dest_dir, os.path.basename(out_files))
    else:
        raise ValueError("in_files must either be a sequence of filenames "
                         "or a string")

def which(program):
    """ returns the path to an executable or None if it can't be found"""

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def reservoir_sample(stream, num_items, item_parser=lambda x: x):
    """
    samples num_items from the stream keeping each with equal probability
    """
    kept = []
    for index, item in enumerate(stream):
        if index < num_items:
            kept.append(item_parser(item))
        else:
            r = random.randint(0, index)
            if r < num_items:
                kept[r] = item_parser(item)
    return kept


def compose(f, g):
    return lambda x: f(g(x))

def dictapply(d, fn):
    """
    apply a function to all non-dict values in a dictionary
    """
    for k, v in d.items():
        if isinstance(v, dict):
            v = dictapply(v, fn)
        else:
            d[k] = fn(v)
    return d

def Rscript_cmd():
    """Retrieve path to locally installed Rscript or first in PATH.

    Prefers Rscript version installed via conda to a system version.
    """
    rscript = which(os.path.join(os.path.dirname(os.path.realpath(sys.executable)), "Rscript"))
    if rscript:
        return rscript
    else:
        return which("Rscript")

def R_sitelib():
    """Retrieve the R site-library installed with the bcbio installer.
    """
    from bcbio import install
    return os.path.join(install.get_defaults().get("tooldir", "/usr/local"),
                        "lib", "R", "site-library")

def R_package_path(package):
    """
    return the path to an installed R package
    """
    local_sitelib = R_sitelib()
    rscript = Rscript_cmd()
    cmd = """{rscript} -e '.libPaths(c("{local_sitelib}")); find.package("{package}")'"""
    try:
        output = subprocess.check_output(cmd.format(**locals()), shell=True)
    except subprocess.CalledProcessError, e:
        return None
    for line in output.split("\n"):
        if "[1]" not in line:
            continue
        dirname = line.split("[1]")[1].replace("\"", "").strip()
        if os.path.exists(dirname):
            return dirname
    return None

def perl_cmd():
    """Retrieve path to locally installed conda Perl or first in PATH.
    """
    perl = which(os.path.join(os.path.dirname(os.path.realpath(sys.executable)), "perl"))
    if perl:
        return perl
    else:
        return which("perl")

def get_perl_exports(tmpdir=None):
    """Environmental exports to use conda installed perl.
    """
    perl_path = os.path.dirname(perl_cmd())
    out = "unset PERL5LIB && export PATH=%s:$PATH" % (perl_path)
    if tmpdir:
        out += " && export TMPDIR=%s" % (tmpdir)
    return out

def is_gzipped(fname):
    _, ext = os.path.splitext(fname)
    return ext in [".gz", "gzip"]

def is_bzipped(fname):
    _, ext = os.path.splitext(fname)
    return ext in [".bz2", "bzip2"]

def open_possible_gzip(fname, flag="r"):
    if is_gzipped(fname):
        if "b" not in flag:
            flag += "b"
        return gzip.open(fname, flag)
    else:
        return open(fname, flag)

def filter_missing(xs):
    """
    remove items from a list if they evaluate to False
    """
    return filter(lambda x: x, xs)

def rbind(dfs):
    """
    acts like rbind for pandas dataframes
    """
    if len(dfs) == 1:
        return dfs[0]
    df = dfs[0]
    for d in dfs[1:]:
        df = df.append(d)
    return df

def max_command_length():
    """
    get the maximum length of the command line, in bytes, defaulting
    to a conservative number if not set
    http://www.in-ulm.de/~mascheck/various/argmax/
    """
    DEFAULT_MAX_LENGTH = 150000 # lowest seen so far is 200k
    try:
        arg_max = os.sysconf('SC_ARG_MAX')
        env_lines = len(os.environ) * 4
        env_chars = sum([len(x) + len(y) for x, y in os.environ.iteritems()])
        arg_length = arg_max - env_lines - 2048
    except ValueError:
        arg_length = DEFAULT_MAX_LENGTH
    return arg_length if arg_length > 0 else DEFAULT_MAX_LENGTH

# LazyImport from NIPY
# https://github.com/nipy/nitime/blob/master/nitime/lazyimports.py

class LazyImport(types.ModuleType):
    """
    This class takes the module name as a parameter, and acts as a proxy for
    that module, importing it only when the module is used, but effectively
    acting as the module in every other way (including inside IPython with
    respect to introspection and tab completion) with the *exception* of
    reload()- reloading a :class:`LazyImport` raises an :class:`ImportError`.
    >>> mlab = LazyImport('matplotlib.mlab')
    No import happens on the above line, until we do something like call an
    ``mlab`` method or try to do tab completion or introspection on ``mlab``
    in IPython.
    >>> mlab
    <module 'matplotlib.mlab' will be lazily loaded>
    Now the :class:`LazyImport` will do an actual import, and call the dist
    function of the imported module.
    >>> mlab.dist(1969,2011)
    42.0
    """
    def __getattribute__(self, x):
        # This method will be called only once, since we'll change
        # self.__class__ to LoadedLazyImport, and __getattribute__ will point
        # to module.__getattribute__
        name = object.__getattribute__(self, '__name__')
        __import__(name)
        # if name above is 'package.foo.bar', package is returned, the docs
        # recommend that in order to get back the full thing, that we import
        # and then lookup the full name is sys.modules, see:
        # http://docs.python.org/library/functions.html#__import__
        module = sys.modules[name]
        # Now that we've done the import, cutout the middleman and make self
        # act as the imported module
        class LoadedLazyImport(types.ModuleType):
            __getattribute__ = module.__getattribute__
            __repr__ = module.__repr__
        object.__setattr__(self, '__class__', LoadedLazyImport)
        # The next line will make "reload(l)" a silent no-op
        # sys.modules[name] = self
        return module.__getattribute__(x)
    def __repr__(self):
        return "<module '%s' will be lazily loaded>" %\
                object.__getattribute__(self,'__name__')
