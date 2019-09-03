"""Helpful utilities for building analysis pipelines.
"""
import glob
import gzip
import os
import tempfile
import time
import shutil
import contextlib
import itertools
import functools
import random
import fnmatch
import subprocess
import sys
import types

import six
import toolz as tz
import yaml

from collections import Mapping, OrderedDict


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
        return f(*args, **kwargs)
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
    # On busy filesystems can have issues accessing main directory. Allow retries
    num_tries = 0
    max_tries = 5
    cur_dir = None
    while cur_dir is None:
        try:
            cur_dir = os.getcwd()
        except OSError:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
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


def get_size(path):
    """ Returns the size in bytes if `path` is a file,
        or the size of all files in `path` if it's a directory.
        Analogous to `du -s`.
    """
    if os.path.isfile(path):
        return os.path.getsize(path)
    return sum(get_size(os.path.join(path, f)) for f in os.listdir(path))


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
        for ext in ["", ".bai"]:
            if os.path.exists(fname + ext):
                with open(fname + ext, "w") as out_handle:
                    out_handle.write("File removed to save disk space: %s" % reason)

def read_galaxy_amqp_config(galaxy_config, base_dir):
    """Read connection information on the RabbitMQ server from Galaxy config.
    """
    galaxy_config = add_full_path(galaxy_config, base_dir)
    config = six.moves.configparser.ConfigParser()
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

def remove_plus(orig):
    """Remove a fils, including biological index files.
    """
    for ext in ["", ".idx", ".gbi", ".tbi", ".bai"]:
        if os.path.exists(orig + ext):
            remove_safe(orig + ext)

def copy_plus(orig, new):
    """Copy a fils, including biological index files.
    """
    for ext in ["", ".idx", ".gbi", ".tbi", ".bai"]:
        if os.path.exists(orig + ext) and (not os.path.lexists(new + ext) or not os.path.exists(new + ext)):
            shutil.copyfile(orig + ext, new + ext)

def symlink_plus(orig, new):
    """Create relative symlinks and handle associated biological index files.
    """
    orig = os.path.abspath(orig)
    if not os.path.exists(orig):
        raise RuntimeError("File not found: %s" % orig)
    for ext in ["", ".idx", ".gbi", ".tbi", ".bai", ".fai"]:
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
    for sub_ext in [".bai", ".dict"]:
        if os.path.exists(orig_noext + sub_ext) and not os.path.lexists(new_noext + sub_ext):
            with chdir(os.path.dirname(new_noext)):
                os.symlink(os.path.relpath(orig_noext + sub_ext), os.path.basename(new_noext + sub_ext))

def open_gzipsafe(f, is_gz=False):
    if f.endswith(".gz") or is_gz:
        if six.PY3:
            return gzip.open(f, "rt", encoding="utf-8", errors="ignore")
        else:
            return gzip.open(f)
    else:
        if six.PY3:
            return open(f, encoding="utf-8", errors="ignore")
        else:
            return open(f)

def is_empty_gzipsafe(f):
    h = open_gzipsafe(f)
    is_empty = len(h.read(1)) > 0
    h.close()
    return is_empty

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
                x.append(next(it))
            except StopIteration:
                yield x
                # Omitting this StopIteration results in a segfault!
                raise StopIteration
        yield x

def partition(pred, iterable, tolist=False):
    'Use a predicate to partition entries into false entries and true entries'
    # partition(is_odd, range(10)) --> 0 2 4 6 8   and  1 3 5 7 9
    t1, t2 = itertools.tee(iterable)
    ifalse = six.moves.filterfalse(pred, t1)
    itrue = six.moves.filter(pred, t2)
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
            config = yaml.safe_load(in_handle)
        return config
    out = _load_yaml(fnames[0])
    for fname in fnames[1:]:
        cur = _load_yaml(fname)
        for k, v in cur.items():
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
    for k, v in org.items():
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

def safe_to_float(x):
    """Convert to float, handling None and non-float inputs.

    Useful for cleaning complicated output from variant callers.
    """
    if x is None:
        return None
    else:
        try:
            return float(x)
        except ValueError:
            return None

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
        if isinstance(el, (list, tuple)):
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
    return (not is_string(arg) and
            (hasattr(arg, "__getitem__") or
             hasattr(arg, "__iter__")))


def is_pair(arg):
    """
    check if 'arg' is a two-item sequence

    """
    return is_sequence(arg) and len(arg) == 2

def is_string(arg):
    return isinstance(arg, six.string_types)


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


def which(program, env=None):
    """ returns the path to an executable or None if it can't be found"""
    if env is None:
        env = os.environ.copy()

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in env["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    for path in get_all_conda_bins():
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
    rscript = which(os.path.join(get_bcbio_bin(), "Rscript"))
    if rscript:
        return rscript
    else:
        return which("Rscript")

def R_sitelib():
    """Retrieve the R site-library installed with the bcbio installer.
    """
    return os.path.join(os.path.dirname(get_bcbio_bin()), "lib", "R", "library")

def R_package_path(package):
    """
    return the path to an installed R package
    """
    local_sitelib = R_sitelib()
    rscript = Rscript_cmd()
    cmd = """{rscript} --vanilla -e '.libPaths(c("{local_sitelib}")); find.package("{package}")'"""
    try:
        output = subprocess.check_output(cmd.format(**locals()), shell=True)
    except subprocess.CalledProcessError as e:
        return None
    for line in output.decode().split("\n"):
        if "[1]" not in line:
            continue
        dirname = line.split("[1]")[1].replace("\"", "").strip()
        if os.path.exists(dirname):
            return dirname
    return None

def R_package_resource(package, resource):
    """
    return a path to an R package resource, if it is available
    """
    package_path = R_package_path(package)
    if not package_path:
        return None
    package_resource = os.path.join(package_path, resource)
    if not file_exists(package_resource):
        return None
    else:
        return package_resource

def get_java_binpath(cmd=None):
    """Retrieve path for java to use, handling custom BCBIO_JAVA_HOME

    Defaults to the dirname of cmd, or local anaconda directory
    """
    if os.environ.get("BCBIO_JAVA_HOME"):
        test_cmd = os.path.join(os.environ["BCBIO_JAVA_HOME"], "bin", "java")
        if os.path.exists(test_cmd):
            cmd = test_cmd
    if not cmd:
        cmd = Rscript_cmd()
    return os.path.dirname(cmd)

def clear_java_home():
    """Clear JAVA_HOME environment or reset to BCBIO_JAVA_HOME.

    Avoids accidental java injection but respects custom BCBIO_JAVA_HOME
    command.
    """
    if os.environ.get("BCBIO_JAVA_HOME"):
        test_cmd = os.path.join(os.environ["BCBIO_JAVA_HOME"], "bin", "java")
        if os.path.exists(test_cmd):
            return "export JAVA_HOME=%s" % os.environ["BCBIO_JAVA_HOME"]
    return "unset JAVA_HOME"

def get_java_clprep(cmd=None):
    """Correctly prep command line for java commands, setting PATH and unsetting JAVA_HOME.
    """
    return "%s && export PATH=%s:\"$PATH\"" % (clear_java_home(), get_java_binpath(cmd))

def get_R_exports():
    return "unset R_HOME && unset R_LIBS && export PATH=%s:\"$PATH\"" % (os.path.dirname(Rscript_cmd()))

def perl_cmd():
    """Retrieve path to locally installed conda Perl or first in PATH.
    """
    perl = which(os.path.join(get_bcbio_bin(), "perl"))
    if perl:
        return perl
    else:
        return which("perl")

def get_perl_exports(tmpdir=None):
    """Environmental exports to use conda installed perl.
    """
    perl_path = os.path.dirname(perl_cmd())
    out = "unset PERL5LIB && export PATH=%s:\"$PATH\"" % (perl_path)
    if tmpdir:
        out += " && export TMPDIR=%s" % (tmpdir)
    return out


def get_bcbio_env():
    env = os.environ.copy()
    env["PATH"] = append_path(get_bcbio_bin(), env['PATH'])
    return env


def append_path(bin, path, at_start=True):
    if at_start:
        tmpl = "{bin}:{path}"
    else:
        tmpl = "{path}:{bin}"
    return tmpl.format(bin=bin, path=path)


def get_bcbio_bin():
    return os.path.dirname(os.path.realpath(sys.executable))

def get_all_conda_bins():
    """Retrieve all possible conda bin directories, including environments.
    """
    bcbio_bin = get_bcbio_bin()
    conda_dir = os.path.dirname(bcbio_bin)
    if os.path.join("anaconda", "envs") in conda_dir:
        conda_dir = os.path.join(conda_dir[:conda_dir.rfind(os.path.join("anaconda", "envs"))], "anaconda")
    return [bcbio_bin] + list(glob.glob(os.path.join(conda_dir, "envs", "*", "bin")))

def get_program_python(cmd):
    """Get the full path to a python version linked to the command.

    Allows finding python based programs in python 2 versus python 3
    environments.
    """
    full_cmd = os.path.realpath(which(cmd))
    cmd_python = os.path.join(os.path.dirname(full_cmd), "python")
    env_python = None
    if "envs" in cmd_python:
        parts = cmd_python.split(os.sep)
        env_python = os.path.join(os.sep.join(parts[:parts.index("envs") + 2]), "bin", "python")
    if os.path.exists(cmd_python):
        return cmd_python
    elif env_python and os.path.exists(env_python):
        return env_python
    else:
        return os.path.realpath(sys.executable)

def local_path_export(at_start=True, env_cmd=None):
    """Retrieve paths to local install, also including environment paths if env_cmd included.
    """
    paths = [get_bcbio_bin()]
    if env_cmd:
        env_path = os.path.dirname(get_program_python(env_cmd))
        if env_path not in paths:
            paths.insert(0, env_path)
    if at_start:
        return "export PATH=%s:\"$PATH\" && " % (":".join(paths))
    else:
        return "export PATH=\"$PATH\":%s && " % (":".join(paths))

def locale_export():
    """Exports for dealing with Click-based programs and ASCII/Unicode errors.

    RuntimeError: Click will abort further execution because Python 3 was
    configured to use ASCII as encoding for the environment.
    Consult https://click.palletsprojects.com/en/7.x/python3/ for mitigation steps.
    """
    locale_to_use = get_locale()
    return "export LC_ALL=%s && export LANG=%s && " % (locale_to_use, locale_to_use)

def get_locale():
    """
    Looks up available locales on the system to find an appropriate one to pick,
    defaulting to C.UTF-8 which is globally available on newer systems. Prefers
    C.UTF-8 and en_US encodings, if available
    """
    default_locale = "C.UTF-8"
    preferred_locales = {"c.utf-8", "c.utf8", "en_us.utf-8", "en_us.utf8"}
    locale_to_use = None
    try:
        locales = subprocess.check_output(["locale", "-a"]).decode(errors="ignore").split("\n")
    except subprocess.CalledProcessError:
        locales = []
    # check for preferred locale
    for locale in locales:
        if locale.lower() in preferred_locales:
            locale_to_use = locale
            break
    # if preferred locale not available take first UTF-8 locale
    if not locale_to_use:
        for locale in locales:
            if locale.lower().endswith(("utf-8", "utf8")):
                locale_to_use = locale
                break
    # if locale listing not available, try using the default locale
    if not locale_to_use:
        locale_to_use = default_locale
    return locale_to_use

def java_freetype_fix():
    """Provide workaround for issues FreeType library symbols.

    libfontconfig.so.1: undefined symbol: FT_Done_MM_Var

    Cheap workaround with LD_PRELOAD, I don't know a better one.
    """
    return "export LD_PRELOAD=%s/lib/libfreetype.so && " % os.path.dirname(get_bcbio_bin())

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
        env_chars = sum([len(x) + len(y) for x, y in os.environ.items()])
        arg_length = arg_max - env_lines - 2048
    except ValueError:
        arg_length = DEFAULT_MAX_LENGTH
    return arg_length if arg_length > 0 else DEFAULT_MAX_LENGTH


def get_abspath(path, pardir=None):
    if pardir is None:
        pardir = os.getcwd()
    path = os.path.expandvars(path)
    return os.path.normpath(os.path.join(pardir, path))

def sort_filenames(filenames):
    """
    sort a list of files by filename only, ignoring the directory names
    """
    basenames = [os.path.basename(x) for x in filenames]
    indexes = [i[0] for i in sorted(enumerate(basenames), key=lambda x:x[1])]
    return [filenames[x] for x in indexes]

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

def walk_json(d, func):
    """ Walk over a parsed JSON nested structure `d`, apply `func` to each leaf element and replace it with result
    """
    if isinstance(d, Mapping):
        return OrderedDict((k, walk_json(v, func)) for k, v in d.items())
    elif isinstance(d, list):
        return [walk_json(v, func) for v in d]
    else:
        return func(d)
