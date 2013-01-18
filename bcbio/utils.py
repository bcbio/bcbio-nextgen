"""Helpful utilities for building analysis pipelines.
"""
import os
import tempfile
import time
import shutil
import contextlib
import itertools
import functools
import ConfigParser
try:
    import multiprocessing
    from multiprocessing.pool import IMapIterator
except ImportError:
    multiprocessing = None
import collections
import yaml

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
        if multiprocessing is None:
            raise ImportError("multiprocessing not available")
        # Fix to allow keyboard interrupts in multiprocessing: https://gist.github.com/626518
        def wrapper(func):
            def wrap(self, timeout=None):
                return func(self, timeout=timeout if timeout is not None else 1e100)
            return wrap
        IMapIterator.next = wrapper(IMapIterator.next)
        # recycle threads on Python 2.7; remain compatible with Python 2.6
        try:
            pool = multiprocessing.Pool(int(cores), maxtasksperchild=1)
        except TypeError:
            pool = multiprocessing.Pool(int(cores))
        yield pool.imap_unordered
        pool.terminate()

def map_wrap(f):
    """Wrap standard function to easily pass into 'map' processing.
    """
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return apply(f, *args, **kwargs)
    return wrapper

def memoize_outfile(ext):
    """Creates outfile from input file and ext, running if outfile not present.

    This requires a standard function usage. The first arg, or kwarg 'in_file', needs
    to be the input file that is being processed. The output name is created with the
    provided ext relative to the input. The function is only run if the created
    out_file is not present.
    """
    def decor(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            if len(args) > 0:
                in_file = args[0]
            else:
                in_file = kwargs['in_file']
            out_file = "%s%s" % (os.path.splitext(in_file)[0], ext)
            if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
                kwargs['out_file'] = out_file
                f(*args, **kwargs)
            return out_file
        return wrapper
    return decor

def safe_makedir(dname):
    """Make a directory if it doesn't exist, handling concurrent race conditions.
    """
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
def curdir_tmpdir(remove=True, base_dir=None):
    """Context manager to create and remove a temporary directory.

    This can also handle a configured temporary directory to use.
    """
    if base_dir is not None:
        tmp_dir_base = os.path.join(base_dir, "bcbiotmp")
    else:
        tmp_dir_base = os.path.join(os.getcwd(), "tmp")
    safe_makedir(tmp_dir_base)
    tmp_dir = tempfile.mkdtemp(dir=tmp_dir_base)
    safe_makedir(tmp_dir)
    try :
        yield tmp_dir
    finally :
        if remove:
            try:
                shutil.rmtree(tmp_dir)
            except:
                pass

@contextlib.contextmanager
def chdir(new_dir):
    """Context manager to temporarily change to a new directory.

    http://lucentbeing.com/blog/context-managers-and-the-with-statement-in-python/
    """
    cur_dir = os.getcwd()
    safe_makedir(new_dir)
    os.chdir(new_dir)
    try :
        yield
    finally :
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
    return os.path.exists(fname) and os.path.getsize(fname) > 0

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


def append_stem(filename, word, delim="_"):
    """
    returns a filename with 'word' appended to the stem
    example: append_stem("/path/to/test.sam", "filtered") ->
    "/path/to/test_filtered.sam"

    """
    (base, ext) = os.path.splitext(filename)
    return "".join([base, delim, word, ext])


def replace_suffix(filename, suffix):
    """
    replace the suffix of filename with suffix
    example: replace_suffix("/path/to/test.sam", ".bam") ->
    "/path/to/test.bam"

    """
    (base, _) = os.path.splitext(filename)
    return base + suffix

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
            if out.has_key(k) and isinstance(out[k], dict):
                out[k].update(v)
            else:
                out[k] = v
    return out


def get_in(d, t, default=None):
    """
    look up if you can get a tuple of values from a nested dictionary,
    each item in the tuple a deeper layer

    example: get_in({1: {2: 3}}, (1, 2)) -> 3
    example: get_in({1: {2: 3}}, (2, 3)) -> {}
    """
    result = reduce(lambda d, t: d.get(t, {}), t, d)
    if not result:
        return default
    else:
        return result


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
