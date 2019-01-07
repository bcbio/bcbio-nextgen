"""Update the PATH environment variable to reflect the value of BCBIOPATH.

"""
from __future__ import print_function
import contextlib
import os

from bcbio import utils

def _prepend(original, to_prepend):
    """Prepend paths in a string representing a list of paths to another.

    original and to_prepend are expected to be strings representing
    os.pathsep-separated lists of filepaths.

    If to_prepend is None, original is returned.

    The list of paths represented in the returned value consists of the first of
    occurrences of each non-empty path in the list obtained by prepending the
    paths in to_prepend to the paths in original.

    examples:
    # Unix
    _prepend('/b:/d:/a:/d', '/a:/b:/c:/a') -> '/a:/b:/c:/d'
    _prepend('/a:/b:/a', '/a:/c:/c')       -> '/a:/c:/b'
    _prepend('/c', '/a::/b:/a')            -> '/a:/b:/c'
    _prepend('/a:/b:/a', None)             -> '/a:/b:/a'
    _prepend('/a:/b:/a', '')               -> '/a:/b'

    """

    if to_prepend is None:
        return original

    sep = os.pathsep

    def split_path_value(path_value):
        return [] if path_value == '' else path_value.split(sep)

    seen = set()

    components = []
    for path in split_path_value(to_prepend) + split_path_value(original):
        if path not in seen and path != '':
            components.append(path)
            seen.add(path)

    return sep.join(components)


def prepend_bcbiopath():
    """Prepend paths in the BCBIOPATH environment variable (if any) to PATH.

    Uses either a pre-sent global environmental variable (BCBIOPATH) or the
    local anaconda directory.
    """
    if os.environ.get('BCBIOPATH'):
        os.environ['PATH'] = _prepend(os.environ.get('PATH', ''),
                                      os.environ.get('BCBIOPATH', None))
    else:
        os.environ['PATH'] = _prepend(os.environ.get('PATH', ''),
                                      utils.get_bcbio_bin())

def remove_bcbiopath():
    """Remove bcbio internal path from first element in PATH.

    Useful when we need to access remote programs, like Java 7 for older
    installations.
    """
    to_remove = os.environ.get("BCBIOPATH", utils.get_bcbio_bin()) + ":"
    if os.environ["PATH"].startswith(to_remove):
        os.environ["PATH"] = os.environ["PATH"][len(to_remove):]

@contextlib.contextmanager
def orig_paths():
    """Run using original paths in a `with` block.
    """
    remove_bcbiopath()
    yield None
    prepend_bcbiopath()

if __name__ == '__main__':

    # check examples (poor-man's doctest)
    print(_prepend('/b:/d:/a:/d', '/a:/b:/c:/a'))
    print(_prepend('/a:/b:/a', '/a:/c:/c'))
    print(_prepend('/c', '/a::/b:/a'))
    print(_prepend('/a:/b:/a', None))
    print(_prepend('/a:/b:/a', ''))
