"""Update the PATH environment variable to reflect the value of BCBIOPATH.

"""

import os

def _maybe_prepend(original, to_prepend):
    """Prepend unique non-empty paths unless already present.

    original and to_prepend are expected to be strings corresponding to
    os.pathsep-separated lists of filepaths.

    If to_prepend is the empty sting, original is returned.

    examples:
    # Unix
    _maybe_prepend('/foo:/bar:/foo', '/foo:/baz:/baz') -> '/baz:/foo:/bar:/foo'
    _maybe_prepend('/baz', '/foo::/bar:/foo')          -> '/foo:/bar:/baz'
    _maybe_prepend('/foo:/bar:/foo', '')               -> '/foo:/bar:/foo'

    """

    sep = os.pathsep

    def split_path_value(path_value):
        return [] if path_value == '' else path_value.split(sep)

    new_paths = split_path_value(to_prepend)

    if new_paths == []:
        return original

    paths = split_path_value(original)

    paths_lookup = set(paths)

    # collect in prefix all unique non-empty paths in new_paths that are not
    # already in paths
    prefix = []
    for path in new_paths:
        if path not in paths_lookup and path != '':
            prefix.append(path)
            paths_lookup.add(path)

    return sep.join(prefix + paths)


if __name__ == '__main__':

    # check examples (poor-man's doctest)
    print _maybe_prepend('/foo:/bar:/foo', '/foo:/baz:/baz')
    print _maybe_prepend('/baz', '/foo::/bar:/foo')
    print _maybe_prepend('/foo:/bar:/foo', '')

else:

    # possibly update PATH
    os.environ['PATH'] = _maybe_prepend(os.environ.get('PATH', ''),
                                        os.environ.get('BCBIOPATH', ''))
    del _maybe_prepend
