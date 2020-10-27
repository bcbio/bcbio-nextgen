
class DummyCM(object):
    """Explicit structure of a context manager.
    Allows to monkey-patch context managers defined through
    @contextmanager decorator.

    value: holds whatever the context manager should yield.
    """
    value = None

    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        pass

    def __getattr__(self, attr):
        return self.value.__getattribute__(attr)


class BaseTx(DummyCM, str):
    def __iadd__(self, other):
        """support for += operator to make os.path.join work
        with the value yielded by the dummy context manager"""
        return self.value + other

    def __str__(self):
        return self.value


class DummyTxTmpdir(BaseTx):
    value = 'foo'  # path that dummy tx_tmpdir "yields"


class DummyFileTransaction(BaseTx):
    value = 'tx_file'
