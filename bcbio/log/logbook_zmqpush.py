"""Enable multiple processes to send logs to a central server through ZeroMQ.

Thanks to Zachary Voase: https://github.com/zacharyvoase/logbook-zmqpush
Slightly modified to support Logbook 0.4.1.
"""
import errno
import json
import socket

import logbook.queues
from logbook.base import LogRecord
import zmq
from zmq.utils.garbage import gc

MAX_SOCKETS = 32000

def _increase_gc_sockets():
    """Increase default sockets for zmq gc. Avoids scaling issues.
    https://github.com/zeromq/pyzmq/issues/471
    """
    ctx = zmq.Context()
    ctx.max_sockets = MAX_SOCKETS
    gc.context = ctx

class ZeroMQPushHandler(logbook.queues.ZeroMQHandler):

    """
    A handler that pushes JSON log records over a ZMQ socket.

    Specifically, this handler opens a ``zmq.PUSH`` socket and connects to a
    ``zmq.PULL`` socket at the specified address. You can use
    :class:`ZeroMQPullSubscriber` to receive the log record.

    Example:

        >>> import logbook
        >>> handler = ZeroMQPushHandler('tcp://127.0.0.1:5501')
        >>> with handler.applicationbound():
        ...     logbook.debug("Something happened")

    Switch off hostname injection with `hostname=False`:

        >>> handler = ZeroMQPushHandler('tcp://127.0.0.1:5501', hostname=False)
        >>> with handler.applicationbound():
        ...     logbook.debug("No hostname info")
    """

    def __init__(self, addr=None, level=logbook.NOTSET, filter=None,
                 bubble=False, context=None, hostname=True):
        logbook.Handler.__init__(self, level, filter, bubble)

        self.hostname = hostname
        if context is None:
            context = zmq.Context()
            self._context = context
        else:
            self._context = None
        self._context.max_sockets = MAX_SOCKETS
        _increase_gc_sockets()
        self.socket = context.socket(zmq.PUSH)
        if addr is not None:
            self.socket.connect(addr)

    def emit(self, record):
        if self.hostname:
            inject_hostname.process(record)
        return super(ZeroMQPushHandler, self).emit(record)

    def close(self):
        if self._context:
            self._context.destroy(linger=0)
        self.socket.close(linger=0)

class ZeroMQPullSubscriber(logbook.queues.ZeroMQSubscriber):

    """
    A subscriber which listens on a PULL socket for log records.

    This subscriber opens a ``zmq.PULL`` socket and binds to the specified
    address. You should probably use this in conjunction with
    :class:`ZeroMQPushHandler`.

    Example:

        >>> subscriber = ZeroMQPullSubscriber('tcp://*:5501')
        >>> log_record = subscriber.recv()
    """

    def __init__(self, addr=None, context=None):
        self._zmq = zmq
        self.context = context or zmq.Context.instance()
        self.context.max_sockets = MAX_SOCKETS
        _increase_gc_sockets()
        self.socket = self.context.socket(zmq.PULL)
        if addr is not None:
            self.socket.bind(addr)

    def recv(self, timeout=None):
        """Overwrite standard recv for timeout calls to catch interrupt errors.
        """
        if timeout:
            try:
                testsock = self._zmq.select([self.socket], [], [], timeout)[0]
            except zmq.ZMQError as e:
                if e.errno == errno.EINTR:
                    testsock = None
                else:
                    raise
            if not testsock:
                return
            rv = self.socket.recv(self._zmq.NOBLOCK)
            return LogRecord.from_dict(json.loads(rv))
        else:
            return super(ZeroMQPullSubscriber, self).recv(timeout)

@logbook.Processor
def inject_hostname(log_record):
    """A Logbook processor to inject the current hostname into log records."""
    log_record.extra['source'] = socket.gethostname()

def inject(**params):

    """
    A Logbook processor to inject arbitrary information into log records.

    Simply pass in keyword arguments and use as a context manager:

        >>> with inject(identifier=str(uuid.uuid4())).applicationbound():
        ...     logger.debug('Something happened')
    """

    def callback(log_record):
        log_record.extra.update(params)
    return logbook.Processor(callback)
