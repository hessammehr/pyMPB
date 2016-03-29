import logging
import sys
from defaults import log_format, log_datefmt

logbuffer = []
logger = None
errlogger = None
original_stderr = sys.stderr

def setup_logger(name, filename, quiet=False, redirect_stderr=False):
    """Set up and return a logger with a StreamHandler that prints to stdout 
    and a FileHandler that prints the same to filename. If a logger with the 
    same name exists already, it will be reused (i.e. logs will be printed to 
    the same file). If quiet, output to stdout will be disabled.
    If redirect_stderr, all errors going to sys.stderr from within Python
    will also be written to filename.

    """
    frm = logging.Formatter(fmt=log_format, datefmt=log_datefmt)

    global logger
    logger = logging.getLogger(name)

    if not quiet:
        stdouthandler = logging.StreamHandler(sys.stdout)
        stdouthandler.setFormatter(frm)
        # don't print DEBUG stuff to screen:
        stdouthandler.setLevel(logging.INFO)
        logger.addHandler(stdouthandler)

    if filename:
        fhandler = logging.FileHandler(filename, "a")
        fhandler.setFormatter(frm)
        fhandler.setLevel(logging.DEBUG)
        logger.addHandler(fhandler)

    logger.setLevel(logging.DEBUG)

    if redirect_stderr and filename:
        # don't bother with stderr redirection if we don't use fhandler

        global errlogger
        stderrhandler = logging.StreamHandler(original_stderr)
        stderrhandler.setFormatter(logging.Formatter(fmt='%(message)s'))
        stderrhandler.setLevel(logging.ERROR)
        #errlogger = logging.getLogger(name + '.STDERR')
        errlogger = logging.getLogger('STDERR')
        # remove any previously added handlers:
        for h in errlogger.handlers:
            errlogger.removeHandler(h)
        # errors should still go to original system stderr:
        errlogger.addHandler(stderrhandler)
        # but also to file:
        errlogger.addHandler(fhandler)
        # from now on, all errors raised from within Python will also be logged
        # to fhandler:
        sys.stderr = StreamToLogger(errlogger, logging.ERROR)

    push_buffered()

    return logger

def reset_logger():
    global logger, logbuffer, errlogger
    logbuffer = []
    # flush to stdout and to file and close file:
    if logger:
        for h in logger.handlers:
            h.close()
            logger.removeHandler(h)
    if errlogger:
        for h in errlogger.handlers:
            h.close()
            errlogger.removeHandler(h)
        sys.stderr = original_stderr
    logger = None
    errlogger = None

def check_initialized(level=logging.NOTSET, msg='', ):
    if logger is None:
        #print("WARNING: message submitted to logger, "
        #      "but logger is not initialized yet. Call log.setup_logger!")

        # logger is not initialized yet. 
        # Make a simple buffer with messages, which can be posted later:
        logbuffer.append((level, msg))
        return False
    else:
        return True

def push_buffered():
    global logbuffer
    if logbuffer:
        for level, msg in logbuffer:
            log(level, msg)
    logbuffer = []

def log(level, msg, *args, **kwargs):
    """Send a log message to the logger, specify the level."""
    if check_initialized(level, msg):
        logger.log(level, msg, *args, **kwargs)

def debug(msg, *args, **kwargs):
    """Send a debug message to the logger."""
    log(logging.DEBUG, msg, *args, **kwargs)

def info(msg, *args, **kwargs):
    """Send a info message to the logger."""
    log(logging.INFO, msg, *args, **kwargs)

def warning(msg, *args, **kwargs):
    """Send a warning message to the logger."""
    log(logging.WARNING, msg, *args, **kwargs)

def error(msg, *args, **kwargs):
    """Send an error message to the logger."""
    log(logging.ERROR, msg, *args, **kwargs)

def exception(msg, *args, **kwargs):
    """Send an exception message to the logger."""
    kwargs["exc_info"] = True
    log(logging.ERROR, msg, *args, **kwargs)

def critical(msg, *args, **kwargs):
    """Send a critical message to the logger."""
    log(logging.CRITICAL, msg, *args, **kwargs)


class StreamToLogger(object):
    """File-like stream object that redirects writes to a logger instance.
    Writes are line buffered, meaning the writes go the the logger only when
    a newline is written; everything not terminated with a newline goes into
    a buffer. This buffer is written the next time a newline is written, or
    when flush is called.

    (Buffering stderr seems not a good design choice, since it might be
    possible that the interpreter stops due to an error and quits before the
    buffer could be emptied. But unfortunately, Python errors are send word-by-
    word to the stream, so would result in a lot of individual calls to the
    logger otherwise. On the other hand, line-buffering of stderr is standard
    in Python3, and is not really a problem, on the grounds that errors outputs
    are always terminated with a newline, so the buffer will always be emptied.)

    """
    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        parts = buf.rpartition('\n')
        if parts[1] == '\n':
            # flush everything upto last newline, buffer everything afterwards:
            self.logger.log(self.log_level, self.linebuf + parts[0])
            self.linebuf = parts[2]
        else:
            self.linebuf += parts[2]

    def flush(self):
        if self.linebuf:
            self.logger.log(self.log_level, self.linebuf)
        self.linebuf = ''
        for h in self.logger.handlers:
            h.flush()