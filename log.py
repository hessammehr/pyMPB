import logging
import sys
from defaults import log_format, log_datefmt

logbuffer = []
logger = None

def setup_logger(name, filename, quiet=False):
    """Set up and return a logger with a StreamHandler that prints to stdout 
    and a FileHandler that prints the same to filename. If a logger with the 
    same name exists already, it will be reused (i.e. logs will be printed to 
    the same file). If quiet, output to stdout will be disabled.
    
    """
    frm = logging.Formatter(fmt=log_format, datefmt=log_datefmt)
    if not quiet:
        handler1 = logging.StreamHandler(sys.stdout)
        handler1.setFormatter(frm)
        # don't print DEBUG stuff to screen:
        handler1.setLevel(logging.INFO)
    if filename:
        handler2 = logging.FileHandler(filename, "w")
        handler2.setFormatter(frm)
        handler2.setLevel(logging.DEBUG)

    global logger
    logger = logging.getLogger(name)
    if not quiet:
        logger.addHandler(handler1)
    if filename:
        logger.addHandler(handler2)
    logger.setLevel(logging.DEBUG)
    
    push_buffered()
    
    return logger

def reset_logger():
    global logger, logbuffer
    logbuffer = []
    # flush to stdout and to file and close file:
    for h in logger.handlers:
        h.close()
    logger = None

def check_initialized(func=None, text='', ):
    if logger is None:
        #print("WARNING: message submitted to logger, "
        #      "but logger is not initialized yet. Call log.setup_logger!")
        
        # logger is not initialized yet. 
        # Make a simple buffer with messages, which can be posted later:
        logbuffer.append((func, text))
        return False
    else:
        return True
    
def push_buffered():
    global logbuffer
    if logbuffer:
        for func, text in logbuffer:
            func(text)
    logbuffer = []

def debug(text):
    if check_initialized(debug, text):
        logger.debug(text)

def info(text):
    if check_initialized(info, text):
        logger.info(text)
    
def warning(text):
    if check_initialized(warning, text):
        logger.warning(text)  
        
def error(text):
    if check_initialized(error, text):
        logger.error(text) 
       
def critical(text):
    if check_initialized(critical, text):
        logger.critical(text)        
    
    
def log_info_parameters(text, parameters):
    #TODO is this used?
    lstr = text + '\n    parameters: '
    for name, val in parameters.__dict__.iteritems():
        if not name.startswith('__'):
            lstr += name + "=" + repr(val) + ", "
    logger.info(lstr)