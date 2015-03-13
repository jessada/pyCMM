import logging
#from pycmm import settings
import sys

#def init(file_name=settings.DFLT_LOG_FILE,
#         level=DEBUG,
#         ):
#    root = logging.getLogger('abcd')
#    root.setLevel(logging.DEBUG)
#
#    ch = logging.StreamHandler(sys.stdout)
#    ch.setLevel(logging.DEBUG)
#    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
#    ch.setFormatter(formatter)
#    root.addHandler(ch)
#    #logging.basicConfig(filename=file_name,
#    #                    level=level,
#    #                    )

mylog_file = None
logging.basicConfig(level=logging.DEBUG,
                    format='## [%(levelname)s]\t%(asctime)s - %(name)-12s - %(message)s',
                    )
lg = logging.getLogger('root')

#INFO_FMT  = "## [INFO]    {msg}"
#WARN_FMT  = "## [WARNING] {msg}"
#DEBUG_FMT = "## [DEBUG]   {msg}"
#ERROR_FMT = "## [ERROR]   {msg}"

def getLogger(name):
    global lg
    lg = logging.getLogger(name)

#def init(log_file=None,
#         dbg_mode=False,
#         ):
#        mylog_file = log_file
#        debug_mode=dbg_mode

#def basicConfig(level=logging.DEBUG):
#    logging.basicConfig(level=level)
    #logging.basicConfig(filename=filename,level=level)

## **************  defining basic functions  **************
#def write_log(msg):
##    if 'mylog_file' in globals():
#    if mylog_file is not None:
#        f = open(mylog_file, "a+")
#        print >> f, msg
#        f.close()

#def output_msg(msg):
#    print >> sys.stderr, msg
#    write_log(msg)

def info(msg):
    lg.info(msg)
#    formated_msg=INFO_FMT.format(msg=msg)
#    output_msg(formated_msg)

def warning(msg):
    lg.warning(msg)
#def warn(msg):
#    formated_msg=WARN_FMT.format(msg=msg)
#    output_msg(formated_msg)

def debug(msg):
    lg.debug(msg)
#    formated_msg=DEBUG_FMT.format(msg=msg)
#    if debug_mode:
#        output_msg(formated_msg)
#    else:
#        write_log(formated_msg)

def throw(err_msg):
    formated_msg=ERROR_FMT.format(msg=err_msg)
    write_log(formated_msg)
    raise Exception(formated_msg)
