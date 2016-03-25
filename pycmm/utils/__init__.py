import subprocess
import datetime
from os.path import dirname
from os.path import realpath
from pycmm.utils import mylogger

def get_path(file_name):
    return dirname(realpath(file_name))

def get_file_prefix(file_name, file_ext):
    if file_name.endswith(file_ext):
        return file_name[:-len(file_ext)]
    return file_name

def set_log_file(raw_file):
    if raw_file is not None:
        time_stamp = datetime.datetime.now()
        log_file = log_file_with_time_stamp(raw_file,
                                            time_stamp,
                                            )
        mylogger.set_log_file(log_file)
        return log_file
    else:
        return None

def log_file_with_time_stamp(raw_file, time_stamp):
    log_file = get_file_prefix(raw_file, '.log')
    log_file += '_'
    log_file += time_stamp.strftime("%Y%m%d%H%M%S")
    log_file += '.log'
    return log_file

def exec_sh(cmd):
    mylogger.debug("executing: " + repr(cmd))
    p = subprocess.Popen(cmd,
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         )
    stdout_data, stderr_data = p.communicate()
    return_code = p.returncode
    print stdout_data
    if return_code:
        mylogger.throw("Error found during execute command '%s' with error code: %d, %s" % (cmd, return_code, stderr_data))
    elif stderr_data:
        print stderr_data
    return p, stdout_data

def concat_files(in_files,
                 out_file):
    cmd = "cat"
    if isinstance(in_files, str):
        cmd += " " + in_files
    elif isinstance(in_files, list):
        cmd += " " + " ".join(in_files)
    cmd += " > " + out_file
    exec_sh(cmd)

def count_lines(file_name):
    with open(file_name) as f:
        return sum(1 for _ in f)

def check_equal(var1, var2):
    return var1 == var2

def check_in(var1, var2):
    return var1 in var2

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

