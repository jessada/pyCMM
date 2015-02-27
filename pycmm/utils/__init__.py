import subprocess
from pycmm.utils import mylogger

def get_file_prefix(file_name, file_ext):
    if file_name.endswith(file_ext):
        return file_name[:-len(file_ext)]
    return file_name

def log_file_with_time_stamp(raw_file, time_stamp):
    log_file = get_file_prefix(raw_file, '.log')
    log_file += '_'
    log_file += time_stamp.strftime("%Y%m%d%H%M%S")
    log_file += '.log'
    return log_file

def exec_sh(cmd):
    p = subprocess.Popen(cmd, shell=True)
    mylogger.info("executing: " + repr(cmd))
    error = p.wait()
    if error:
        mylogger.throw("Error found during execute command '%s' with error code %d" % (cmd, error))

def concat_files(in_files,
                 out_file):
    cmd = "cat"
    if isinstance(in_files, str):
        cmd += " " + in_files
    elif isinstance(in_files, list):
        cmd += " " + " ".join(in_files)
    cmd += " > " + out_file
    exec_sh(cmd)

