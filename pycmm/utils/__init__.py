import subprocess
import datetime
import re
from os.path import dirname
from os.path import realpath
from pycmm.utils import mylogger
from collections import OrderedDict
from collections import Callable

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

def exec_sh(cmd, silent=False):
    mylogger.debug("executing: " + repr(cmd))
    p = subprocess.Popen(cmd,
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         )
    stdout_data, stderr_data = p.communicate()
    return_code = p.returncode
    if not silent:
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

def is_version(s):
    # this allow version number like
    # 2
    # 2.1
    # 2.2.3
    # 2.11.0.4
    result = re.match("^(\d+\.){0,3}(\d+)$", s)
    return result is not None

class DefaultOrderedDict(OrderedDict):
    # Source: http://stackoverflow.com/a/6190500/562769
    def __init__(self, default_factory=None, *args, **kwargs):
        if (default_factory is not None and
           not isinstance(default_factory, Callable)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *args, **kwargs)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))

    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory,
                                               OrderedDict.__repr__(self))

