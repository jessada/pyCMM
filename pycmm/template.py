import unittest
import sys
import os
import shutil
import subprocess
import inspect
import fileinput
import gc
import string
import tempfile
import datetime
from random import choice
from os.path import join as join_path
from os.path import dirname
from pycmm.utils import mylogger
from pycmm.utils import exec_sh
from pycmm.settings import ENV_TEST_DIR
from pycmm.settings import DEBUG_MODE

ENV_TMPDIR = "TMPDIR"


class pyCMMBase(object):
    """ pyCMM base class """

    def __init__(self, *args, **kwargs):
        self.__time_stamp = datetime.datetime.now()
        self.pkg_root_dir = dirname(dirname(__file__))
        super(pyCMMBase, self).__init__(*args, **kwargs)
        self.__local_scratch_dir = None

    def __str__(self):
        return '<' + self.__class__.__name__ + ' Object> ' + str(self.get_raw_obj_str())

    def get_raw_obj_str(self):
        return "Not yet implemented"

    @property
    def time_stamp(self):
        return self.__time_stamp

    def get_tmp_file_name(self):
        chars = string.ascii_letters
        return "tmp_" + ''.join([choice(chars) for i in range(6)])

    @property
    def local_scratch_dir(self):
        if self.__local_scratch_dir is None:
            self.__local_scratch_dir = tempfile.mkdtemp()
        return self.__local_scratch_dir

    def new_local_tmp_file(self):
        return join_path(self.local_scratch_dir,
                         self.get_tmp_file_name())

    def remove_dir(self, dir_name):
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name, ignore_errors=True)

    def create_dir(self, dir_name):
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

    def delete_file(self, file_name):
        if os.path.exists(file_name):
            if os.path.islink(file_name):
                os.unlink(file_name)
            else:
                os.remove(file_name)

    def concat_files(self, srcs, dst):
        with open(dst, 'w') as fout:
            for line in fileinput.input(srcs):
                fout.write(line)

    def copy_file(self, src, dst):
        if not os.path.isdir(dst):
            self.delete_file(dst)
        if os.path.islink(src):
            linkto = os.readlink(src)
            os.symlink(linkto, dst)
        else:
            shutil.copy(src, dst)

    def copytree(src, dst, symlinks=False, ignore=None):
        for item in os.listdir(src):
            src_item = os.path.join(src, item)
            dst_item = os.path.join(dst, item)
            if os.path.isdir(s):
                shutil.copytree(src_item, dst_item, symlinks, ignore)
            else:
                shutil.copy2(src_item, dst_item)

    def dbg(self, dbg_msg=""):
        if DEBUG_MODE:
            frm = inspect.stack()[1]
            mod = inspect.getmodule(frm[0])
            mylogger.getLogger(mod.__name__)
            mylogger.debug(dbg_msg)

    def info(self, info_msg=""):
        frm = inspect.stack()[1]
        mod = inspect.getmodule(frm[0])
        mylogger.getLogger(mod.__name__)
        mylogger.info(info_msg)

    def warning(self, warning_msg=""):
        frm = inspect.stack()[1]
        mod = inspect.getmodule(frm[0])
        mylogger.getLogger(mod.__name__)
        mylogger.warning(warning_msg)

    def throw(self, err_msg=""):
        frm = inspect.stack()[1]
        mod = inspect.getmodule(frm[0])
        mylogger.getLogger(mod.__name__)
        mylogger.throw(err_msg)

    # having this function built-in for only informative logging purpose
    def exec_sh(self, cmd, silent=False):
        frm = inspect.stack()[1]
        mod = inspect.getmodule(frm[0])
        mylogger.getLogger(mod.__name__)
        return exec_sh(cmd, silent)

    @property
    def current_func_name(self):
        frame = inspect.currentframe(1)
        code  = frame.f_code
        globs = frame.f_globals
        functype = type(lambda: 0)
        funcs = []
        for func in gc.get_referrers(code):
            if type(func) is functype:
                if getattr(func, "func_code", None) is code:
                    if getattr(func, "func_globals", None) is globs:
                        funcs.append(func)
                        if len(funcs) > 1:
                            return None
        return funcs[0].__name__ if funcs else None


class Tester(unittest.TestCase, pyCMMBase):
    """ general pyCMM template for testing """

    individual_debug = False

    def __init__(self, test_module_name, *args, **kwargs):
        self.test_module_name = test_module_name
        super(Tester, self).__init__(*args, **kwargs)
        pyCMMBase.__init__(self)

    def remove_dir(self, dir_name):
        self.assertTrue(dir_name, '"None" is not a valid directory')
        pyCMMBase.remove_dir(self, dir_name)

    def create_dir(self, dir_name):
        self.assertTrue(dir_name, '"None" is not a valid directory')
        pyCMMBase.create_dir(self, dir_name)

    def empty_working_dir(self):
        if not self.individual_debug:
            self.remove_dir(self.working_dir)
        self.create_dir(self.working_dir)

    def remove_working_dir(self):
        if not self.individual_debug:
            self.remove_dir(self.working_dir)

    def set_dir(self):
        working_subdir = "/".join(self.test_module_name.split('.')[:-2])
        working_subdir = join_path(working_subdir,
                                   self.test_module_name.split('.')[-1][5:])
        self.working_dir = os.getenv(ENV_TEST_DIR,
                                     join_path(self.pkg_root_dir,
                                               "tmp"))
        self.working_dir = join_path(self.working_dir,
                                     working_subdir)
        self.working_dir = join_path(join_path(self.working_dir,
                                               self.__class__.__name__[4:]),
                                     self.test_function)
        data_subdir = "/".join(self.test_module_name.split('.')[:-1])
        data_subdir = join_path(data_subdir,
                                "data")
        data_subdir = join_path(data_subdir,
                                self.test_module_name.split('.')[-1][5:])
        self.data_dir = join_path(self.pkg_root_dir,
                                  data_subdir)
        self.data_dir = join_path(join_path(self.data_dir,
                                            self.__class__.__name__[4:]),
                                  self.test_function)

    def init_test(self,
                  test_function,
                  ):
        self.test_function = test_function
        self.set_dir()
        self.empty_working_dir()

    def tearDown(self):
        self.remove_working_dir()


class SafeTester(Tester):
    """

    General template for testing
    that can be run in both dev and production environment
    The purpose of this template is to test general functionality

    """

    def __init__(self, *args, **kwargs):
        super(SafeTester, self).__init__(*args, **kwargs)


class RiskyTester(Tester):
    """

    General template for testing
    that can be run only in devevelopment environment

    The purpose of this template is to test
    if the modules can function properly in real environment

    """

    def __init__(self, test_name):
        Tester.__init__(self, test_name)

    def remove_user_dir(self):
        if not self.individual_debug:
            self.remove_dir(cbv_const.USER_DATA_ROOT)
