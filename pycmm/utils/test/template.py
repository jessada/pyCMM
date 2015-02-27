import unittest
import os
import pycmm.settings as pycmm_settings
from pycmm.template import SafeTester
from pycmm.template import RiskyTester


class SafeUtilsTester(SafeTester):
    """ General template for "utils" testing """

    def __init__(self, test_name):
        SafeTester.__init__(self, test_name)

    def set_dir(self):
        self.working_dir = os.path.join(os.path.join(os.path.join(os.path.dirname(__file__),
                                                                  'tmp'),
                                                     self.test_class),
                                        self.test_function)
        self.data_dir = os.path.join(os.path.join(os.path.dirname(__file__),
                                                  'data'),
                                     self.test_class)


class RiskyUtilsTester(RiskyTester):
    """ General template for "utils" testing """

    def __init__(self, test_name):
        RiskyTester.__init__(self, test_name)

    def set_dir(self):
        self.working_dir = linkana_settings.linkana_WORKING_DIR
        self.data_dir = os.path.join(os.path.join(os.path.dirname(__file__),
                                                  'big_data'),
                                     self.test_class)
