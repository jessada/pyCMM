import unittest
from pycmm.utils import is_version
from pycmm.utils.ver import VersionManager
from pycmm.template import SafeTester

class TestVersionManager(SafeTester):

    def __init__(self, methodName):
        super(TestVersionManager, self).__init__(methodName=methodName,
                                                 test_module_name=__name__,
                                                 )

    def setUp(self):
        pass

    def test_pyaml_version(self):
        """ check version of pyaml module """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        vm = VersionManager()
        self.assertTrue(is_version(vm.pyaml_version),
                        "Cannot identify version of pyaml package")

    def tearDown(self):
        self.remove_working_dir()
