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

    def test_python_pkgs_version(self):
        """ check version of python packages module """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        vm = VersionManager()
        self.assertTrue(is_version(vm.pysam_version),
                        "Cannot identify version of pysam package")
        self.assertTrue(is_version(vm.pyvcf_version),
                        "Cannot identify version of pyvcf package")
        self.assertTrue(is_version(vm.pyaml_version),
                        "Cannot identify version of pyaml package")
        self.assertTrue(is_version(vm.openpyxl_version),
                        "Cannot identify version of openpyxl package")
        self.assertTrue(is_version(vm.xlsxwriter_version),
                        "Cannot identify version of xlsxwriter package")

    def test_non_python_pkgs_version(self):
        """ check version of non-python packages module """
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vm = VersionManager()
        self.assertTrue(vm.gatk_version is not None,
                        "Cannot identify version of GATK")
        self.assertTrue(is_version(vm.gatk_version),
                        "Cannot identify version of GATK")
        self.assertTrue(is_version(vm.plink_version),
                        "Cannot identify version of PLINK")
        self.assertTrue(vm.table_annovar_version > 0,
                        "Cannot identify version of table_annovar")
        self.assertTrue(len(vm.table_annovar_version.split('\n')) == 1,
                        "Cannot identify version of table_annovar")
        self.assertTrue(vm.vcftools_version > 0,
                        "Cannot identify version of vcftools")
        self.assertTrue(len(vm.vcftools_version.split('\n')) == 1,
                        "Cannot identify version of vcftools")

    def tearDown(self):
        self.remove_working_dir()
