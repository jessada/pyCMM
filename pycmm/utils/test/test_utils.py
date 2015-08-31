import unittest
import filecmp
from os.path import join as join_path
from os.path import dirname
from pycmm import settings
from pycmm.template import SafeTester
from pycmm.utils import mylogger
from pycmm.utils import get_file_prefix
from pycmm.utils import concat_files


class TestFunctions(SafeTester):


    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            test_module_name=__name__,
                            )

    def setUp(self):
        mylogger.getLogger(__name__)
        pass

    def __create_db_instance(self):
        return None

    def test_get_file_prefix(self):
        """ to check if the file prefix is correctly retrieved """

        self.init_test(self.current_func_name)
        self.assertEqual(get_file_prefix("logvcf2avdb.log", ".log"), "logvcf2avdb", "Incorrect file prefix")
        self.assertEqual(get_file_prefix("/home/jessada/vcf2avdb.log", ".log"), "/home/jessada/vcf2avdb", "Incorrect file prefix")
        self.assertEqual(get_file_prefix("/home/jessada/abcab", "ab"), "/home/jessada/abc", "Incorrect file prefix")
        self.assertEqual(get_file_prefix("/home/jessada/abcabab", "ab"), "/home/jessada/abcab", "Incorrect file prefix")
        self.assertEqual(get_file_prefix("./test.txt", ".txt"), "./test", "Incorrect file prefix")

    def test_concat_files_1(self):
        """ to check if small avinput files can be concantenated """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file_wildcards = join_path(self.data_dir,
                                         'input_*')
        out_file = join_path(self.working_dir,
                                'out_'+self.current_func_name)
        concat_files(in_file_wildcards, out_file)
        exp_file = join_path(self.data_dir,
                                'expected_result')
        self.assertTrue(filecmp.cmp(out_file,
                                    exp_file),
                        "concat_files cannot concat small avinput files correctly")

    def tearDown(self):
        self.remove_working_dir()
