import unittest
import os
import filecmp
from pycmm import settings
from pycmm.utils.test.template import SafeUtilsTester
from pycmm.utils import get_file_prefix
from pycmm.utils import concat_files


class TestFunctions(SafeUtilsTester):


    def __init__(self, test_name):
        SafeUtilsTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'NoClass'

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
        in_file_wildcards = os.path.join(self.data_dir,
                                         'in_'+self.current_func_name+'*')
        out_file = os.path.join(self.working_dir,
                                'out_'+self.current_func_name)
        concat_files(in_file_wildcards, out_file)
        exp_file = os.path.join(self.data_dir,
                                'exp_'+self.current_func_name)
        self.assertTrue(filecmp.cmp(out_file, 
                                    exp_file),
                        "concat_files cannot concat small avinput files correctly")

    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_concat_files_2(self):
        """ to check if 400K-lined avinput files can be concantenated """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file_wildcards = os.path.join(self.data_dir,
                                         'in_'+self.current_func_name+'*')
        out_file = os.path.join(self.working_dir,
                                'out_'+self.current_func_name)
        concat_files(in_file_wildcards, out_file)
        exp_file = os.path.join(self.data_dir,
                                'exp_'+self.current_func_name)
        self.assertTrue(filecmp.cmp(out_file, 
                                    exp_file),
                        "concat_files cannot concat small avinput files correctly")

    def tearDown(self):
        self.remove_working_dir()
