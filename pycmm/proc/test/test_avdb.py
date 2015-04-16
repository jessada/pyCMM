import filecmp
import unittest
from os.path import join as join_path
from os.path import dirname
from pycmm import settings
from pycmm.template import SafeTester
from pycmm.proc.avdb import uniq_avdb


class TestAvdbFunctions(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            )

    def setUp(self):
        self.module_name = 'avdb'

    def __create_db_instance(self):
        return None

    def test_uniq_avdb_1(self):
        """
        to validate if the funcion work with small number of records
        """

        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input')
        out_file = join_path(self.working_dir,
                                'output.avdb')
        uniq_avdb(in_file, out_file)
        exp_file = join_path(self.data_dir,
                                'expected_result')
        self.assertTrue(filecmp.cmp(out_file,
                                    exp_file),
                        "uniq_avdb doesn't work correctly with small amount of avdb records")

    @unittest.skip("obsolete")
    def test_uniq_avdb_2(self):
        """
        to validate if the funcion work with 1.2M number of records
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'in_'+self.current_func_name)
        out_file = join_path(self.working_dir,
                                'out_'+self.current_func_name+'.avdb')
        uniq_avdb(in_file, out_file)
        exp_file = join_path(self.data_dir,
                                'exp_'+self.current_func_name)
        self.assertTrue(filecmp.cmp(out_file, 
                                    exp_file),
                        "uniq_avdb doesn't work correctly with large amount of avdb records")
