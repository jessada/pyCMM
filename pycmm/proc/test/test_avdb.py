import os
import filecmp
from pycmm.proc.test.template import SafeProcTester
from pycmm.proc.avdb import uniq_avdb


class TestAvdbFunctions(SafeProcTester):

    def __init__(self, test_name):
        SafeProcTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'NoClass'

    def __create_db_instance(self):
        return None

    def test_uniq_avdb_1(self):
        """
        to validate if the funcion work with small number of records
        """

        self.init_test(self.current_func_name)
        in_file = os.path.join(self.data_dir,
                               'in_'+self.current_func_name)
        out_file = os.path.join(self.working_dir,
                                'out_'+self.current_func_name+'.avdb')
        uniq_avdb(in_file, out_file)
        exp_file = os.path.join(self.data_dir,
                                'exp_'+self.current_func_name)
        self.assertTrue(filecmp.cmp(out_file, 
                                    exp_file),
                        "uniq_avdb doesn't work correctly with small amount of avdb records")

    def test_uniq_avdb_2(self):
        """
        to validate if the funcion work with 1.2M number of records
        """

        self.init_test(self.current_func_name)
        in_file = os.path.join(self.data_dir,
                               'in_'+self.current_func_name)
        out_file = os.path.join(self.working_dir,
                                'out_'+self.current_func_name+'.avdb')
        uniq_avdb(in_file, out_file)
        exp_file = os.path.join(self.data_dir,
                                'exp_'+self.current_func_name)
        self.assertTrue(filecmp.cmp(out_file, 
                                    exp_file),
                        "uniq_avdb doesn't work correctly with large amount of avdb records")
