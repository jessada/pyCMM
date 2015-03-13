import os
import filecmp
from pycmm.proc.test.template import SafeProcTester
from pycmm.proc.vcf import cal_zygo 


class TestVcfFunctions(SafeProcTester):

    def __init__(self, test_name):
        SafeProcTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'NoClass'

    def __create_db_instance(self):
        return None

    def test_cal_zygo_1(self):
        """
        to check if a vcf file can be parsed into the temporary file
        later to be converted into avdb file
        """

        self.init_test(self.current_func_name)
#        in_file = os.path.join(self.data_dir,
#                               'in_'+self.current_func_name+'.vcf.gz')
#        out_file = os.path.join(self.working_dir,
#                                'out_'+self.current_func_name+'.pavdb')
#        vcf2pavdb(in_file, out_file)
#        exp_file = os.path.join(self.data_dir,
#                                'exp_'+self.current_func_name)
#        self.assertTrue(filecmp.cmp(out_file, 
#                                    exp_file),
#                        "vcf2pavdb doesn't funciton correctly")
