import os
import filecmp
import unittest
from pycmm import settings
from pycmm.proc.test.template import SafeProcTester
from pycmm.proc.converter import vcf2pavdb 
from pycmm.proc.converter import pavdb2avinputs
from pycmm.proc.converter import raw_avdb2key_avdb
from pycmm.proc.converter import vcf2avdb 
from pycmm.proc.converter import avdb2bed


class TestVcfFunctions(SafeProcTester):

    def __init__(self, test_name):
        SafeProcTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'NoClass'

    def __create_db_instance(self):
        return None

    def test_vcf2pavdb_1(self):
        """
        to check if a vcf file can be parsed into the temporary file
        later to be converted into avdb file
        """

        self.init_test(self.current_func_name)
        in_file = os.path.join(self.data_dir,
                               'in_'+self.current_func_name+'.vcf.gz')
        out_file = os.path.join(self.working_dir,
                                'out_'+self.current_func_name+'.pavdb')
        vcf2pavdb(in_file, out_file)
        exp_file = os.path.join(self.data_dir,
                                'exp_'+self.current_func_name)
        self.assertTrue(filecmp.cmp(out_file, 
                                    exp_file),
                        "vcf2pavdb doesn't funciton correctly")

    def test_vcf2pavdb_2(self):
        """
        to check if a vcf file can be parsed into the temporary file
        later to be converted into avdb file
        """

        self.init_test(self.current_func_name)
        in_file = os.path.join(self.data_dir,
                               'in_'+self.current_func_name+'.vcf.gz')
        out_file = os.path.join(self.working_dir,
                                'out_'+self.current_func_name+'.pavdb')
        vcf2pavdb(in_file, out_file)
        exp_file = os.path.join(self.data_dir,
                                'exp_'+self.current_func_name)
        self.assertTrue(filecmp.cmp(out_file, 
                                    exp_file),
                        "vcf2pavdb doesn't funciton correctly")

    def test_pavdb2avinputs(self):
        """
        to check if convert2annovar.pl can correctly generate avdb files
        """

        self.init_test(self.current_func_name)
        in_file = os.path.join(self.data_dir,
                               'in_'+self.current_func_name+'.pavdb')
        out_file_prefix = os.path.join(self.working_dir,
                                       'out_'+self.current_func_name)
        pavdb2avinputs(in_file, out_file_prefix)
        out_1_file = out_file_prefix + '.dummy1.avinput'
        out_2_file = out_file_prefix + '.dummy2.avinput'
        exp_1_file = os.path.join(self.data_dir,
                                  'exp_1_' + self.current_func_name)
        exp_2_file = os.path.join(self.data_dir,
                                  'exp_2_' + self.current_func_name)
        self.assertTrue(filecmp.cmp(out_1_file, 
                                    exp_1_file),
                        "pavdb2avinputs doesn't funciton correctly")
        self.assertTrue(filecmp.cmp(out_2_file, 
                                    exp_2_file),
                        "pavdb2avinputs doesn't funciton correctly")

    def test_raw_avdb2key_avdb(self):
        """
        to check if the function can correctly convert the raw avdb (avinput)
        file into avdb file with vcf key
        """

        self.init_test(self.current_func_name)
        in_file = os.path.join(self.data_dir,
                               'in_'+self.current_func_name+'.avinput')
        out_file = os.path.join(self.working_dir,
                                'out_'+self.current_func_name)
        raw_avdb2key_avdb(in_file, out_file)
        exp_file = os.path.join(self.data_dir,
                                  'exp_' + self.current_func_name)
        self.assertTrue(filecmp.cmp(out_file, 
                                    exp_file),
                        "raw_avdb2key_avdb doesn't funciton correctly")

    def test_vcf2avdb_1(self):
        """ to check if vcf file can be correctly converted in avdb file """

        self.init_test(self.current_func_name)
        in_file = os.path.join(self.data_dir,
                               'in_'+self.current_func_name+'.vcf.gz')
        out_file = os.path.join(self.working_dir,
                                'out_'+self.current_func_name+'.avdb')
        vcf2avdb(in_file, out_file, self.working_dir)
        exp_file = os.path.join(self.data_dir,
                                'exp_'+self.current_func_name)
        self.assertTrue(filecmp.cmp(out_file, 
                                    exp_file),
                        "vcf2avdb doesn't funciton correctly with small vcf file")

    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_vcf2avdb_2(self):
        """ to check if vcf file can be correctly converted in avdb file """

        self.init_test(self.current_func_name)
        in_file = os.path.join(self.data_dir,
                               'in_'+self.current_func_name+'.vcf.gz')
        out_file = os.path.join(self.working_dir,
                                'out_'+self.current_func_name+'.avdb')
        vcf2avdb(in_file, out_file, self.working_dir)
        exp_file = os.path.join(self.data_dir,
                                'exp_'+self.current_func_name)
        self.assertTrue(filecmp.cmp(out_file, 
                                    exp_file),
                        "vcf2avdb doesn't funciton correctly with small vcf file")

    def test_avdb2bed_1(self):
        """ to check if avdb2bed can correctly convert in small avdb file """

#        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = os.path.join(self.data_dir,
                               'in_'+self.current_func_name+'.avdb')
        out_file = os.path.join(self.working_dir,
                                'out_'+self.current_func_name+'.bed')
        avdb2bed(in_file, out_file, self.working_dir)
        exp_file = os.path.join(self.data_dir,
                                'exp_'+self.current_func_name)
        self.assertTrue(filecmp.cmp(out_file, 
                                    exp_file),
                        "avdb2bed doesn't funciton correctly with small avdb file")

    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_avdb2bed_2(self):
        """ to check if avdb2bed can correctly convert in small avdb file """

#        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = os.path.join(self.data_dir,
                               'in_'+self.current_func_name+'.avdb')
        out_file = os.path.join(self.working_dir,
                                'out_'+self.current_func_name+'.bed')
        avdb2bed(in_file, out_file, self.working_dir)
        exp_file = os.path.join(self.data_dir,
                                'exp_'+self.current_func_name)
        self.assertTrue(filecmp.cmp(out_file, 
                                    exp_file),
                        "avdb2bed doesn't funciton correctly with large avdb file")
