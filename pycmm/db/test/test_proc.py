import filecmp
import unittest
import vcf
from os.path import join as join_path
from os.path import dirname
from pycmm.template import SafeTester
from pycmm.db.proc import create_empty_vcf


class TestProcFunctions(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            )

    def setUp(self):
        self.module_name = 'proc'

    def __create_db_instance(self):
        return None

    def test_create_empty_vcf_1(self):
        """
        to check if an empty vcf created by copying meta info from 
        a real vcf file is readable
        """

        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                            'input.vcf.gz')
        vcf_reader = vcf.Reader(filename=in_file)
        for record in vcf_reader:
            print record
        self.assertEqual(len(list(vcf_reader)),
                         0,
                         "The file is not empty")

    def test_create_empty_vcf_2(self):
        """ to check if an empty vcf with on vcf version is readable """

#        self.individual_debug = True
        self.init_test(self.current_func_name)
        out_file = join_path(self.working_dir,
                                'output.vcf')
        create_empty_vcf(out_file)
        vcf_reader = vcf.Reader(filename=out_file+'.gz')
        for record in vcf_reader:
            print record
        self.assertEqual(len(list(vcf_reader)),
                         0,
                         "The file is not empty")
#        in_file = join_path(self.data_dir,
#                               'in_'+self.current_func_name+'.vcf.gz')
#        out_file = join_path(self.working_dir,
#                                'out_'+self.current_func_name+'.pavdb')
#        vcf2pavdb(in_file, out_file)
#        exp_file = join_path(self.data_dir,
#                                'exp_'+self.current_func_name)
#        self.assertTrue(filecmp.cmp(out_file,
#                                    exp_file),
#                        "vcf2pavdb doesn't funciton correctly")

#    def test_vcf2pavdb_1(self):
#        """
#        to check if a vcf file can be parsed into the temporary file
#        later to be converted into avdb file
#        """
#
#        self.init_test(self.current_func_name)
#        in_file = join_path(self.data_dir,
#                               'in_'+self.current_func_name+'.vcf.gz')
#        out_file = join_path(self.working_dir,
#                                'out_'+self.current_func_name+'.pavdb')
#        vcf2pavdb(in_file, out_file)
#        exp_file = join_path(self.data_dir,
#                                'exp_'+self.current_func_name)
#        self.assertTrue(filecmp.cmp(out_file, 
#                                    exp_file),
#                        "vcf2pavdb doesn't funciton correctly")
#
#    def test_vcf2pavdb_2(self):
#        """
#        to check if a vcf file can be parsed into the temporary file
#        later to be converted into avdb file
#        """
#
#        self.init_test(self.current_func_name)
#        in_file = join_path(self.data_dir,
#                               'in_'+self.current_func_name+'.vcf.gz')
#        out_file = join_path(self.working_dir,
#                                'out_'+self.current_func_name+'.pavdb')
#        vcf2pavdb(in_file, out_file)
#        exp_file = join_path(self.data_dir,
#                                'exp_'+self.current_func_name)
#        self.assertTrue(filecmp.cmp(out_file, 
#                                    exp_file),
#                        "vcf2pavdb doesn't funciton correctly")
#
#    def test_pavdb2avinputs(self):
#        """
#        to check if convert2annovar.pl can correctly generate avdb files
#        """
#
#        self.init_test(self.current_func_name)
#        in_file = join_path(self.data_dir,
#                               'in_'+self.current_func_name+'.pavdb')
#        out_file_prefix = join_path(self.working_dir,
#                                       'out_'+self.current_func_name)
#        pavdb2avinputs(in_file, out_file_prefix)
#        out_1_file = out_file_prefix + '.dummy1.avinput'
#        out_2_file = out_file_prefix + '.dummy2.avinput'
#        exp_1_file = join_path(self.data_dir,
#                                  'exp_1_' + self.current_func_name)
#        exp_2_file = join_path(self.data_dir,
#                                  'exp_2_' + self.current_func_name)
#        self.assertTrue(filecmp.cmp(out_1_file, 
#                                    exp_1_file),
#                        "pavdb2avinputs doesn't funciton correctly")
#        self.assertTrue(filecmp.cmp(out_2_file, 
#                                    exp_2_file),
#                        "pavdb2avinputs doesn't funciton correctly")
#
#    def test_raw_avdb2key_avdb(self):
#        """
#        to check if the function can correctly convert the raw avdb (avinput)
#        file into avdb file with vcf key
#        """
#
#        self.init_test(self.current_func_name)
#        in_file = join_path(self.data_dir,
#                               'in_'+self.current_func_name+'.avinput')
#        out_file = join_path(self.working_dir,
#                                'out_'+self.current_func_name)
#        raw_avdb2key_avdb(in_file, out_file)
#        exp_file = join_path(self.data_dir,
#                                  'exp_' + self.current_func_name)
#        self.assertTrue(filecmp.cmp(out_file, 
#                                    exp_file),
#                        "raw_avdb2key_avdb doesn't funciton correctly")
#
#    def test_vcf2avdb_1(self):
#        """ to check if vcf file can be correctly converted in avdb file """
#
#        self.init_test(self.current_func_name)
#        in_file = join_path(self.data_dir,
#                               'in_'+self.current_func_name+'.vcf.gz')
#        out_file = join_path(self.working_dir,
#                                'out_'+self.current_func_name+'.avdb')
#        vcf2avdb(in_file, out_file, self.working_dir)
#        exp_file = join_path(self.data_dir,
#                                'exp_'+self.current_func_name)
#        self.assertTrue(filecmp.cmp(out_file, 
#                                    exp_file),
#                        "vcf2avdb doesn't funciton correctly with small vcf file")
#
#    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
#    def test_vcf2avdb_2(self):
#        """ to check if vcf file can be correctly converted in avdb file """
#
#        self.init_test(self.current_func_name)
#        in_file = join_path(self.data_dir,
#                               'in_'+self.current_func_name+'.vcf.gz')
#        out_file = join_path(self.working_dir,
#                                'out_'+self.current_func_name+'.avdb')
#        vcf2avdb(in_file, out_file, self.working_dir)
#        exp_file = join_path(self.data_dir,
#                                'exp_'+self.current_func_name)
#        self.assertTrue(filecmp.cmp(out_file, 
#                                    exp_file),
#                        "vcf2avdb doesn't funciton correctly with small vcf file")
#
#    def test_avdb2bed_1(self):
#        """ to check if avdb2bed can correctly convert in small avdb file """
#
##        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        in_file = join_path(self.data_dir,
#                               'in_'+self.current_func_name+'.avdb')
#        out_file = join_path(self.working_dir,
#                                'out_'+self.current_func_name+'.bed')
#        avdb2bed(in_file, out_file, self.working_dir)
#        exp_file = join_path(self.data_dir,
#                                'exp_'+self.current_func_name)
#        self.assertTrue(filecmp.cmp(out_file, 
#                                    exp_file),
#                        "avdb2bed doesn't funciton correctly with small avdb file")
#
#    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
#    def test_avdb2bed_2(self):
#        """ to check if avdb2bed can correctly convert in small avdb file """
#
##        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        in_file = join_path(self.data_dir,
#                               'in_'+self.current_func_name+'.avdb')
#        out_file = join_path(self.working_dir,
#                                'out_'+self.current_func_name+'.bed')
#        avdb2bed(in_file, out_file, self.working_dir)
#        exp_file = join_path(self.data_dir,
#                                'exp_'+self.current_func_name)
#        self.assertTrue(filecmp.cmp(out_file, 
#                                    exp_file),
#                        "avdb2bed doesn't funciton correctly with large avdb file")
#
