import filecmp
from os.path import join as join_path
from os.path import dirname
from pycmm.template import SafeTester
import vcf


class TestPyVCFFunctions(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            )

    def setUp(self):
        self.module_name = 'PyVCF'

    def __create_db_instance(self):
        return None

    def test_cal_zygo_1(self):
        """
        to check if a vcf file can be parsed into the temporary file
        later to be converted into avdb file
        """

        self.init_test(self.current_func_name)

    def test_multiallelic_vcf_1(self):
        """ test if PyVCF can read multiallelic vcf file correctly """
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = vcf.Reader(filename=in_file)

#    def test_write_vcf_1(self):
#        """ test if PyVCF can correctly write VCF file """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        in_file = join_path(self.data_dir,
#                               'input.vcf.gz')
#                               #'in_'+self.current_func_name+'.vcf.gz')
#        out_file = join_path(self.working_dir,
#                               'output.vcf')
##        write_vcf(in_file, out_file)
#
#    def test_load_vcf_1(self):
#        """ test if PyVCF can correctly load VCF file """
#
#        self.init_test(self.current_func_name)
#        in_file = join_path(self.data_dir,
#                               'input.vcf')
#        load_vcf(in_file)
#
#    def test_read_vcf_1(self):
#        """ test if PyVCF can correctly read VCF file """
#
#        self.init_test(self.current_func_name)
#        in_file = join_path(self.data_dir,
#                               'input.vcf.gz')
#        read_vcf(in_file)
