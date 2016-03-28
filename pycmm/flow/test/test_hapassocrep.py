# This Python file uses the following encoding: utf-8
import unittest
from os.path import join as join_path
from pycmm.settings import FULL_SYSTEM_TEST
from pycmm.template import SafeTester
from pycmm.flow.plink import create_jobs_setup_file
from pycmm.flow.hapassocrep import HapAssocRepPipeline

class TestHapAssocRepPipeline(SafeTester):

    def __init__(self, methodName):
        super(TestHapAssocRepPipeline, self).__init__(methodName=methodName,
                                                      test_module_name=__name__,
                                                      )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self,
                                 project_name=None,
                                 project_out_dir=None,
                                 input_file_prefix=None,
                                 input_binary=None,
                                 input_dna_regions=None,
                                 phenotype_file=None,
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        if project_name is None:
            project_name = self.test_function
        if project_out_dir is None:
            project_out_dir = self.working_dir
        if input_file_prefix is None:
            input_file_prefix = join_path(self.data_dir,
                                          "input")
        if input_binary is None:
            input_binary = True
        if phenotype_file is None:
            phenotype_file = join_path(self.data_dir,
                                       "phenotype.txt")
        create_jobs_setup_file(project_name=project_name,
                               project_out_dir=project_out_dir,
                               input_file_prefix=input_file_prefix,
                               phenotype_file=phenotype_file,
                               input_binary=input_binary,
                               input_dna_regions=input_dna_regions,
                               )
        return jobs_setup_file

    def test_gen_report_1(self):
        """ test generating basic haplotype associaiton study report """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        hap_assoc_file = join_path(self.data_dir,
                                   "input.assoc.hap")
        snp_info_file = join_path(self.data_dir,
                                  "input.snp.info")
        pl = HapAssocRepPipeline(hap_assoc_file=hap_assoc_file,
                                 snp_info_file=snp_info_file,
                                 jobs_setup_file=jobs_setup_file,
                                 )
        rpt_file = pl.gen_report()
        self.dbg("************** not yet tested *****************")

    def tearDown(self):
        self.remove_working_dir()
