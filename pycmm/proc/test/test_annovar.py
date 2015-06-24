import unittest
import filecmp
import datetime
from os.path import exists as path_exists
from os.path import join as join_path
from os.path import dirname
from os.path import isdir
from pycmm import settings
from pycmm.template import SafeTester
from pycmm.utils import mylogger

FAST_PROJECT_CODE = 'b2011158'
SLOW_PROJECT_CODE = 'b2011097'

class TestGATKBPPipeline(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            )

    def setUp(self):
        self.module_name = 'gatkbp'

#    def __create_jobs_setup_file(self,
#                                 dataset_name=None,
#                                 sample_group='test_group',
#                                 project_code=SLOW_PROJECT_CODE,
#                                 variants_calling="YES",
#                                 targets_interval_list=None,
#                                 dataset_usage_mail="NO",
#                                 sample_usage_mail={},
#                                 ):
#        jobs_setup_file = join_path(self.working_dir,
#                                    self.test_function+'_jobs_setup.txt')
#        if dataset_name is None:
#            dataset_name = self.test_function
#        known_indels = []
#        known_indels.append(join_path(self.data_dir,
#                                      '1000G_phase1.indels.b37.vcf'))
#        known_indels.append(join_path(self.data_dir,
#                                      'Mills_and_1000G_gold_standard.indels.b37.vcf'))
#        dbsnp_file = join_path(self.data_dir,
#                               'dbsnp_138.b37.vcf')
#        reference_file = join_path(self.data_dir,
#                                   'ref.fa')
#        create_jobs_setup_file(dataset_name=dataset_name,
#                               sample_group=sample_group,
#                               project_code=project_code,
#                               reference_file=reference_file,
#                               project_out_dir=self.working_dir,
#                               samples_root_dir=self.data_dir,
#                               known_indels_file=known_indels,
#                               dbsnp_file=dbsnp_file,
#                               targets_interval_list=targets_interval_list,
#                               dataset_usage_mail=dataset_usage_mail,
#                               sample_usage_mail=sample_usage_mail,
#                               out_jobs_setup_file=jobs_setup_file,
#                               )
#        return jobs_setup_file
#
    def tearDown(self):
        self.remove_working_dir()

class TestFunctions(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            )

    def setUp(self):
        self.module_name = 'gatkbp'


    def tearDown(self):
        self.remove_working_dir()
