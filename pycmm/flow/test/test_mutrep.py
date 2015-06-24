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
from pycmm.flow.mutrep import MutRepPipeline

from pycmm.proc.annovar import ANNOVAR_PARAMS_INPUT_FILE_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_DB_FOLDER_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_BUILDVER_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_OUT_PREFIX_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_DB_LIST_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_DB_NAME_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_DB_OP_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_NASTRING_KEY

FAST_PROJECT_CODE = 'b2011158'
SLOW_PROJECT_CODE = 'b2011097'

class TestMutRepPipeline(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            )

    def setUp(self):
        self.module_name = 'mutrep'

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
    def test_load_jobs_info_1(self):
        """ test if job configurations are loaded correctly """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        jobs_setup_file = join_path(self.data_dir,
                                    self.test_function+'.txt')
        pl = MutRepPipeline(jobs_setup_file)
        exp_dataset_name = "test-dataset"
        self.assertEqual(pl.dataset_name,
                         exp_dataset_name,
                         "GATKBPPipeline cannot correctly read meta info 'dataset name' from jobs setup file")
        self.assertEqual(pl.project_code,
                         "b1234567",
                         "GATKBPPipeline cannot correctly read meta info 'project code' from jobs setup file")
        self.assertEqual(pl.input_vcf_tabix,
                         "/home/jesthu/test.vcf.gz",
                         "GATKBPPipeline cannot correctly read meta info 'input vcf tabix' from jobs setup file")
        exp_jobs_report_file = join_path(pl.output_dir,
                                         exp_dataset_name+"_rpt.txt")
        self.assertEqual(pl.jobs_report_file,
                         exp_jobs_report_file,
                         "GATKBPPipeline cannot correctly read meta info 'jobs report file' from jobs setup file")
        self.assertEqual(pl.output_dir,
                         "/glob/jessada/private/projects/",
                         "GATKBPPipeline cannot correctly read meta info 'output dir' from jobs setup file")
        self.assertEqual(pl.annovar_config[ANNOVAR_PARAMS_INPUT_FILE_KEY],
                         "/home/jesthu/test.vcf.gz",
                         "GATKBPPipeline cannot correctly read meta info 'annovar input file' from jobs setup file")
        self.assertEqual(pl.annovar_config[ANNOVAR_PARAMS_DB_FOLDER_KEY],
                         "/home/jesthu/annovar",
                         "GATKBPPipeline cannot correctly read meta info 'annovar db folder' from jobs setup file")
        self.assertEqual(pl.annovar_config[ANNOVAR_PARAMS_BUILDVER_KEY],
                         "hg19",
                         "GATKBPPipeline cannot correctly read meta info 'annovar buildver' from jobs setup file")
        exp_annovar_out_prefix = join_path(pl.rpts_out_dir,
                                           exp_dataset_name)
        self.assertEqual(pl.annovar_config[ANNOVAR_PARAMS_OUT_PREFIX_KEY],
                         exp_annovar_out_prefix,
                         "GATKBPPipeline cannot correctly read meta info 'annovar out prefix' from jobs setup file")
        self.assertEqual(pl.annovar_config[ANNOVAR_PARAMS_NASTRING_KEY],
                         ".",
                         "GATKBPPipeline cannot correctly read meta info 'annovar nastring' from jobs setup file")
        self.assertEqual(pl.annovar_config[ANNOVAR_PARAMS_DB_LIST_KEY][1][ANNOVAR_PARAMS_DB_NAME_KEY],
                         "cytoBand",
                         "GATKBPPipeline cannot correctly read meta info 'annovar db name' from jobs setup file")
        self.assertEqual(pl.annovar_config[ANNOVAR_PARAMS_DB_LIST_KEY][4][ANNOVAR_PARAMS_DB_NAME_KEY],
                         "1000g2014oct_eur",
                         "GATKBPPipeline cannot correctly read meta info 'annovar db name' from jobs setup file")
        self.assertEqual(pl.annovar_config[ANNOVAR_PARAMS_DB_LIST_KEY][2][ANNOVAR_PARAMS_DB_OP_KEY],
                         "r",
                         "GATKBPPipeline cannot correctly read meta info 'annovar db name' from jobs setup file")
        self.assertEqual(pl.annovar_config[ANNOVAR_PARAMS_DB_LIST_KEY][5][ANNOVAR_PARAMS_DB_OP_KEY],
                         "f",
                         "GATKBPPipeline cannot correctly read meta info 'annovar db name' from jobs setup file")

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
