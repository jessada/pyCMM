import unittest
import filecmp
import datetime
from os.path import exists as path_exists
from os.path import join as join_path
from os.path import dirname
from os.path import isdir
from pycmm import settings
from pycmm.flow.mutrep import MutRepPipeline
from pycmm.flow.mutrep import create_jobs_setup_file
from pycmm.settings import DFLT_ANNOVAR_DB_FOLDER
from pycmm.settings import FAST_PROJECT_CODE
from pycmm.settings import SLOW_PROJECT_CODE
from pycmm.template import SafeTester
from pycmm.utils import mylogger

DFLT_ANNOVAR_TEST_DB_FOLDER = DFLT_ANNOVAR_DB_FOLDER 
DFLT_ANNOVAR_TEST_DB_NAMES = "refGene" 
DFLT_ANNOVAR_TEST_DB_OPS = "g" 
DFLT_ANNOVAR_TEST_DB_NAMES += ",cytoBand" 
DFLT_ANNOVAR_TEST_DB_OPS += ",r" 
DFLT_ANNOVAR_TEST_DB_NAMES += ",genomicSuperDups" 
DFLT_ANNOVAR_TEST_DB_OPS += ",r" 
DFLT_ANNOVAR_TEST_DB_NAMES += ",exac03" 
DFLT_ANNOVAR_TEST_DB_OPS += ",f" 
DFLT_ANNOVAR_TEST_DB_NAMES += ",custom" 
DFLT_ANNOVAR_TEST_DB_OPS += ",f" 

class TestAnnovar(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            )

    def setUp(self):
        self.module_name = 'annovar'

    def __create_jobs_setup_file(self,
                                 vcf_tabix_file,
                                 dataset_name=None,
                                 project_code=FAST_PROJECT_CODE,
                                 annovar_human_db_dir=DFLT_ANNOVAR_TEST_DB_FOLDER,
                                 annovar_buildver="hg19",
                                 annovar_db_names=DFLT_ANNOVAR_TEST_DB_NAMES,
                                 annovar_db_ops=DFLT_ANNOVAR_TEST_DB_OPS,
                                 annovar_nastring=".",
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        if dataset_name is None:
            dataset_name = self.test_function
        create_jobs_setup_file(dataset_name=dataset_name,
                               project_out_dir=self.working_dir,
                               vcf_tabix_file=vcf_tabix_file,
                               project_code=project_code,
                               annovar_human_db_dir=annovar_human_db_dir,
                               annovar_buildver=annovar_buildver,
                               annovar_db_names=annovar_db_names,
                               annovar_db_ops=annovar_db_ops,
                               annovar_nastring=annovar_nastring,
                               out_jobs_setup_file=jobs_setup_file,
                               )
        return jobs_setup_file

#    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_table_annovar_offline(self):
#        """ test offline version (w/o slurm) of table_annovar """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        vcf_tabix_file = join_path(self.data_dir,
#                                   "input.vcf.gz")
#        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
#                                                        project_code=FAST_PROJECT_CODE)
#        pl = MutRepPipeline(jobs_setup_file)
#        cmd = pl.table_annovar()
#        mylogger.debug(cmd)
#
    def tearDown(self):
        self.remove_working_dir()
