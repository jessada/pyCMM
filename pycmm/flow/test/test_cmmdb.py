import unittest
import filecmp
import datetime
from os.path import exists as path_exists
from os.path import join as join_path
from os.path import dirname
from os.path import isdir
from pycmm.template import SafeTester
from pycmm.settings import SLURM_CMMDB_TEST
from pycmm.settings import FULL_SYSTEM_TEST
from pycmm.settings import DFLT_ANNOVAR_DB_FOLDER
from pycmm.settings import DFLT_ANNOVAR_DB_NAMES
from pycmm.settings import DFLT_ANNOVAR_DB_OPS
from pycmm.settings import DFLT_CMMDB_ALLOC_TIME
from pycmm.settings import FAST_PROJECT_CODE
from pycmm.settings import SLOW_PROJECT_CODE
from pycmm.flow.cmmdb import CMMDBPipeline
from pycmm.flow.cmmdb import create_jobs_setup_file
from pycmm.flow.cmmdb import JOBS_SETUP_SAMPLE_INFOS_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_FAMILY_ID_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_MEMBERS_LIST_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_SAMPLE_ID_KEY
from pycmm.proc.annovarlib import ANNOVAR_PARAMS_INPUT_FILE_KEY
from pycmm.proc.annovarlib import ANNOVAR_PARAMS_DB_FOLDER_KEY
from pycmm.proc.annovarlib import ANNOVAR_PARAMS_BUILDVER_KEY
from pycmm.proc.annovarlib import ANNOVAR_PARAMS_OUT_PREFIX_KEY
from pycmm.proc.annovarlib import ANNOVAR_PARAMS_DB_NAMES_KEY
from pycmm.proc.annovarlib import ANNOVAR_PARAMS_DB_OPS_KEY
from pycmm.proc.annovarlib import ANNOVAR_PARAMS_NASTRING_KEY

DFLT_ANNOVAR_TEST_DB_FOLDER = DFLT_ANNOVAR_DB_FOLDER
DFLT_ANNOVAR_TEST_DB_NAMES = "refGene"
DFLT_ANNOVAR_TEST_DB_OPS = "g"
DFLT_ANNOVAR_TEST_DB_NAMES += ",cytoBand"
DFLT_ANNOVAR_TEST_DB_OPS += ",r"
DFLT_ANNOVAR_TEST_DB_NAMES += ",genomicSuperDups"
DFLT_ANNOVAR_TEST_DB_OPS += ",r"


class TestCMMDBPipeline(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            test_module_name=__name__,
                            )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self,
                                 vcf_tabix_file,
                                 dataset_name=None,
                                 project_code=SLOW_PROJECT_CODE,
                                 db_alloc_time=DFLT_CMMDB_ALLOC_TIME,
                                 db_region="18",
                                 sample_infos=None,
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
                               db_region=db_region,
                               sample_infos=sample_infos,
                               project_code=project_code,
                               db_alloc_time=db_alloc_time,
                               annovar_human_db_dir=annovar_human_db_dir,
                               annovar_buildver=annovar_buildver,
                               annovar_db_names=annovar_db_names,
                               annovar_db_ops=annovar_db_ops,
                               annovar_nastring=annovar_nastring,
                               out_jobs_setup_file=jobs_setup_file,
                               )
        return jobs_setup_file

    @unittest.skipUnless(isdir("/proj/b2011117"), "This can only run in UPPMAX")
    def test_load_jobs_info_1(self):
        """ test if job configurations are loaded correctly """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        dummy_vcf_tabix_file = "/path/to/vcf_tabix_file"
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=dummy_vcf_tabix_file)
        pl = CMMDBPipeline(jobs_setup_file)
        exp_dataset_name = self.test_function
        self.assertEqual(pl.dataset_name,
                         exp_dataset_name,
                         "CMMDBPipeline cannot correctly read meta info 'dataset name' from jobs setup file")
        self.assertEqual(pl.project_code,
                         SLOW_PROJECT_CODE,
                         "CMMDBPipeline cannot correctly read meta info 'project code' from jobs setup file")
        self.assertEqual(pl.db_alloc_time,
                         DFLT_CMMDB_ALLOC_TIME,
                         "CMMDBPipeline cannot correctly read meta info 'db allocation time' from jobs setup file")
        self.assertEqual(pl.input_vcf_tabix,
                         dummy_vcf_tabix_file,
                         "CMMDBPipeline cannot correctly read meta info 'input vcf tabix' from jobs setup file")
        self.assertEqual(pl.db_region,
                         "18",
                         "CMMDBPipeline cannot correctly read meta info 'vcf region' from jobs setup file")
        self.assertEqual(pl.samples_list,
                         None,
                         "CMMDBPipeline cannot correctly read meta info 'samples list' from jobs setup file")
        exp_jobs_report_file = join_path(pl.output_dir,
                                         exp_dataset_name+"_rpt.txt")
        self.assertEqual(pl.jobs_report_file,
                         exp_jobs_report_file,
                         "CMMDBPipeline cannot correctly read meta info 'jobs report file' from jobs setup file")
        self.assertEqual(pl.output_dir,
                         self.working_dir,
                         "CMMDBPipeline cannot correctly read meta info 'output dir' from jobs setup file")
        self.assertEqual(pl.annovar_config.input_file,
                         dummy_vcf_tabix_file,
                         "CMMDBPipeline cannot corretly read annovar configs 'annovar input file' from jobs setup file")
        self.assertEqual(pl.annovar_config.db_folder,
                         DFLT_ANNOVAR_TEST_DB_FOLDER,
                         "CMMDBPipeline cannot corretly read annovar configs 'annovar db folder' from jobs setup file")
        self.assertEqual(pl.annovar_config.buildver,
                         "hg19",
                         "CMMDBPipeline cannot corretly read annovar configs 'annovar buildver' from jobs setup file")
        exp_annovar_out_prefix = join_path(pl.data_out_dir,
                                           exp_dataset_name)
        self.assertEqual(pl.annovar_config.out_prefix,
                         exp_annovar_out_prefix,
                         "CMMDBPipeline cannot corretly read annovar configs 'annovar out prefix' from jobs setup file")
        self.assertEqual(pl.annovar_config.nastring,
                         ".",
                         "CMMDBPipeline cannot corretly read annovar configs 'annovar nastring' from jobs setup file")
        self.assertEqual(pl.annovar_config.protocols,
                         DFLT_ANNOVAR_TEST_DB_NAMES,
                         "CMMDBPipeline cannot corretly read annovar configs 'annovar db names' from jobs setup file")
        self.assertEqual(pl.annovar_config.operations,
                         DFLT_ANNOVAR_TEST_DB_OPS,
                         "CMMDBPipeline cannot corretly read annovar configs 'annovar db ops' from jobs setup file")
        exp_annotated_vcf = exp_annovar_out_prefix + ".hg19_multianno.vcf"
        self.assertEqual(pl.annovar_config.annotated_vcf,
                         exp_annotated_vcf,
                         "CMMDBPipeline cannot corretly determine annovar configs 'annovar annotated vcf' file")

    def test_load_jobs_info_2(self):
        """ test if modified job configurations are loaded correctly """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        dummy_vcf_tabix_file = "/path/to/vcf_tabix_file"
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=dummy_vcf_tabix_file,
                                                        project_code=None,
                                                        db_alloc_time="120:00:00",
                                                        db_region=None,
                                                        sample_infos="Co-441,Co-666",
                                                        )
        pl = CMMDBPipeline(jobs_setup_file)
        exp_dataset_name = self.test_function
        self.assertEqual(pl.project_code,
                         None,
                         "CMMDBPipeline cannot correctly read meta info 'project code' from jobs setup file")
        self.assertEqual(pl.db_alloc_time,
                         "120:00:00",
                         "CMMDBPipeline cannot correctly read meta info 'db allocation time' from jobs setup file")
        self.assertEqual(pl.db_region,
                         None,
                         "CMMDBPipeline cannot correctly read meta info 'vcf region' from jobs setup file")
        self.assertEqual(pl.samples_list,
                         ["Co-441","Co-666"],
                         "CMMDBPipeline cannot correctly read meta info 'samples list' from jobs setup file")
        self.assertEqual(pl.family_infos,
                         None,
                         "CMMDBPipeline cannot correctly read meta info 'families list' from jobs setup file")

    def test_load_jobs_info_3(self):
        """ test if families structures are loaded correctly """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        dummy_vcf_tabix_file = "/path/to/vcf_tabix_file"
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=dummy_vcf_tabix_file,
                                                        sample_infos="242:Co-441:Co-666:Co-771,6067:Br-397:884-99D",
                                                        )
        pl = CMMDBPipeline(jobs_setup_file)
        self.assertEqual(pl.samples_list,
                         ["Co-441","Co-666","Co-771","Br-397","884-99D"],
                         "CMMDBPipeline cannot correctly read 'samples list' from jobs setup file")
        self.assertEqual(len(pl.family_infos),
                         2,
                         "CMMDBPipeline cannot correctly read meta info 'families list' from jobs setup file")
        fam_info_1 = pl.family_infos['242']
        fam_info_2 = pl.family_infos['6067']
        self.assertEqual(fam_info_2.fam_id,
                         "6067",
                         "CMMDBPipeline cannot correctly read meta info 'families list' from jobs setup file")
        self.assertEqual(len(fam_info_1.members),
                         3,
                         "CMMDBPipeline cannot correctly read meta info 'families list' from jobs setup file")
        self.assertEqual(len(fam_info_2.members),
                         2,
                         "CMMDBPipeline cannot correctly read meta info 'families list' from jobs setup file")
        self.assertEqual(fam_info_2.members[1].sample_id,
                         "884-99D",
                         "CMMDBPipeline cannot correctly read meta info 'families list' from jobs setup file")

    def test_load_jobs_info_4(self):
        """ test if can load sample infos as a file """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        dummy_vcf_tabix_file = "/path/to/vcf_tabix_file"
        sample_infos_file = join_path(self.data_dir,
                                      "sample_infos")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=dummy_vcf_tabix_file,
                                                        sample_infos=sample_infos_file,
                                                        )
        pl = CMMDBPipeline(jobs_setup_file)
        self.assertEqual(pl.samples_list,
                         ["Co-35", "Co-37", "Co-135", "Co-131", "1322-11D"],
                         "CMMDBPipeline cannot correctly read 'samples list' from jobs setup file")

    def test_load_jobs_info_5(self):
        """ test if can load sample infos in family structur format as a file """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        dummy_vcf_tabix_file = "/path/to/vcf_tabix_file"
        sample_infos_file = join_path(self.data_dir,
                                      "sample_infos")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=dummy_vcf_tabix_file,
                                                        sample_infos=sample_infos_file,
                                                        )
        pl = CMMDBPipeline(jobs_setup_file)
        self.assertEqual(pl.samples_list,
                         ["Alb-31", "Br-466", "Br-432", "Al-161", "Br-504", "Al-65"],
                         "CMMDBPipeline cannot correctly read 'samples list' from jobs setup file")

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_cal_mut_stat_offline_1(self):
        """ test basic offline version (w/o slurm) of cal_mut_stat (one chrom)"""

        self.individual_debug = True
        self.init_test(self.current_func_name)
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        exp_result = join_path(self.data_dir,
                               "exp_stat")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        project_code=None,
                                                        db_region=None,
                                                        )
        pl = CMMDBPipeline(jobs_setup_file)
        pl.cal_mut_stat()
        self.assertTrue(filecmp.cmp(pl.out_stat_file,
                                    exp_result),
                        "cal_mut_stat doesn't function correctly")


    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_cal_mut_stat_offline_2(self):
        """ test a little advance offline version (w/o slurm) of cal_mut_stat (all chroms) """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        vcf_tabix_file = join_path(self.data_dir,
                                   "multi_chrs.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        project_code=None,
                                                        db_region=None,
                                                        )
        pl = CMMDBPipeline(jobs_setup_file)
        pl.cal_mut_stat()

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_cal_mut_stat_offline_3(self):
        """ to offline version (w/o slurm) if it can handle multi allelic stat correctly """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        exp_result = join_path(self.data_dir,
                               "exp_stat")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        project_code=None,
                                                        db_region=None,
                                                        )
        pl = CMMDBPipeline(jobs_setup_file)
        pl.cal_mut_stat()
        self.assertTrue(filecmp.cmp(pl.out_stat_file,
                                    exp_result),
                        "cal_mut_stat doesn't function correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_cal_mut_stat_offline_4(self):
        """ to offline version (w/o slurm) if it can calculate stat of subpopulation """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        exp_result = join_path(self.data_dir,
                               "exp_stat")
        sample_infos = "Br-176"
        sample_infos += ",Br-120"
        sample_infos += ",Alb-31"
        sample_infos += ",Br-331"
        sample_infos += ",Al-77"
        sample_infos += ",Br-472"
        sample_infos += ",Br-706"
        sample_infos += ",Br-366"
        sample_infos += ",Al-140"
        sample_infos += ",Al-169"
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        project_code=None,
                                                        db_region=None,
                                                        sample_infos=sample_infos,
                                                        )
        pl = CMMDBPipeline(jobs_setup_file)
        pl.cal_mut_stat()
        self.assertTrue(filecmp.cmp(pl.out_stat_file,
                                    exp_result),
                        "cal_mut_stat doesn't function correctly")

    @unittest.skipUnless(SLURM_CMMDB_TEST, "taking too much UPPMAX cpu-core hours")
    def test_cal_mut_stat_slurm_1(self):
        """ test basic slurm version of cal_mut_stat (indicate vcf region)"""

        self.individual_debug = True
        self.init_test(self.current_func_name)
        vcf_tabix_file = join_path(self.data_dir,
                                   "multi_chrs.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(db_alloc_time="5:00:00",
                                                        dataset_name="test_cal_stat_1",
                                                        vcf_tabix_file=vcf_tabix_file,
                                                        project_code=FAST_PROJECT_CODE,
                                                        db_region="X",
                                                        )
        pl = CMMDBPipeline(jobs_setup_file)
        pl.cal_mut_stat()

    @unittest.skipUnless(SLURM_CMMDB_TEST, "taking too much UPPMAX cpu-core hours")
    def test_cal_mut_stat_slurm_2(self):
        """ test a little advance slurm version of cal_mut_stat (all chroms)"""

        self.individual_debug = True
        self.init_test(self.current_func_name)
        vcf_tabix_file = join_path(self.data_dir,
                                   "multi_chrs.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name="test_cal_stat_2",
                                                        vcf_tabix_file=vcf_tabix_file,
                                                        project_code=FAST_PROJECT_CODE,
                                                        db_region=None,
                                                        )
        pl = CMMDBPipeline(jobs_setup_file)
        pl.cal_mut_stat()

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_table_annovar_offline_1(self):
        """ test offline version (w/o slurm) of table_annovar """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        project_code=None)
        pl = CMMDBPipeline(jobs_setup_file)
        pl.table_annovar()

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_table_annovar_offline_2(self):
        """ test offline version (w/o slurm) of table_annovar together with custom cal_stat """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        annovar_db_names = DFLT_ANNOVAR_TEST_DB_NAMES
        annovar_db_names += ",test_pyCMM"
        annovar_db_ops = DFLT_ANNOVAR_TEST_DB_OPS
        annovar_db_ops += ",f"
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annovar_db_names=annovar_db_names,
                                                        annovar_db_ops=annovar_db_ops,
                                                        project_code=None)
        pl = CMMDBPipeline(jobs_setup_file)
        pl.table_annovar()

    @unittest.skipUnless(SLURM_CMMDB_TEST, "taking too much UPPMAX cpu-core hours")
    def test_table_annovar_slurm(self):
        """ test slurm version of table_annovar """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        project_code=FAST_PROJECT_CODE)
        pl = CMMDBPipeline(jobs_setup_file)
        pl.table_annovar()

    def tearDown(self):
        self.remove_working_dir()
