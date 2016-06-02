import unittest
import filecmp
import datetime
from os.path import exists as path_exists
from os.path import join as join_path
from os.path import dirname
from os.path import isdir
from pycmm.template import SafeTester
from pycmm.settings import FULL_SYSTEM_TEST
#from pycmm.settings import FAST_PROJECT_CODE
from pycmm.settings import SLOW_PROJECT_CODE
from pycmm.utils import count_lines
from pycmm.flow.cmmdb import CMMDBPipeline
from pycmm.flow.cmmdb import create_jobs_setup_file
from pycmm.cmmlib.test.test_annovarlib import DFLT_ANV_TEST_DB_DIR
from pycmm.cmmlib.test.test_annovarlib import DFLT_ANV_TEST_DB_NAMES
from pycmm.cmmlib.test.test_annovarlib import DFLT_ANV_TEST_DB_OPS
from pycmm.cmmlib.familylib import NO_FAMILY

CMMDB_TEST=False


class TestCMMDBPipeline(SafeTester):

    def __init__(self, methodName):
        super(TestCMMDBPipeline, self).__init__(methodName=methodName,
                                                test_module_name=__name__,
                                                )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self,
                                 vcf_tabix_file,
                                 project_name=None,
                                 project_code=None,
                                 job_alloc_time=None,
                                 db_region="18",
                                 sample_info=None,
                                 annovar_human_db_dir=DFLT_ANV_TEST_DB_DIR,
                                 annovar_buildver="hg19",
                                 annovar_db_names=DFLT_ANV_TEST_DB_NAMES,
                                 annovar_db_ops=DFLT_ANV_TEST_DB_OPS,
                                 annovar_nastring=".",
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        if project_name is None:
            project_name = self.test_function
        create_jobs_setup_file(project_name=project_name,
                               project_out_dir=self.working_dir,
                               vcf_tabix_file=vcf_tabix_file,
                               db_region=db_region,
                               sample_info=sample_info,
                               project_code=project_code,
                               job_alloc_time=job_alloc_time,
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
        pl = CMMDBPipeline(jobs_setup_file=jobs_setup_file)
        exp_project_name = self.test_function
        self.assertEqual(pl.project_name,
                         exp_project_name,
                         "CMMDBPipeline cannot correctly read meta info 'dataset name' from jobs setup file")
        self.assertEqual(pl.mutstat_params.input_vcf_tabix,
                         dummy_vcf_tabix_file,
                         "CMMDBPipeline cannot correctly read meta info 'input vcf tabix' from jobs setup file")
        self.assertEqual(pl.mutstat_params.db_region,
                         "18",
                         "CMMDBPipeline cannot correctly read meta info 'db region' from jobs setup file")
        self.assertEqual(pl.samples_dict,
                         None,
                         "CMMDBPipeline cannot correctly read meta info 'samples' from jobs setup file")
        exp_jobs_report_file = join_path(pl.project_out_dir,
                                         exp_project_name+"_rpt.txt")
        self.assertEqual(pl.jobs_report_file,
                         exp_jobs_report_file,
                         "CMMDBPipeline cannot correctly read meta info 'jobs report file' from jobs setup file")
        self.assertEqual(pl.project_out_dir,
                         self.working_dir,
                         "CMMDBPipeline cannot correctly read meta info 'output dir' from jobs setup file")
        self.assertEqual(pl.annovar_params.input_vcf_tabix,
                         dummy_vcf_tabix_file,
                         "CMMDBPipeline cannot corretly read annovar configs 'annovar input file' from jobs setup file")
        self.assertEqual(pl.annovar_params.db_dir,
                         DFLT_ANV_TEST_DB_DIR,
                         "CMMDBPipeline cannot corretly read annovar configs 'annovar db folder' from jobs setup file")
        self.assertEqual(pl.annovar_params.buildver,
                         "hg19",
                         "CMMDBPipeline cannot corretly read annovar configs 'annovar buildver' from jobs setup file")
        self.assertEqual(pl.annovar_params.nastring,
                         ".",
                         "CMMDBPipeline cannot corretly read annovar configs 'annovar nastring' from jobs setup file")
        self.assertEqual(pl.annovar_params.protocols,
                         DFLT_ANV_TEST_DB_NAMES,
                         "CMMDBPipeline cannot corretly read annovar configs 'annovar db names' from jobs setup file")
        self.assertEqual(pl.annovar_params.operations,
                         DFLT_ANV_TEST_DB_OPS,
                         "CMMDBPipeline cannot corretly read annovar configs 'annovar db ops' from jobs setup file")

    def test_load_jobs_info_2(self):
        """ test if modified job configurations are loaded correctly """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        dummy_vcf_tabix_file = "/path/to/vcf_tabix_file"
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=dummy_vcf_tabix_file,
                                                        project_code=SLOW_PROJECT_CODE,
                                                        job_alloc_time="120:00:00",
                                                        db_region=None,
                                                        sample_info="Co-441,Co-666",
                                                        )
        pl = CMMDBPipeline(jobs_setup_file=jobs_setup_file)
        exp_project_name = self.test_function
        self.assertEqual(pl.project_code,
                         SLOW_PROJECT_CODE,
                         "CMMDBPipeline cannot correctly read meta info 'project code' from jobs setup file")
        self.assertEqual(pl.job_alloc_time,
                         "120:00:00",
                         "CMMDBPipeline cannot correctly read meta info 'job allocation time' from jobs setup file")
        self.assertEqual(pl.mutstat_params.db_region,
                         None,
                         "CMMDBPipeline cannot correctly read meta info 'vcf region' from jobs setup file")
        self.assertEqual(pl.samples_id,
                         ["Co-441","Co-666"],
                         "CMMDBPipeline cannot correctly read meta info 'samples' from jobs setup file")
        self.assertEqual(len(pl.families_info),
                         1,
                         "CMMDBPipeline cannot correctly read meta info 'families list' from jobs setup file")
        self.assertTrue(NO_FAMILY in pl.families_info,
                        "CMMDBPipeline cannot correctly read meta info 'families list' from jobs setup file")

    def test_load_jobs_info_3(self):
        """ test if families structures are loaded correctly """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        dummy_vcf_tabix_file = "/path/to/vcf_tabix_file"
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=dummy_vcf_tabix_file,
                                                        sample_info="242:Co-441:Co-666:Co-771,6067:Br-397:884-99D",
                                                        )
        pl = CMMDBPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(pl.samples_id,
                         ["Co-441","Co-666","Co-771","Br-397","884-99D"],
                         "CMMDBPipeline cannot correctly read 'samples' from jobs setup file")
        self.assertEqual(len(pl.families_info),
                         2,
                         "CMMDBPipeline cannot correctly read meta info 'families list' from jobs setup file")
        fam_info_1 = pl.families_info['242']
        fam_info_2 = pl.families_info['6067']
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
        sample_info_file = join_path(self.data_dir,
                                      "sample_info")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=dummy_vcf_tabix_file,
                                                        sample_info=sample_info_file,
                                                        )
        pl = CMMDBPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(len(pl.samples_id),
                         5,
                         "CMMDBPipeline cannot correctly read 'samples' from jobs setup file")
        self.assertTrue("Co-35" in pl.samples_id,
                        "CMMDBPipeline cannot correctly read 'samples' from jobs setup file")
        self.assertTrue("Co-37" in pl.samples_id,
                        "CMMDBPipeline cannot correctly read 'samples' from jobs setup file")
        self.assertTrue("Co-135" in pl.samples_id,
                        "CMMDBPipeline cannot correctly read 'samples' from jobs setup file")
        self.assertTrue("Co-131" in pl.samples_id,
                        "CMMDBPipeline cannot correctly read 'samples' from jobs setup file")
        self.assertTrue("1322-11D" in pl.samples_id,
                        "CMMDBPipeline cannot correctly read 'samples' from jobs setup file")

    def test_load_jobs_info_5(self):
        """ test if can load sample infos in family structured format as a file """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        dummy_vcf_tabix_file = "/path/to/vcf_tabix_file"
        sample_info_file = join_path(self.data_dir,
                                      "sample_info")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=dummy_vcf_tabix_file,
                                                        sample_info=sample_info_file,
                                                        )
        pl = CMMDBPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(pl.samples_id,
                         ["Alb-31", "Br-466", "Br-432", "Al-161", "Br-504", "Al-65"],
                         "CMMDBPipeline cannot correctly read 'samples' from jobs setup file")

    @unittest.skipUnless(FULL_SYSTEM_TEST or CMMDB_TEST, "taking too long time to test")
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
        pl = CMMDBPipeline(jobs_setup_file=jobs_setup_file)
        pl.cal_mut_stat()
        self.assertTrue(filecmp.cmp(pl.out_stat_files,
                                    exp_result),
                        "cal_mut_stat doesn't function correctly")


    @unittest.skipUnless(FULL_SYSTEM_TEST or CMMDB_TEST, "taking too long time to test")
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
        pl = CMMDBPipeline(jobs_setup_file=jobs_setup_file)
        pl.cal_mut_stat()

    @unittest.skipUnless(FULL_SYSTEM_TEST or CMMDB_TEST, "taking too long time to test")
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
        pl = CMMDBPipeline(jobs_setup_file=jobs_setup_file)
        pl.cal_mut_stat()
        self.assertTrue(filecmp.cmp(pl.out_stat_files,
                                    exp_result),
                        "cal_mut_stat doesn't function correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST or CMMDB_TEST, "taking too long time to test")
    def test_cal_mut_stat_offline_4(self):
        """ to offline version (w/o slurm) if it can calculate stat of subpopulation """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        exp_result = join_path(self.data_dir,
                               "exp_stat")
        sample_info = "Br-176"
        sample_info += ",Br-120"
        sample_info += ",Alb-31"
        sample_info += ",Br-331"
        sample_info += ",Al-77"
        sample_info += ",Br-472"
        sample_info += ",Br-706"
        sample_info += ",Br-366"
        sample_info += ",Al-140"
        sample_info += ",Al-169"
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        project_code=None,
                                                        db_region=None,
                                                        sample_info=sample_info,
                                                        )
        pl = CMMDBPipeline(jobs_setup_file=jobs_setup_file)
        pl.cal_mut_stat()
        self.assertTrue(filecmp.cmp(pl.out_stat_files,
                                    exp_result),
                        "cal_mut_stat doesn't function correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST or CMMDB_TEST, "taking too long time to test")
    def test_table_annovar_offline_1(self):
        """ test offline version (w/o slurm) of table_annovar """

        self.init_test(self.current_func_name)
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        project_code=None)
        pl = CMMDBPipeline(jobs_setup_file=jobs_setup_file)
        pl.table_annovar()
        out_annotated_vcf = join_path(self.working_dir,
                                      self.current_func_name+"_annotated.vcf")
        self.assertEqual(count_lines(out_annotated_vcf),
                         151,
                         "table annovar result is incorrect ")

    @unittest.skipUnless(FULL_SYSTEM_TEST or CMMDB_TEST, "taking too long time to test")
    def test_table_annovar_offline_2(self):
        """ test offline version (w/o slurm) of table_annovar together with custom cal_stat """

        self.init_test(self.current_func_name)
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        annovar_db_names = DFLT_ANV_TEST_DB_NAMES
        annovar_db_names += ",test_pyCMM"
        annovar_db_ops = DFLT_ANV_TEST_DB_OPS
        annovar_db_ops += ",f"
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annovar_db_names=annovar_db_names,
                                                        annovar_db_ops=annovar_db_ops,
                                                        project_code=None)
        pl = CMMDBPipeline(jobs_setup_file=jobs_setup_file)
        pl.table_annovar()
        out_annotated_vcf = join_path(self.working_dir,
                                      self.current_func_name+"_annotated.vcf")
        self.assertEqual(count_lines(out_annotated_vcf),
                         215,
                         "table annovar result is incorrect ")
