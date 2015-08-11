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
from pycmm.settings import FAST_PROJECT_CODE
from pycmm.settings import SLOW_PROJECT_CODE
from pycmm.settings import DFLT_MUTREP_FREQ_RATIOS
from pycmm.settings import DFLT_MUTREP_ANNO_COLS
from pycmm.flow.mutrep import MutRepPipeline
from pycmm.flow.cmmdb import create_jobs_setup_file

DFLT_TEST_MUTREP_COLS = []
DFLT_TEST_MUTREP_COLS.append("Func.refGene")
DFLT_TEST_MUTREP_COLS.append("Gene.refGene")
DFLT_TEST_MUTREP_COLS.append("GeneDetail.refGene")
DFLT_TEST_MUTREP_COLS.append("ExonicFunc.refGene")
DFLT_TEST_MUTREP_COLS.append("cytoBand")
DFLT_TEST_MUTREP_COLS.append("1000g2014oct_all")

DFLT_TEST_REPORT_REGIONS = "18:12512255-14542551"
DFLT_TEST_FREQ_RATIOS = DFLT_MUTREP_FREQ_RATIOS
DFLT_TEST_FREQ_RATIOS += ",ExAC_ALL:0.1"

class TestMutRepPipeline(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            test_module_name=__name__,
                            )

    def setUp(self):
        self.module_name = 'mutrep'

    def __create_jobs_setup_file(self,
                                 vcf_tabix_file,
                                 dataset_name=None,
                                 project_code=SLOW_PROJECT_CODE,
                                 vcf_region="18",
                                 samples_info=None,
                                 anno_cols=DFLT_TEST_MUTREP_COLS,
                                 report_regions=DFLT_TEST_REPORT_REGIONS,
                                 call_info="NO",
                                 frequency_ratios=DFLT_TEST_FREQ_RATIOS,
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        if dataset_name is None:
            dataset_name = self.test_function
        create_jobs_setup_file(dataset_name=dataset_name,
                               project_out_dir=self.working_dir,
                               vcf_tabix_file=vcf_tabix_file,
                               vcf_region=vcf_region,
                               samples_info=samples_info,
                               project_code=project_code,
                               anno_cols=",".join(anno_cols),
                               report_regions=report_regions,
                               call_info=call_info,
                               frequency_ratios=frequency_ratios,
                               out_jobs_setup_file=jobs_setup_file,
                               )
        return jobs_setup_file

    def test_load_jobs_info_1(self):
        """ test if layout configurations are loaded correctly """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        dummy_vcf_tabix_file = "/path/to/vcf_tabix_file"
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=dummy_vcf_tabix_file,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        self.assertEqual(pl.report_layout.anno_cols[2],
                         "GeneDetail.refGene",
                         "MutRepPipeline cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(pl.report_layout.anno_cols[4],
                         "cytoBand",
                         "MutRepPipeline cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(pl.report_layout.report_regions[0].chrom,
                         "18",
                         "MutRepPipeline cannot correctly read report layout info 'report regions' from jobs setup file")
        self.assertEqual(pl.report_layout.call_info,
                         False,
                         "MutRepPipeline cannot correctly read report layout info 'call info' from jobs setup file")
        mylogger.debug(pl.report_layout.freq_ratios)
        self.assertEqual(pl.report_layout.freq_ratios["esp6500siv2_all"],
                         0.2,
                         "MutRepPipeline cannot correctly read report layout info 'frequency ratios' from jobs setup file")
        self.assertEqual(pl.report_layout.freq_ratios["ExAC_ALL"],
                         0.1,
                         "MutRepPipeline cannot correctly read report layout info 'frequency ratios' from jobs setup file")

    def test_load_jobs_info_2(self):
        """ test if non-default layout configurations are loaded correctly """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        dummy_vcf_tabix_file = "/path/to/vcf_tabix_file"
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=dummy_vcf_tabix_file,
                                                        report_regions="6:78161823-78164117,"+DFLT_TEST_REPORT_REGIONS+",22",
                                                        call_info="YES",
                                                        frequency_ratios="ExAC:0.5",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        self.assertEqual(pl.report_layout.report_regions[1].end_pos,
                         "14542551",
                         "MutRepPipeline cannot correctly read report layout info 'report regions' from jobs setup file")
        self.assertEqual(pl.report_layout.report_regions[2].chrom,
                         "22",
                         "MutRepPipeline cannot correctly read report layout info 'report regions' from jobs setup file")
        self.assertEqual(pl.report_layout.report_regions[2].start_pos,
                         None,
                         "MutRepPipeline cannot correctly read report layout info 'report regions' from jobs setup file")
        self.assertEqual(pl.report_layout.call_info,
                         True,
                         "MutRepPipeline cannot correctly read report layout info 'call info' from jobs setup file")
        self.assertEqual(pl.report_layout.freq_ratios["ExAC"],
                         0.5,
                         "MutRepPipeline cannot correctly read report layout info 'frequency ratios' from jobs setup file")

    def test_load_jobs_info_3(self):
        """ test if non-default (None) layout configurations are loaded correctly """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        dummy_vcf_tabix_file = "/path/to/vcf_tabix_file"
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=dummy_vcf_tabix_file,
                                                        report_regions=None,
                                                        frequency_ratios=None,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        self.assertEqual(pl.report_layout.report_regions,
                         None,
                         "MutRepPipeline cannot correctly read report layout info 'report regions' from jobs setup file")
        self.assertEqual(pl.report_layout.freq_ratios,
                         None,
                         "MutRepPipeline cannot correctly read report layout info 'frequency ratios' from jobs setup file")

#    @unittest.skip("Disable for temporary")
    def test_summary_report_2(self):
        """ test if summary reprot with custom configuration can be correctly generated """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "chr18.vcf.gz")
        rpt_out_file = join_path(self.working_dir,
                                 self.current_func_name + ".xlsx")
        anno_cols=DFLT_TEST_MUTREP_COLS
        anno_cols.remove("1000g2014oct_all")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        anno_cols=anno_cols,
                                                        call_info="YES",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_report(out_file=rpt_out_file)

#    @unittest.skip("Disable for temporary")
    def test_summary_report_4(self):
        """ test if summary report with default configuration can be correctly generated """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        rpt_out_file = join_path(self.working_dir,
                                 self.current_func_name + ".xlsx")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        anno_cols=DFLT_MUTREP_ANNO_COLS,
                                                        report_regions="6",
                                                        call_info="YES",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_report(out_file=rpt_out_file)

    def tearDown(self):
        self.remove_working_dir()
