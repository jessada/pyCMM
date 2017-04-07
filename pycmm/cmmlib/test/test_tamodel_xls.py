import unittest
from os.path import join as join_path
from os.path import dirname
from pycmm.template import SafeTester
from pycmm.settings import ALL_MUTREP_ANNO_COLS
from pycmm.settings import PRIMARY_MAF_VAR
from pycmm.settings import EXAC_ALL_COL_NAME
from pycmm.settings import PATHOGENIC_COUNT_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_SYN_Z_COL_NAME
from pycmm.settings import INTERVAR_CLASS_COL_NAME
from pycmm.settings import INTERVAR_EVIDENCE_COL_NAME
from pycmm.settings import FULL_SYSTEM_TEST
from pycmm.settings import AXEQ_CHR9_COLS_TAG
from pycmm.settings import AXEQ_CHR3_6_14_18_COLS_TAG
from pycmm.settings import AXEQ_CHR5_19_COLS_TAG
from pycmm.settings import EXAC_CONSTRAINT_COLS_TAG
from pycmm.settings import LJB_SCORE_COLS_TAG
from pycmm.cmmlib.xlslib import XlsUtils
from pycmm.cmmlib.intervarlib import INTERVAR_CLASS_BENIGN
from pycmm.cmmlib.intervarlib import INTERVAR_CLASS_LIKELY_BENIGN
from pycmm.cmmlib.intervarlib import INTERVAR_CLASS_UNCERTAIN_SIGNIFICANCE
from pycmm.cmmlib.intervarlib import INTERVAR_CLASS_LIKELY_PATHOGENIC
from pycmm.cmmlib.intervarlib import INTERVAR_CLASS_PATHOGENIC
from pycmm.flow.mutrep import MutRepPipeline
from pycmm.flow.mutrep import create_jobs_setup_file
from pycmm.flow.mutrep import JOBS_SETUP_RPT_FILTER_RARE
from pycmm.flow.mutrep import JOBS_SETUP_RPT_FILTER_PASS_VQSR
from pycmm.flow.mutrep import JOBS_SETUP_RPT_FILTER_NON_INTRONIC
from pycmm.flow.mutrep import JOBS_SETUP_RPT_FILTER_NON_INTERGENIC
from pycmm.flow.mutrep import JOBS_SETUP_RPT_FILTER_NON_UPSTREAM
from pycmm.flow.mutrep import JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM
from pycmm.flow.mutrep import JOBS_SETUP_RPT_FILTER_NON_UTR
from pycmm.flow.mutrep import JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS
from pycmm.flow.mutrep import JOBS_SETUP_RPT_FILTER_HAS_MUTATION
from pycmm.flow.mutrep import JOBS_SETUP_RPT_FILTER_HAS_SHARED
from pycmm.flow.test.test_mutrep import DFLT_TEST_MUTREP_COLS
from pycmm.flow.test.test_mutrep import DFLT_TEST_ANNO_EXCL_TAGS


class TestTAVcfCallXls(SafeTester):

    def __init__(self, methodName):
        super(TestTAVcfCallXls, self).__init__(methodName=methodName,
                                               test_module_name=__name__,
                                               )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self,
                                 project_name=None,
                                 sample_info=None,
                                 anno_cols=DFLT_TEST_MUTREP_COLS,
                                 anno_excl_tags=None,
                                 annotated_vcf_tabix=None,
                                 summary_families_sheet=False,
                                 frequency_ratios=None,
                                 rows_filter_actions=None,
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        if project_name is None:
            project_name = self.test_function
        create_jobs_setup_file(project_name=project_name,
                               project_out_dir=self.working_dir,
                               sample_info=sample_info,
                               anno_cols=",".join(anno_cols),
                               anno_excl_tags=anno_excl_tags,
                               annotated_vcf_tabix=annotated_vcf_tabix,
                               summary_families_sheet=summary_families_sheet,
                               frequency_ratios=frequency_ratios,
                               rows_filter_actions=rows_filter_actions,
                               out_jobs_setup_file=jobs_setup_file,
                               )
        return jobs_setup_file

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_parse_mutated_xls_3_2(self):
        """
        - test if shared mutations in an xls sheet with a family with 2 members 
        can be detected
        - this test is corresponding to 'test_parse_mutated_3_2' in 'test_tamodel'
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        project_name = self.test_function
        sample_info = []
        sample_info.append("24:Co-166:Co-213")
        frequency_ratios = PRIMARY_MAF_VAR + ":0.2"
        frequency_ratios += "," + ILL_BR_PF_COL_NAME + ":0.3"
        frequency_ratios += "," + EXAC_ALL_COL_NAME + ":0.3"
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=DFLT_TEST_MUTREP_COLS,
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        frequency_ratios=frequency_ratios,
                                                        sample_info=",".join(sample_info),
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_family_report('24', pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_fam24.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         1,
                         "shared mutations cannot be correctly determined")

class TestTAVcfRecordXls(SafeTester):

    def __init__(self, methodName):
        super(TestTAVcfRecordXls, self).__init__(methodName=methodName,
                                                 test_module_name=__name__,
                                                 )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self,
                                 project_name=None,
                                 sample_info=None,
                                 anno_cols=DFLT_TEST_MUTREP_COLS,
                                 anno_excl_tags=None,
                                 annotated_vcf_tabix=None,
                                 summary_families_sheet=False,
                                 frequency_ratios=None,
                                 rows_filter_actions=None,
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        if project_name is None:
            project_name = self.test_function
        create_jobs_setup_file(project_name=project_name,
                               project_out_dir=self.working_dir,
                               sample_info=sample_info,
                               anno_cols=",".join(anno_cols),
                               anno_excl_tags=anno_excl_tags,
                               annotated_vcf_tabix=annotated_vcf_tabix,
                               summary_families_sheet=summary_families_sheet,
                               frequency_ratios=frequency_ratios,
                               rows_filter_actions=rows_filter_actions,
                               out_jobs_setup_file=jobs_setup_file,
                               )
        return jobs_setup_file

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_has_mutation_xls_4(self):
        """
        - test report that show only mutation
        - this test doesn't corresponding to any other tests
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        project_name = self.test_function+'_has_no_mutation'
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        sample_info="24:Co-166:Co-213:Co-648,8:Co-37,275:Co-618,478:Co-1274",
                                                        anno_cols=DFLT_TEST_MUTREP_COLS,
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        summary_families_sheet=True,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=1),
                         34,
                         "shared mutation cannot be correctly determined")
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_HAS_MUTATION
        project_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        sample_info="24:Co-166:Co-213:Co-648,8:Co-37,275:Co-618,478:Co-1274",
                                                        anno_cols=DFLT_TEST_MUTREP_COLS,
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        summary_families_sheet=True,
                                                        rows_filter_actions=rows_filter_actions,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=1),
                         13,
                         "shared mutation cannot be correctly determined")

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_has_shared_xls_4(self):
        """
        - test if shared mutations in an xls sheet with a family with 2 members 
        can be detected
        - this test doesn't corresponding to any other tests
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        project_name = self.test_function
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_RARE
        rows_filter_actions += "," + JOBS_SETUP_RPT_FILTER_NON_INTRONIC
        rows_filter_actions += "," + JOBS_SETUP_RPT_FILTER_NON_INTERGENIC
        rows_filter_actions += "," + JOBS_SETUP_RPT_FILTER_NON_UPSTREAM
        rows_filter_actions += "," + JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM
        rows_filter_actions += "," + JOBS_SETUP_RPT_FILTER_NON_UTR
        rows_filter_actions += "," + JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS
        rows_filter_actions += "," + JOBS_SETUP_RPT_FILTER_HAS_MUTATION
        sample_info = []
        sample_info.append("24:Co-166:Co-213")
        frequency_ratios = PRIMARY_MAF_VAR + ":0.2"
        frequency_ratios += "," + ILL_BR_PF_COL_NAME + ":0.3"
        frequency_ratios += "," + EXAC_ALL_COL_NAME + ":0.3"
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=DFLT_TEST_MUTREP_COLS,
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        rows_filter_actions=rows_filter_actions,
                                                        frequency_ratios=frequency_ratios,
                                                        sample_info=",".join(sample_info),
                                                        summary_families_sheet=True,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_family_report('24', pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_fam24.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         5,
                         "shared mutations cannot be correctly determined")


    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_is_rare_xls_3(self):
        """
        - test filter rare feature with xls records
        - this test doesn't corresponding to any other tests
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        project_name = self.test_function + '_with_common'
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        anno_cols=ALL_MUTREP_ANNO_COLS,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         88,
                         "rare mutations cannot be correctly determined")
        project_name = self.test_function
        frequency_ratios = PRIMARY_MAF_VAR + ":0.2"
        frequency_ratios += "," + ILL_BR_PF_COL_NAME + ":0.3"
        frequency_ratios += "," + EXAC_ALL_COL_NAME + ":0.3"
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_RARE
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=ALL_MUTREP_ANNO_COLS,
                                                        frequency_ratios=frequency_ratios,
                                                        rows_filter_actions=rows_filter_actions,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         54,
                         "rare mutations cannot be correctly determined")

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_is_pass_vqsr_xls_2(self):
        """
        - test filter SNPs that pass QC (VQSR)
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        project_name = self.test_function + '_all'
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        anno_cols=ALL_MUTREP_ANNO_COLS,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        )
        pl = MutRepPipeline(jobs_setup_file=jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         7,
                         "number of mutations that pass VQSR cannot be correctly determined")
        project_name = self.test_function + '_pass_vqsr'
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_PASS_VQSR
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=ALL_MUTREP_ANNO_COLS,
                                                        rows_filter_actions=rows_filter_actions,
                                                        )
        pl = MutRepPipeline(jobs_setup_file=jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         5,
                         "number of mutations that pass VQSR cannot be correctly determined")

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_is_intergenic_xls_3(self):
        """
        - test filter non-intergenic feature with xls record
        - this test doesn't corresponding to any other tests
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        project_name = self.test_function + '_with_intergenic'
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         22,
                         "intergenic mutations cannot be correctly determined")
        project_name = self.test_function
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_NON_INTERGENIC
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        rows_filter_actions=rows_filter_actions,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         11,
                         "intergenic mutations cannot be correctly determined")

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_is_intronic_xls_3(self):
        """
        - test filter non-intronic feature with xls records
        - this test doesn't corresponding to any other tests
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        project_name = self.test_function + '_with_intronic'
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=DFLT_TEST_MUTREP_COLS,
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         26,
                         "intronic mutations cannot be correctly determined")
        project_name = self.test_function
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_NON_INTRONIC
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=DFLT_TEST_MUTREP_COLS,
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        rows_filter_actions=rows_filter_actions,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         9,
                         "intronic mutations cannot be correctly determined")

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_is_upstream_xls_3(self):
        """
        - test filter non-upstream feature with xls records
        - this test doesn't corresponding to any other tests
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        project_name = self.test_function + '_with_upstream'
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=DFLT_TEST_MUTREP_COLS,
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         32,
                         "upstream mutations cannot be correctly determined")
        project_name = self.test_function
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_NON_UPSTREAM
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=DFLT_TEST_MUTREP_COLS,
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        rows_filter_actions=rows_filter_actions,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         28,
                         "upstream mutations cannot be correctly determined")

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_is_downstream_xls_3(self):
        """
        - test filter non-downstream feature with xls records
        - this test doesn't corresponding to any other tests
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        project_name = self.test_function + '_with_downstream'
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=DFLT_TEST_MUTREP_COLS,
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         32,
                         "downstream mutations cannot be correctly determined")
        project_name = self.test_function
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=DFLT_TEST_MUTREP_COLS,
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        rows_filter_actions=rows_filter_actions,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         31,
                         "downstream mutations cannot be correctly determined")

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_is_utr_xls_3(self):
        """
        - test filter non-UTR feature with xls records
        - this test doesn't corresponding to any other tests
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        project_name = self.test_function + '_with_utr'
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=DFLT_TEST_MUTREP_COLS,
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         32,
                         "UTR mutations cannot be correctly determined")
        project_name = self.test_function
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_NON_UTR
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=DFLT_TEST_MUTREP_COLS,
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        rows_filter_actions=rows_filter_actions,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         20,
                         "UTR mutations cannot be correctly determined")

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_is_synonymous_xls_3(self):
        """
        - test filter non-synonymous feature with xls records
        - this test doesn't corresponding to any other tests
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        project_name = self.test_function + '_with_synonymous'
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=DFLT_TEST_MUTREP_COLS,
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         32,
                         "synonymous mutations cannot be correctly determined")
        project_name = self.test_function
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=DFLT_TEST_MUTREP_COLS,
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        rows_filter_actions=rows_filter_actions,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         28,
                         "synonymous mutations cannot be correctly determined")

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_pathogenic_count_xls_2(self):
        """
        - test if number of pathogenic count can be correctly displayed
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        project_name = self.test_function
        custom_excl_tags = DFLT_TEST_ANNO_EXCL_TAGS
        custom_excl_tags += "," + AXEQ_CHR3_6_14_18_COLS_TAG
        custom_excl_tags += "," + AXEQ_CHR5_19_COLS_TAG
        custom_excl_tags += "," + LJB_SCORE_COLS_TAG
        custom_excl_tags += "," + EXAC_CONSTRAINT_COLS_TAG
        custom_excl_tags += "," + LJB_SCORE_COLS_TAG
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_NON_INTERGENIC
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_INTRONIC
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_UPSTREAM
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_UTR
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=ALL_MUTREP_ANNO_COLS,
                                                        anno_excl_tags=custom_excl_tags,
                                                        rows_filter_actions=rows_filter_actions,
                                                        )
        pl = MutRepPipeline(jobs_setup_file=jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        pc_col_idx = xu.get_col_idx(PATHOGENIC_COUNT_COL_NAME)
        self.assertEqual(xu.get_cell_value(4, pc_col_idx),
                         1,
                         "Incorect number of harmful pathogenic predictions"
                         )
        self.assertEqual(xu.get_cell_value(5, pc_col_idx),
                         0,
                         "Incorect number of harmful pathogenic predictions"
                         )
        self.assertEqual(xu.get_cell_value(9, pc_col_idx),
                         5,
                         "Incorect number of harmful pathogenic predictions"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_parse_exac_constraint_xls_3(self):
        """
        - test if exac constraints can be correctly parsed and displayed
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        project_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=ALL_MUTREP_ANNO_COLS,
                                                        )
        pl = MutRepPipeline(jobs_setup_file=jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        syn_z_col_idx = xu.get_col_idx(pl.correct_header(EXAC03_CONSTRAINT_SYN_Z_COL_NAME))
        self.assertEqual(xu.get_cell_value(2, syn_z_col_idx),
                         None,
                         "Incorect number of ExAC constraint value"
                         )
        self.assertEqual(round(xu.get_cell_value(4, syn_z_col_idx), 2),
                         0.53,
                         "Incorect number of ExAC constraint value"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_intervar_1(self):
        """
        - test if intervar classification and evidence can be parsed
        """

        self.init_test(self.current_func_name)
        anno_cols = list(DFLT_TEST_MUTREP_COLS)
        anno_cols.append(INTERVAR_CLASS_COL_NAME)
        anno_cols.append(INTERVAR_EVIDENCE_COL_NAME)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        project_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=anno_cols,
                                                        )
        pl = MutRepPipeline(jobs_setup_file=jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             project_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        intervar_class_col_idx = xu.get_col_idx(pl.correct_header(INTERVAR_CLASS_COL_NAME))
        self.assertEqual(xu.get_cell_value(2, intervar_class_col_idx),
                         INTERVAR_CLASS_BENIGN,
                         "Incorect intervar value"
                         )
        self.assertEqual(xu.get_cell_value(3, intervar_class_col_idx),
                         INTERVAR_CLASS_LIKELY_BENIGN,
                         "Incorect intervar value"
                         )
        self.assertEqual(xu.get_cell_value(4, intervar_class_col_idx),
                         INTERVAR_CLASS_UNCERTAIN_SIGNIFICANCE,
                         "Incorect intervar value"
                         )
        self.assertEqual(xu.get_cell_value(5, intervar_class_col_idx),
                         INTERVAR_CLASS_BENIGN,
                         "Incorect intervar value"
                         )
        self.assertEqual(xu.get_cell_value(6, intervar_class_col_idx),
                         INTERVAR_CLASS_PATHOGENIC,
                         "Incorect intervar value"
                         )
        self.assertEqual(xu.get_cell_value(7, intervar_class_col_idx),
                         INTERVAR_CLASS_LIKELY_PATHOGENIC,
                         "Incorect intervar value"
                         )
        intervar_evidence_col_idx = xu.get_col_idx(pl.correct_header(INTERVAR_EVIDENCE_COL_NAME))
        self.assertEqual(xu.get_cell_value(2, intervar_evidence_col_idx),
                         "BA1, BS1, BP4, BP7",
                         "Incorect intervar value"
                         )
        self.assertEqual(xu.get_cell_value(3, intervar_evidence_col_idx),
                         "PM1, BS2, BP4",
                         "Incorect intervar value"
                         )
        self.assertEqual(xu.get_cell_value(4, intervar_evidence_col_idx),
                         None,
                         "Incorect intervar value"
                         )
        self.assertEqual(xu.get_cell_value(5, intervar_evidence_col_idx),
                         "BS1, BS2",
                         "Incorect intervar value"
                         )
        self.assertEqual(xu.get_cell_value(6, intervar_evidence_col_idx),
                         "PVS1, PM2, PP5",
                         "Incorect intervar value"
                         )
        self.assertEqual(xu.get_cell_value(7, intervar_evidence_col_idx),
                         "PVS1, PM2",
                         "Incorect intervar value"
                         )


    def tearDown(self):
        self.remove_working_dir()
