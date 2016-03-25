# This Python file uses the following encoding: utf-8
import unittest
from os.path import join as join_path
from os.path import dirname
from collections import OrderedDict
from pycmm.template import SafeTester
from pycmm.utils.xlsutils import XlsUtils
from pycmm.settings import DFLT_MUTREP_FREQ_RATIOS
from pycmm.settings import ALL_MUTREP_ANNO_COLS
from pycmm.settings import MT_COLS_TAG
from pycmm.settings import AXEQ_CHR9_COLS_TAG
from pycmm.settings import CRC_CRC_COLS_TAG
from pycmm.settings import NK64_COLS_TAG
from pycmm.settings import MUTSTAT_DETAILS_COLS_TAG
from pycmm.settings import EXAC_OTH_COLS_TAG
from pycmm.settings import UNKNOWN_COLS_TAG
from pycmm.settings import FUNC_REFGENE_COL_NAME
from pycmm.settings import EXONICFUNC_REFGENE_COL_NAME
from pycmm.settings import GENE_REFGENE_COL_NAME
from pycmm.settings import GENEDETAIL_REFGENE_COL_NAME
from pycmm.settings import CYTOBAND_COL_NAME
from pycmm.settings import KG2014OCT_ALL_COL_NAME
from pycmm.settings import AXEQ_CHR9_HET_COL_NAME
from pycmm.settings import AXEQ_CHR3_6_14_18_PF_COL_NAME
from pycmm.settings import CRC_CAFAM_WT_COL_NAME
from pycmm.settings import AXEQ_CHR5_19_GF_COL_NAME
from pycmm.settings import CRC_CRC_PF_COL_NAME
from pycmm.settings import PRIMARY_MAF_VAR
from pycmm.settings import ILL_BR_PF_COL_NAME
from pycmm.settings import EXAC_ALL_COL_NAME
from pycmm.settings import FULL_SYSTEM_TEST
from pycmm.flow.mutrep import MutRepPipeline
from pycmm.flow.cmmdb import create_jobs_setup_file
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FILTER_RARE
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FILTER_NON_INTERGENIC
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FILTER_NON_INTRONIC
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FILTER_HAS_MUTATION
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FILTER_HAS_SHARED
from pycmm.flow.mutrep import ACTION_DELETE_ROW
from pycmm.flow.mutrep import ACTION_COLOR_ROW
from pycmm.flow.mutrep import ACTION_COLOR_COL

DFLT_TEST_MUTREP_COLS = OrderedDict()
DFLT_TEST_MUTREP_COLS[FUNC_REFGENE_COL_NAME] = ALL_MUTREP_ANNO_COLS[FUNC_REFGENE_COL_NAME]
DFLT_TEST_MUTREP_COLS[EXONICFUNC_REFGENE_COL_NAME] = ALL_MUTREP_ANNO_COLS[EXONICFUNC_REFGENE_COL_NAME]
DFLT_TEST_MUTREP_COLS[GENE_REFGENE_COL_NAME] = ALL_MUTREP_ANNO_COLS[GENE_REFGENE_COL_NAME]
DFLT_TEST_MUTREP_COLS[GENEDETAIL_REFGENE_COL_NAME] = ALL_MUTREP_ANNO_COLS[GENEDETAIL_REFGENE_COL_NAME]
DFLT_TEST_MUTREP_COLS[CYTOBAND_COL_NAME] = ALL_MUTREP_ANNO_COLS[CYTOBAND_COL_NAME]
DFLT_TEST_MUTREP_COLS[KG2014OCT_ALL_COL_NAME] = ALL_MUTREP_ANNO_COLS[KG2014OCT_ALL_COL_NAME]
DFLT_TEST_MUTREP_COLS[AXEQ_CHR9_HET_COL_NAME] = ALL_MUTREP_ANNO_COLS[AXEQ_CHR9_HET_COL_NAME]
DFLT_TEST_MUTREP_COLS[AXEQ_CHR3_6_14_18_PF_COL_NAME] = ALL_MUTREP_ANNO_COLS[AXEQ_CHR3_6_14_18_PF_COL_NAME]
DFLT_TEST_MUTREP_COLS[CRC_CAFAM_WT_COL_NAME] = ALL_MUTREP_ANNO_COLS[CRC_CAFAM_WT_COL_NAME]
DFLT_TEST_MUTREP_COLS[AXEQ_CHR5_19_GF_COL_NAME] = ALL_MUTREP_ANNO_COLS[AXEQ_CHR5_19_GF_COL_NAME]
DFLT_TEST_MUTREP_COLS[CRC_CRC_PF_COL_NAME] = ALL_MUTREP_ANNO_COLS[CRC_CRC_PF_COL_NAME]

DFLT_TEST_REPORT_REGIONS = "18:12512255-14542551"
DFLT_TEST_FREQ_RATIOS = DFLT_MUTREP_FREQ_RATIOS

DFLT_TEST_ANNO_EXCL_TAGS = MT_COLS_TAG
DFLT_TEST_ANNO_EXCL_TAGS += "," + AXEQ_CHR9_COLS_TAG
DFLT_TEST_ANNO_EXCL_TAGS += "," + CRC_CRC_COLS_TAG
DFLT_TEST_ANNO_EXCL_TAGS += "," + NK64_COLS_TAG
DFLT_TEST_ANNO_EXCL_TAGS += "," + MUTSTAT_DETAILS_COLS_TAG
DFLT_TEST_ANNO_EXCL_TAGS += "," + EXAC_OTH_COLS_TAG
DFLT_TEST_ANNO_EXCL_TAGS += "," + UNKNOWN_COLS_TAG

MUTREP_TEST = False

class TestMutRepPipeline(SafeTester):

    def __init__(self, methodName):
        super(TestMutRepPipeline, self).__init__(methodName=methodName,
                                                 test_module_name=__name__,
                                                 )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self,
                                 dataset_name=None,
                                 sample_infos=None,
                                 anno_cols=DFLT_TEST_MUTREP_COLS,
                                 anno_excl_tags=None,
                                 annotated_vcf_tabix=None,
                                 report_regions=DFLT_TEST_REPORT_REGIONS,
                                 frequency_ratios=DFLT_TEST_FREQ_RATIOS,
                                 expression_patterns=None,
                                 expression_usages=None,
                                 split_chrom=False,
                                 summary_families_sheet=False,
                                 call_detail=False,
                                 rows_filter_actions=None,
                                 only_summary=False,
                                 only_families=False,
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        if dataset_name is None:
            dataset_name = self.test_function
        create_jobs_setup_file(dataset_name=dataset_name,
                               project_out_dir=self.working_dir,
                               sample_infos=sample_infos,
                               anno_cols=",".join(anno_cols),
                               anno_excl_tags=anno_excl_tags,
                               annotated_vcf_tabix=annotated_vcf_tabix,
                               report_regions=report_regions,
                               frequency_ratios=frequency_ratios,
                               expression_patterns=expression_patterns,
                               expression_usages=expression_usages,
                               split_chrom=split_chrom,
                               summary_families_sheet=summary_families_sheet,
                               call_detail=call_detail,
                               rows_filter_actions=rows_filter_actions,
                               only_summary=only_summary,
                               only_families=only_families,
                               out_jobs_setup_file=jobs_setup_file,
                               )
        return jobs_setup_file

    def test_load_jobs_info_1(self):
        """ test if layout configurations are loaded correctly """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(annotated_vcf_tabix=annotated_vcf_tabix)
        pl = MutRepPipeline(jobs_setup_file)
        self.assertEqual(pl.report_layout.anno_cols[3],
                         "GeneDetail.refGene",
                         "MutRepPipeline cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(pl.report_layout.anno_cols[4],
                         "cytoBand",
                         "MutRepPipeline cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(pl.report_layout.anno_cols[6],
                         AXEQ_CHR9_HET_COL_NAME,
                         "MutRepPipeline cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(len(pl.report_layout.anno_cols),
                         11,
                         "MutRepPipeline cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(len(pl.report_layout.anno_excl_tags),
                         0,
                         "MutRepPipeline cannot correctly read report layout info 'annotation excluded tags' from jobs setup file")
        self.assertEqual(pl.report_layout.report_regions[0].chrom,
                         "18",
                         "MutRepPipeline cannot correctly read report layout info 'report regions' from jobs setup file")
        self.assertFalse(pl.report_layout.exprs,
                         "MutRepPipeline cannot correctly read report layout info 'expressions' from jobs setup file")
        self.assertFalse(pl.report_layout.split_chrom,
                         "MutRepPipeline cannot correctly read report layout info 'split chrom' from jobs setup file")
        self.assertFalse(pl.report_layout.summary_families_sheet,
                         "MutRepPipeline cannot correctly read report layout info 'summary families sheet' from jobs setup file")
        self.assertFalse(pl.report_layout.call_detail,
                         "MutRepPipeline cannot correctly read report layout info 'call info' from jobs setup file")
        self.assertFalse(pl.report_layout.filter_rare,
                         "MutRepPipeline cannot correctly read report layout info 'filter rare' from jobs setup file")
        self.assertFalse(pl.report_layout.filter_non_intergenic,
                         "MutRepPipeline cannot correctly read report layout info 'filter non-intergenic' from jobs setup file")
        self.assertFalse(pl.report_layout.filter_non_intronic,
                         "MutRepPipeline cannot correctly read report layout info 'filter non-intronic' from jobs setup file")
        self.assertFalse(pl.report_layout.filter_has_mutation,
                         "MutRepPipeline cannot correctly read report layout info 'filter has-mutation' from jobs setup file")
        self.assertFalse(pl.report_layout.only_summary,
                         "MutRepPipeline cannot correctly read report layout info 'only summary' from jobs setup file")
        self.assertFalse(pl.report_layout.only_families,
                         "MutRepPipeline cannot correctly read report layout info 'only families' from jobs setup file")

    def test_load_jobs_info_2(self):
        """ test if non-default layout configurations are loaded correctly """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        dummy_annotated_vcf_tabix = join_path(self.data_dir,
                                              "input.vcf.gz")
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_RARE
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_INTERGENIC
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_INTRONIC
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_HAS_MUTATION
        jobs_setup_file = self.__create_jobs_setup_file(anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        annotated_vcf_tabix=dummy_annotated_vcf_tabix,
                                                        report_regions="6:78161823-78164117,"+DFLT_TEST_REPORT_REGIONS+",22",
                                                        frequency_ratios=None,
                                                        expression_patterns="test1:1>0",
                                                        expression_usages="test1:del_row;test1:color_col:cytoBand:yellow",
                                                        split_chrom=True,
                                                        summary_families_sheet=True,
                                                        call_detail=True,
                                                        rows_filter_actions=rows_filter_actions,
                                                        only_summary=True,
                                                        only_families=True,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        self.assertEqual(pl.report_layout.anno_cols[3],
                         "GeneDetail.refGene",
                         "MutRepPipeline cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(pl.report_layout.anno_cols[4],
                         "cytoBand",
                         "MutRepPipeline cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(pl.report_layout.anno_cols[7],
                         AXEQ_CHR5_19_GF_COL_NAME,
                         "MutRepPipeline cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(len(pl.report_layout.anno_cols),
                         8,
                         "MutRepPipeline cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(len(pl.report_layout.anno_excl_tags),
                         7,
                         "MutRepPipeline cannot correctly read report layout info 'annotation excluded tags' from jobs setup file")
        self.assertEqual(pl.report_layout.report_regions[1].end_pos,
                         "14542551",
                         "MutRepPipeline cannot correctly read report layout info 'report regions' from jobs setup file")
        self.assertEqual(pl.report_layout.report_regions[2].chrom,
                         "22",
                         "MutRepPipeline cannot correctly read report layout info 'report regions' from jobs setup file")
        self.assertEqual(pl.report_layout.report_regions[2].start_pos,
                         None,
                         "MutRepPipeline cannot correctly read report layout info 'report regions' from jobs setup file")
        self.assertEqual(pl.annotated_vcf_tabix,
                         dummy_annotated_vcf_tabix,
                         "MutRepPipeline cannot correctly determine report layout info 'annotated vcf tabix' file")
        self.assertEqual(pl.report_layout.exprs.patterns["test1"],
                         '1>0',
                         "MutRepPipeline cannot correctly read report layout info 'expression patterns' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_DELETE_ROW][0].pattern,
                         '1>0',
                         "MutRepPipeline cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_COL][0].pattern,
                         '1>0',
                         "MutRepPipeline cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_COL][0].col_name,
                         'cytoBand',
                         "MutRepPipeline cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_COL][0].color,
                         'yellow',
                         "MutRepPipeline cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertTrue(pl.report_layout.split_chrom,
                        "MutRepPipeline cannot correctly read report layout info 'split chrom' from jobs setup file")
        self.assertTrue(pl.report_layout.summary_families_sheet,
                        "MutRepPipeline cannot correctly read report layout info 'summary families sheet' from jobs setup file")
        self.assertTrue(pl.report_layout.call_detail,
                        "MutRepPipeline cannot correctly read report layout info 'call info' from jobs setup file")
        self.assertTrue(pl.report_layout.filter_rare,
                        "MutRepPipeline cannot correctly read report layout info 'filter rare' from jobs setup file")
        self.assertTrue(pl.report_layout.filter_non_intergenic,
                        "MutRepPipeline cannot correctly read report layout info 'filter non-intergenic' from jobs setup file")
        self.assertTrue(pl.report_layout.filter_non_intronic,
                        "MutRepPipeline cannot correctly read report layout info 'filter non-intronic' from jobs setup file")
        self.assertTrue(pl.report_layout.filter_has_mutation,
                        "MutRepPipeline cannot correctly read report layout info 'filter has-mutation' from jobs setup file")
        self.assertTrue(pl.report_layout.only_summary,
                        "MutRepPipeline cannot correctly read report layout info 'only summary' from jobs setup file")
        self.assertTrue(pl.report_layout.only_families,
                        "MutRepPipeline cannot correctly read report layout info 'only families' from jobs setup file")

    def test_load_jobs_info_3(self):
        """ test if non-default (None) layout configurations are loaded correctly """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(report_regions=None,
                                                        frequency_ratios="ExAC:0.5",
                                                        expression_patterns='expr_with_key:"abc" < 4;expr_wo:jkl>5; expr78 : 2>3',
                                                        expression_usages="expr_with_key:del_row;expr_wo:color_col:F_U_8:red;expr_wo:color_row:green;expr78:color_row:orange",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        self.assertEqual(pl.report_layout.report_regions,
                         None,
                         "MutRepPipeline cannot correctly read report layout info 'report regions' from jobs setup file")
        self.assertEqual(pl.report_layout.freq_ratios["ExAC"],
                         0.5,
                         "MutRepPipeline cannot correctly read report layout info 'frequency ratios' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.patterns["expr_with_key"],
                         '"abc" < 4',
                         "MutRepPipeline cannot correctly read report layout info 'expressions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.patterns["expr_wo"],
                         'jkl>5',
                         "MutRepPipeline cannot correctly read report layout info 'expressions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_DELETE_ROW][0].pattern,
                         '"abc" < 4',
                         "MutRepPipeline cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_COL][0].pattern,
                         'jkl>5',
                         "MutRepPipeline cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_COL][0].col_name,
                         'F_U_8',
                         "MutRepPipeline cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_COL][0].color,
                         'red',
                         "MutRepPipeline cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_ROW][0].pattern,
                         'jkl>5',
                         "MutRepPipeline cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_ROW][0].color,
                         'green',
                         "MutRepPipeline cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_ROW][1].pattern,
                         '2>3',
                         "MutRepPipeline cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_ROW][1].color,
                         'orange',
                         "MutRepPipeline cannot correctly read report layout info 'expression actions' from jobs setup file")

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_summary_report_1(self):
        """ test if summary report with default configuration can be correctly generated """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "chr6_18.vcf.gz")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions=None,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         18,
                         "shared mutations cannot be correctly determined")

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_summary_report_2(self):
        """ test summary with multiple report_regions """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "chr6_18.vcf.gz")
        rpt_out_file = join_path(self.working_dir,
                                 self.current_func_name + ".xlsx")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions="6:78171940-78172992,18:28610987-28611790",
                                                        call_detail="YES",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions, out_file=rpt_out_file)
        xu = XlsUtils(rpt_out_file)
        self.assertEqual(xu.count_rows(),
                         10,
                         "shared mutations cannot be correctly determined")

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_summary_report_3(self):
        """
        test if mutations with more than one alternate alleles in summary
        is correctly generated
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions="6",
                                                        call_detail=False,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertFalse(xu.compare_vals("AXEQ_CHR3_6_14_18_PF",
                                         5,
                                         6,
                                         ),
                         "information of mutations with more than one alternate alleles are incorrect"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_summary_report_4(self):
        """ test summary with multiple report_regions and many sample infos """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "chr6_18.vcf.gz")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions="6:78171940-78172992,18:28610987-28611790",
                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
                                                        summary_families_sheet=True,
                                                        call_detail="YES",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=1),
                         10,
                         "Incorrect number of rows"
                         )
        self.assertEqual(xu.count_cols(col_name1="1234-Alb-31",
                                       col_name2="6789-Al-65",
                                       sheet_idx=1),
                         11,
                         "Incorrect number of columns"
                         )
        self.assertEqual(xu.nsheets,
                         2,
                         "Incorrect number of sheets"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_summary_report_5(self):
        """ test summary report of CRC samples """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions=None,
                                                        call_detail=False,
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         15,
                         "CRC report cannot be generated correctly"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_summary_report_6(self):
        """ test if unicode character 'ö' is allowed in the report """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_NON_INTRONIC
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        report_regions=None,
                                                        call_detail=False,
                                                        rows_filter_actions=rows_filter_actions,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         2,
                         "report with swedish unicode character cannot be generated correctly"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_family_report_1(self):
        """ test with only one family which has only one members """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        call_detail="YES",
                                                        report_regions="6",
                                                        sample_infos="6789:Al-65",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_family_report('6789', pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam6789.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.nsheets,
                         1,
                         "invalid number of sheets")
        self.assertEqual(xu.get_sheet_idx("Al-65"),
                         0,
                         "invalid sheet name")

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_family_report_2(self):
        """ test with only one family which has two members """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        call_detail="YES",
                                                        report_regions="6",
                                                        sample_infos="1234:Alb-31:Br-466",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_family_report('1234', pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam1234.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.nsheets,
                         3,
                         "invalid number of sheets")
        self.assertEqual(xu.get_sheet_idx("shared"),
                         0,
                         "invalid sheet name")
        self.assertEqual(xu.get_sheet_idx("Alb-31"),
                         1,
                         "invalid sheet name")
        self.assertEqual(xu.get_sheet_idx("Br-466"),
                         2,
                         "invalid sheet name")

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_family_report_3(self):
        """ test with only one family which has three members """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        call_detail="YES",
                                                        report_regions="6",
                                                        sample_infos="6067:Br-432:Al-161:Br-504",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_family_report('6067', pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam6067.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.nsheets,
                         4,
                         "invalid number of sheets")
        self.assertEqual(xu.get_sheet_idx("shared"),
                         0,
                         "invalid sheet name")
        self.assertEqual(xu.get_sheet_idx("Br-432"),
                         1,
                         "invalid sheet name")
        self.assertEqual(xu.get_sheet_idx("Al-161"),
                         2,
                         "invalid sheet name")
        self.assertEqual(xu.get_sheet_idx("Br-504"),
                         3,
                         "invalid sheet name")

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_family_report_4(self):
        """
        test with number of mutations are correct in each tab
        - shared
        - member 1
        - member 2
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        dataset_name = self.test_function
        frequency_ratios = PRIMARY_MAF_VAR + ":0.2"
        frequency_ratios += "," + ILL_BR_PF_COL_NAME + ":0.3"
        frequency_ratios += "," + EXAC_ALL_COL_NAME + ":0.3"
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        call_detail="YES",
                                                        report_regions=None,
                                                        frequency_ratios=frequency_ratios,
                                                        sample_infos="24:Co-166:Co-213",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_family_report('24', pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam24.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         8,
                         "incorrect number of mutations")
        self.assertEqual(xu.count_rows(sheet_idx=1),
                         11,
                         "incorrect number of mutations")
        self.assertEqual(xu.count_rows(sheet_idx=2),
                         13,
                         "incorrect number of mutations")

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_family_report_5(self):
        """ test if report can be run with missing columns """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions=None,
                                                        anno_cols=ALL_MUTREP_ANNO_COLS,
                                                        sample_infos="24:Co-166:Co-213",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_family_report('24', pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam24.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         8,
                         "incorrect number of mutations")
        self.assertEqual(xu.count_rows(sheet_idx=1),
                         11,
                         "incorrect number of mutations")
        self.assertEqual(xu.count_rows(sheet_idx=2),
                         13,
                         "incorrect number of mutations")

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_families_reports_1(self):
        """ test generating offline families reports (1 family)"""

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        call_detail="YES",
                                                        report_regions="6",
                                                        sample_infos="6067:Br-432:Al-161:Br-504",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_families_reports()
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam6067.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.nsheets,
                         4,
                         "invalid number of sheets")
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         4,
                         "incorrect number of mutations")
        self.assertEqual(xu.count_rows(sheet_idx=1),
                         6,
                         "incorrect number of mutations")
        self.assertEqual(xu.count_rows(sheet_idx=2),
                         7,
                         "incorrect number of mutations")
        self.assertEqual(xu.count_rows(sheet_idx=3),
                         7,
                         "incorrect number of mutations")

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_families_reports_2(self):
        """ test generating offline families reports (3 families)"""

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions="6",
                                                        call_detail="YES",
                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_families_reports()
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam1234.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.nsheets,
                         3,
                         "invalid number of sheets")
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam6067.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.nsheets,
                         4,
                         "invalid number of sheets")
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam6789.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.nsheets,
                         1,
                         "invalid number of sheets")

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_families_reports_3(self):
        """ test generating 101 CRC families reports """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions=None,
                                                        call_detail=False,
                                                        sample_infos="8:Co-35:Co-37,12:Co-89:Co-90,13:Co-95,26:Co-131:Co-135,31:1793-11D:1322-11D,87:Co-218:Co-258,91:Co-454:Co-700,94:Co-238,110:1526-02D:Co-1301,134:Co-460:Co-553,141:Co-305:Co-785,185:Co-603:Co-669,191:Co-384,214:Co-484,216:Co-367:Co-446,221:Co-358,227:Co-364,231:Co-555:Co-572,254:Co-616:Co-1156,275:Co-618:Co-1262,288:Co-1141,296:Co-793:Co-876,301:Co-837:Co-840:Co-1053,306:Co-779,309:Co-783,312:Co-1116,315:1462-01D,325:Co-851:Co-859,348:Co-846:Co-857,350:1104-03D:Co-866,409:Co-1254,415:Co-1031:Co-1037,425:Co-1458:Co-1595,434:Co-1051:Co-1534,445:Co-1157:Co-1158,478:Co-1207:Co-1274,485:Co-1302:Co-1322,532:Co-1583:Co-1584,574:468-04:474-05,578:531-04o:Co-1349,650:398-05o:729-05o,695:Co-1354:Co-1359:Co-1368,739:529-05:Co-1467,740:602-05o:Co-1373:Co-1383,849:Co-1764:Co-1765,869:Co-1685,871:Co-1618:Co-1661,918:134-06:354-06,975:Co-1591:Co-1600,1025:Co-1529,1085:Co-1518,1113:642-06:Co-1538,1206:1052-05D:Co-1552,1207:2818-07D,1213:Co-1666,1252:Co-1719,1290:Co-1723,prostate:P001:P002:P003",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_families_reports()
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam869.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.nsheets,
                         1,
                         "invalid number of sheets")
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam13.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.nsheets,
                         1,
                         "invalid number of sheets")
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam425.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.nsheets,
                         3,
                         "invalid number of sheets")
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam301.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.nsheets,
                         4,
                         "invalid number of sheets")

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_families_reports_4(self):
        """ test generating NK64 fam24 families reports """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions=None,
                                                        call_detail=False,
                                                        sample_infos="24:Co-166:Co-213",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_families_reports()
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam24.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.nsheets,
                         3,
                         "invalid number of sheets")
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         2,
                         "incorrect number of mutations")
        self.assertEqual(xu.count_rows(sheet_idx=1),
                         5,
                         "incorrect number of mutations")
        self.assertEqual(xu.count_rows(sheet_idx=2),
                         3,
                         "incorrect number of mutations")

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_families_reports_5(self):
        """ test generating NK64 PMS2 families reports """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        anno_cols = list(DFLT_TEST_MUTREP_COLS)
        anno_cols.append("LRT_score")
        anno_cols.append("LRT_pred")
        anno_cols.append("Polyphen2_HDIV_score")
        anno_cols.append("Polyphen2_HDIV_pred")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        anno_cols=anno_cols,
                                                        report_regions=None,
                                                        call_detail=False,
                                                        sample_infos="MYUTH:2014-06388-02,119:Co-1209:Co-1220:Co-1222:Co-1285,1244:2310-08D:627-04o:3015-11D,1330:861-10D,1410:930-06o,1680:301-08F",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_families_reports()
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_famMYUTH.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.nsheets,
                         1,
                         "invalid number of sheets")
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam119.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.nsheets,
                         5,
                         "invalid number of sheets")
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam1330.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.nsheets,
                         1,
                         "invalid number of sheets")
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_fam1680.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.nsheets,
                         1,
                         "invalid number of sheets")

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_expression_action_del_row_1(self):
        """
        test if expreesion patterns and expression actions can be used
        for deleting rows with a given pattern
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions="6:78171940-78172992,18:28610987-28611790",
                                                        expression_patterns='exp_ns_snv:"ExonicFunc.refGene" == \'synonymous_SNV\'',
                                                        expression_usages="exp_ns_snv:del_row",
                                                        call_detail="YES",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         8,
                         "Incorrect number of rows"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_expression_action_color_row_1(self):
        """
        test if expreesion patterns and expression actions can be used
        for coloring rows with a given pattern
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions="6:78171940-78172992,18:28610987-28611790",
                                                        expression_patterns='exp_splicing:"Func.refGene" == \'splicing\'',
                                                        expression_usages="exp_splicing:color_row:ROSY_BROWN",
                                                        call_detail="YES",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        col_idx = xu.get_col_idx(col_name="Func.refGene", sheet_idx=0)
        self.assertEqual(xu.get_cell_value(8, col_idx, sheet_idx=0),
                         "splicing",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(8, 2, sheet_idx=0),
                         "FFBC8F8F",
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_rgb(8, col_idx, sheet_idx=0),
                         "FFBC8F8F",
                         "Incorrect color"
                         )
        self.assertNotEqual(xu.get_cell_rgb(7, col_idx, sheet_idx=0),
                            "FFBC8F8F",
                            "Incorrect color"
                            )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_expression_action_color_col_1(self):
        """
        test if expreesion patterns and expression actions can be used
        for coloring rows with a given pattern
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions="6:78171940-78172992,18:28610987-28611790",
                                                        expression_patterns='exp_exonic:"Func.refGene" == \'exonic\';exp_ns_snv:"ExonicFunc.refGene" == \'nonsynonymous_SNV\'',
                                                        expression_usages="exp_exonic:color_row:ROSY_BROWN;exp_ns_snv:color_col:cytoBand:TEAL;exp_exonic:color_col:Gene.refGene:YELLOW",
                                                        call_detail="YES",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        col_idx = xu.get_col_idx(col_name="Func.refGene", sheet_idx=0)
        self.assertEqual(xu.get_cell_value(4, col_idx, sheet_idx=0),
                         "exonic",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(5, 14, sheet_idx=0),
                         "FFBC8F8F",
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_rgb(5, col_idx, sheet_idx=0),
                         "FFBC8F8F",
                         "Incorrect color"
                         )
        self.assertNotEqual(xu.get_cell_rgb(2, col_idx, sheet_idx=0),
                            "FFBC8F8F",
                            "Incorrect color"
                            )
        col_idx = xu.get_col_idx(col_name="cytoBand", sheet_idx=0)
        self.assertEqual(xu.get_cell_rgb(4, col_idx, sheet_idx=0),
                         "FF008080",
                         "Incorrect color"
                         )
        col_idx = xu.get_col_idx(col_name="Gene.refGene", sheet_idx=0)
        self.assertEqual(xu.get_cell_rgb(4, col_idx, sheet_idx=0),
                         "FFFFFF00",
                         "Incorrect color"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST, "taking too long time to test")
    def test_color_shared_mutations_1(self):
        """
        test if shared mutation can be colored correctly
        """

        self.init_test(self.current_func_name)
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        dataset_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file(dataset_name=dataset_name,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions="9:99700709-99702632",
                                                        sample_infos="8:Co-35:Co-37,13:Co-95,275:Co-1262:Co-618,296:Co-793:Co-876",
                                                        summary_families_sheet=True,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)
        xls_file = join_path(self.working_dir,
                             "rpts",
                             dataset_name+"_summary.xlsx")
        xu = XlsUtils(xls_file)
        col_idx = xu.get_col_idx(col_name="275-Co-1262", sheet_idx=1)
        self.assertEqual(xu.get_cell_value(5, col_idx, sheet_idx=1),
                         "het",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(5, col_idx, sheet_idx=1),
                         "FFC0C0C0",
                         "Incorrect color"
                         )

    def tearDown(self):
        self.remove_working_dir()
