# This Python file uses the following encoding: utf-8
import unittest
from os.path import join as join_path
#from os.path import dirname
from collections import OrderedDict
from pycmm.template import SafeTester
from pycmm.cmmlib.xlslib import XlsUtils
from pycmm.cmmlib.colorlib import COLORS_RGB
from pycmm.settings import XLS_TEST
from pycmm.settings import DB_TEST
#from pycmm.settings import DFLT_MUTREP_FREQ_RATIOS
from pycmm.settings import ALL_MUTREP_ANNO_COLS
#from pycmm.settings import MT_COLS_TAG
from pycmm.settings import AXEQ_CHR9_COLS_TAG
from pycmm.settings import AXEQ_CHR3_6_14_18_COLS_TAG
from pycmm.settings import AXEQ_CHR5_19_COLS_TAG
from pycmm.settings import MUTSTAT_DETAILS_COLS_TAG
from pycmm.settings import EXAC_OTH_COLS_TAG
from pycmm.settings import EXAC_CONSTRAINT_COLS_TAG
from pycmm.settings import LJB_SCORE_COLS_TAG
from pycmm.settings import UNKNOWN_COLS_TAG
from pycmm.settings import FUNC_REFGENE_COL_NAME
from pycmm.settings import EXONICFUNC_REFGENE_COL_NAME
from pycmm.settings import GENE_REFGENE_COL_NAME
from pycmm.settings import GENEDETAIL_REFGENE_COL_NAME
from pycmm.settings import CYTOBAND_COL_NAME
from pycmm.settings import SPIDEX_DPSI_ZSCORE_COL_NAME
from pycmm.settings import PATHOGENIC_COUNT_COL_NAME
from pycmm.settings import PREDICTION_COLS
from pycmm.settings import LJB_METASVM_PREDICTION_COL_NAME
from pycmm.settings import LJB_METALR_PREDICTION_COL_NAME
from pycmm.settings import KG2014OCT_ALL_COL_NAME
#from pycmm.settings import KG2014OCT_EUR_COL_NAME
from pycmm.settings import AXEQ_CHR9_HET_COL_NAME
from pycmm.settings import AXEQ_CHR3_6_14_18_PF_COL_NAME
from pycmm.settings import AXEQ_CHR5_19_GF_COL_NAME
#from pycmm.settings import PRIMARY_MAF_VAR
#from pycmm.settings import EXAC_ALL_COL_NAME
from pycmm.settings import FULL_SYSTEM_TEST
from pycmm.settings import WES294_OAF_EARLYONSET_AF_COL_NAME
from pycmm.settings import WES294_OAF_EARLYONSET_GF_COL_NAME
#from pycmm.settings import EST_KVOT_EARLYONSET_VS_BRC_COL_NAME
#from pycmm.settings import EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME
#from pycmm.settings import EST_KVOT_EARLYONSET_VS_KG_EUR_COL_NAME
from pycmm.settings import WES294_OAF_BRCS_AF_COL_NAME
from pycmm.settings import SWEGEN_AF_COL_NAME
from pycmm.settings import MAX_REF_MAF_COL_NAME
from pycmm.settings import REF_MAF_COL_NAMES
from pycmm.settings import COMPOUND_HETEROZYGOTE_AFFECTED_COUNT_COL_NAME
from pycmm.settings import COMPOUND_HETEROZYGOTE_FREQ_RATIO_COL_NAME
from pycmm.settings import HOMOZYGOTE_AFFECTED_COUNT_COL_NAME
from pycmm.settings import HOMOZYGOTE_FREQ_RATIO_COL_NAME
from pycmm.proc.mutrep.mutrep import MutRepController
from pycmm.proc.mutrep.mutrep import create_jobs_setup_file
#from pycmm.proc.mutrep.mutrep import JOBS_SETUP_RPT_FILTER_RARE
from pycmm.proc.mutrep.mutrep import JOBS_SETUP_RPT_FILTER_PASS_VQSR
from pycmm.proc.mutrep.mutrep import JOBS_SETUP_RPT_FILTER_NON_INTERGENIC
from pycmm.proc.mutrep.mutrep import JOBS_SETUP_RPT_FILTER_NON_INTRONIC
from pycmm.proc.mutrep.mutrep import JOBS_SETUP_RPT_FILTER_NON_UPSTREAM
from pycmm.proc.mutrep.mutrep import JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM
from pycmm.proc.mutrep.mutrep import JOBS_SETUP_RPT_FILTER_NON_UTR
from pycmm.proc.mutrep.mutrep import JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS
from pycmm.proc.mutrep.mutrep import JOBS_SETUP_RPT_FILTER_HAS_MUTATION
from pycmm.proc.mutrep.mutrep import JOBS_SETUP_RPT_FILTER_HAS_SHARED
from pycmm.proc.mutrep.mutrep import JOBS_SETUP_RPT_FILTER_NON_RECESSIVE_GENE
from pycmm.proc.mutrep.mutrep import ACTION_DELETE_ROW
from pycmm.proc.mutrep.mutrep import ACTION_COLOR_ROW
from pycmm.proc.mutrep.mutrep import ACTION_COLOR_COL

DFLT_TEST_MUTREP_COLS = OrderedDict()
DFLT_TEST_MUTREP_COLS[FUNC_REFGENE_COL_NAME] = ALL_MUTREP_ANNO_COLS[FUNC_REFGENE_COL_NAME]
DFLT_TEST_MUTREP_COLS[EXONICFUNC_REFGENE_COL_NAME] = ALL_MUTREP_ANNO_COLS[EXONICFUNC_REFGENE_COL_NAME]
DFLT_TEST_MUTREP_COLS[GENE_REFGENE_COL_NAME] = ALL_MUTREP_ANNO_COLS[GENE_REFGENE_COL_NAME]
DFLT_TEST_MUTREP_COLS[GENEDETAIL_REFGENE_COL_NAME] = ALL_MUTREP_ANNO_COLS[GENEDETAIL_REFGENE_COL_NAME]
DFLT_TEST_MUTREP_COLS[CYTOBAND_COL_NAME] = ALL_MUTREP_ANNO_COLS[CYTOBAND_COL_NAME]
DFLT_TEST_MUTREP_COLS[KG2014OCT_ALL_COL_NAME] = ALL_MUTREP_ANNO_COLS[KG2014OCT_ALL_COL_NAME]
DFLT_TEST_MUTREP_COLS[AXEQ_CHR9_HET_COL_NAME] = ALL_MUTREP_ANNO_COLS[AXEQ_CHR9_HET_COL_NAME]
DFLT_TEST_MUTREP_COLS[AXEQ_CHR3_6_14_18_PF_COL_NAME] = ALL_MUTREP_ANNO_COLS[AXEQ_CHR3_6_14_18_PF_COL_NAME]
DFLT_TEST_MUTREP_COLS[AXEQ_CHR5_19_GF_COL_NAME] = ALL_MUTREP_ANNO_COLS[AXEQ_CHR5_19_GF_COL_NAME]

DFLT_TEST_REPORT_REGIONS = "18:12512255-14542551"

DFLT_TEST_ANNO_EXCL_TAGS = AXEQ_CHR9_COLS_TAG
DFLT_TEST_ANNO_EXCL_TAGS += "," + MUTSTAT_DETAILS_COLS_TAG
DFLT_TEST_ANNO_EXCL_TAGS += "," + EXAC_OTH_COLS_TAG
DFLT_TEST_ANNO_EXCL_TAGS += "," + EXAC_CONSTRAINT_COLS_TAG
DFLT_TEST_ANNO_EXCL_TAGS += "," + UNKNOWN_COLS_TAG

MUTREP_TEST = False

RGB_NO_FILL = "00000000"

class TestMutRepController(SafeTester):

    def __init__(self, methodName):
        super(TestMutRepController, self).__init__(methodName=methodName,
                                                   test_module_name=__name__,
                                                   )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self, *args, **kwargs):
        if 'project_name' not in kwargs:
            kwargs['project_name'] = self.test_function
        if 'db_file' not in kwargs:
            kwargs['db_file'] = join_path(self.data_dir, "input.db")
        kwargs['project_out_dir'] = self.working_dir
        kwargs['jobs_setup_file'] = join_path(self.working_dir,
                                              kwargs['project_name']+'_jobs_setup.txt')
        create_jobs_setup_file(*args, **kwargs)
        return kwargs['jobs_setup_file']

    def test_load_jobs_info_1(self):
        """ test if default layout configurations are loaded correctly """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            "input.db")
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file)
        pl = MutRepController(jobs_setup_file)
        self.assertEqual(len(pl.report_layout.anno_excl_tags),
                         0,
                         "MutRepController cannot correctly read report layout info 'annotation excluded tags' from jobs setup file")
        self.assertFalse(pl.report_layout.exprs,
                         "MutRepController cannot correctly read report layout info 'expressions' from jobs setup file")
        self.assertFalse(pl.report_layout.split_chrom,
                         "MutRepController cannot correctly read report layout info 'split chrom' from jobs setup file")
#        self.assertFalse(pl.report_layout.call_detail,
#                         "MutRepController cannot correctly read report layout info 'call info' from jobs setup file")
#        self.assertFalse(pl.report_layout.call_gq,
#                         "MutRepController cannot correctly read report layout info 'call info' from jobs setup file")
        self.assertFalse(pl.report_layout.show_shared_variants,
                         "MutRepController cannot correctly read report layout info 'show shared mutations' from jobs setup file")
#        self.assertFalse(pl.report_layout.filter_rare,
#                         "MutRepController cannot correctly read report layout info 'filter rare' from jobs setup file")
        self.assertFalse(pl.report_layout.filter_pass_vqsr,
                         "MutRepController cannot correctly read report layout info 'filter pass vqsr' from jobs setup file")
        self.assertFalse(pl.report_layout.filter_non_intergenic,
                         "MutRepController cannot correctly read report layout info 'filter non-intergenic' from jobs setup file")
        self.assertFalse(pl.report_layout.filter_non_intronic,
                         "MutRepController cannot correctly read report layout info 'filter non-intronic' from jobs setup file")
        self.assertFalse(pl.report_layout.filter_has_mutation,
                         "MutRepController cannot correctly read report layout info 'filter has-mutation' from jobs setup file")
        self.assertFalse(pl.report_layout.filter_non_recessive_gene,
                         "MutRepController cannot correctly read report layout info 'filter non-recessive-gene' from jobs setup file")
        self.assertTrue(pl.report_layout.filter_genes is None,
                        "MutRepController cannot correctly read report layout info 'filter genes' from jobs setup file")
        self.assertTrue(pl.report_layout.color_genes is None,
                        "MutRepController cannot correctly read report layout info 'color genes' from jobs setup file")
        self.assertFalse(pl.report_layout.coloring_shared,
                         "MutRepController cannot correctly read report layout info 'coloring shared' from jobs setup file")
        self.assertFalse(pl.report_layout.coloring_zygosity,
                         "MutRepController cannot correctly read report layout info 'coloring zygosity' from jobs setup file")

    def test_load_jobs_info_2(self):
        """ test if non-default layout configurations are loaded correctly """

        self.init_test(self.current_func_name)
        dummy_db_file = join_path(self.data_dir,
                                  "input.db")
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_PASS_VQSR
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_INTERGENIC
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_INTRONIC
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_HAS_MUTATION
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_RECESSIVE_GENE
        jobs_setup_file = self.__create_jobs_setup_file(anno_cols=",".join(DFLT_TEST_MUTREP_COLS),
                                                        anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        header_corrections="old_header:new_header",
                                                        db_file=dummy_db_file,
                                                        report_regions="6:78161823-78164117,"+DFLT_TEST_REPORT_REGIONS+",22",
                                                        frequency_ratios=None,
                                                        filter_genes="CHEK2",
                                                        color_genes="BRCA1",
                                                        expression_patterns="test1:1>0",
                                                        expression_usages="test1:DELETE_ROW,test1:COLOR_COLUMN:cytoBand:yellow",
                                                        split_chrom=True,
                                                        call_detail=True,
                                                        call_gq=True,
                                                        show_shared_variants=True,
                                                        rows_filter_actions=rows_filter_actions,
                                                        coloring_shared=True,
                                                        coloring_zygosity=True,
                                                        )
        pl = MutRepController(jobs_setup_file)
        self.assertEqual(pl.report_layout.anno_cols[3],
                         "GeneDetail.refGene",
                         "MutRepController cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(pl.report_layout.anno_cols[4],
                         "cytoBand",
                         "MutRepController cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(pl.report_layout.anno_cols[6],
                         AXEQ_CHR5_19_GF_COL_NAME,
                         "MutRepController cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(len(pl.report_layout.anno_cols),
                         7,
                         "MutRepController cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(len(pl.report_layout.anno_excl_tags),
                         5,
                         "MutRepController cannot correctly read report layout info 'annotation excluded tags' from jobs setup file")
        self.assertEqual(pl.report_layout.report_regions[1].end_pos,
                         "14542551",
                         "MutRepController cannot correctly read report layout info 'report regions' from jobs setup file")
        self.assertEqual(pl.report_layout.report_regions[2].chrom,
                         "22",
                         "MutRepController cannot correctly read report layout info 'report regions' from jobs setup file")
        self.assertEqual(pl.report_layout.report_regions[2].start_pos,
                         None,
                         "MutRepController cannot correctly read report layout info 'report regions' from jobs setup file")
        self.assertEqual(pl.report_layout.db_file,
                         dummy_db_file,
                         "MutRepController cannot correctly determine report layout info 'annotated vcf tabix' file")
        self.assertEqual(pl.report_layout.exprs.patterns["test1"],
                         '1>0',
                         "MutRepController cannot correctly read report layout info 'expression patterns' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_DELETE_ROW][0].pattern,
                         '1>0',
                         "MutRepController cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_COL]['cytoBand'][0].pattern,
                         '1>0',
                         "MutRepController cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_COL]['cytoBand'][0].col_name,
                         'cytoBand',
                         "MutRepController cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_COL]['cytoBand'][0].color,
                         'yellow',
                         "MutRepController cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertTrue(pl.report_layout.split_chrom,
                        "MutRepController cannot correctly read report layout info 'split chrom' from jobs setup file")
#        self.assertTrue(pl.report_layout.call_detail,
#                        "MutRepController cannot correctly read report layout info 'call info' from jobs setup file")
#        self.assertTrue(pl.report_layout.call_gq,
#                        "MutRepController cannot correctly read report layout info 'call info' from jobs setup file")
        self.assertTrue(pl.report_layout.show_shared_variants,
                        "MutRepController cannot correctly read report layout info 'show shared mutations' from jobs setup file")
#        self.assertTrue(pl.report_layout.filter_rare,
#                        "MutRepController cannot correctly read report layout info 'filter rare' from jobs setup file")
        self.assertTrue(pl.report_layout.filter_pass_vqsr,
                        "MutRepController cannot correctly read report layout info 'filter pass vqsr' from jobs setup file")
        self.assertTrue(pl.report_layout.filter_non_intergenic,
                        "MutRepController cannot correctly read report layout info 'filter non-intergenic' from jobs setup file")
        self.assertTrue(pl.report_layout.filter_non_intronic,
                        "MutRepController cannot correctly read report layout info 'filter non-intronic' from jobs setup file")
        self.assertTrue(pl.report_layout.filter_has_mutation,
                        "MutRepController cannot correctly read report layout info 'filter has-mutation' from jobs setup file")
        self.assertTrue(pl.report_layout.filter_non_recessive_gene,
                        "MutRepController cannot correctly read report layout info 'filter non-recessive-gene' from jobs setup file")
        self.assertEqual(pl.report_layout.filter_genes[0],
                         'CHEK2',
                         "MutRepController cannot correctly read report layout info 'filter genes' from jobs setup file")
        self.assertEqual(pl.report_layout.color_genes[0],
                         'BRCA1',
                         "MutRepController cannot correctly read report layout info 'color genes' from jobs setup file")
        self.assertTrue(pl.report_layout.coloring_shared,
                        "MutRepController cannot correctly read report layout info 'coloring shared' from jobs setup file")
        self.assertTrue(pl.report_layout.coloring_zygosity,
                        "MutRepController cannot correctly read report layout info 'coloring zygosity' from jobs setup file")

    def test_load_jobs_info_3(self):
        """ test if non-default (None) layout configurations are loaded correctly """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(report_regions=None,
                                                        header_corrections="old_header1:new_header1,old_header2:new_header2",
                                                        frequency_ratios="ExAC:0.5",
                                                        filter_genes="CHEK2,NOTCH1,NOTCH4",
                                                        color_genes="BRCA1,MSH2,SRC",
                                                        expression_patterns='expr_with_key:"abc" < 4,expr_wo:jkl>5, expr78 : 2>3',
                                                        expression_usages="expr_with_key:DELETE_ROW,expr_wo:COLOR_COLUMN:F_U_8:red,expr_wo:COLOR_ROW:green,expr78:COLOR_ROW:orange",
                                                        )
        pl = MutRepController(jobs_setup_file=jobs_setup_file)
        self.assertEqual(pl.report_layout.report_regions,
                         None,
                         "MutRepController cannot correctly read report layout info 'report regions' from jobs setup file")
#        self.assertEqual(pl.report_layout.freq_ratios["ExAC"],
#                         0.5,
#                         "MutRepController cannot correctly read report layout info 'frequency ratios' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.patterns["expr_with_key"],
                         '\'"abc" < 4\'',
                         "MutRepController cannot correctly read report layout info 'expressions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.patterns["expr_wo"],
                         'jkl>5',
                         "MutRepController cannot correctly read report layout info 'expressions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_DELETE_ROW][0].pattern,
                         '\'"abc" < 4\'',
                         "MutRepController cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_COL]['F_U_8'][0].pattern,
                         'jkl>5',
                         "MutRepController cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_COL]['F_U_8'][0].col_name,
                         'F_U_8',
                         "MutRepController cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_COL]['F_U_8'][0].color,
                         'red',
                         "MutRepController cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_ROW][0].pattern,
                         'jkl>5',
                         "MutRepController cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_ROW][0].color,
                         'green',
                         "MutRepController cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_ROW][1].pattern,
                         '2>3',
                         "MutRepController cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.exprs.actions[ACTION_COLOR_ROW][1].color,
                         'orange',
                         "MutRepController cannot correctly read report layout info 'expression actions' from jobs setup file")
        self.assertEqual(pl.report_layout.filter_genes[1],
                         'NOTCH1',
                         "MutRepController cannot correctly read report layout info 'filter genes' from jobs setup file")
        self.assertEqual(pl.report_layout.filter_genes[2],
                         'NOTCH4',
                         "MutRepController cannot correctly read report layout info 'filter genes' from jobs setup file")
        self.assertEqual(pl.report_layout.color_genes[0],
                         'BRCA1',
                         "MutRepController cannot correctly read report layout info 'color genes' from jobs setup file")
        self.assertEqual(pl.report_layout.color_genes[2],
                         'SRC',
                         "MutRepController cannot correctly read report layout info 'color genes' from jobs setup file")

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_select_dataset_1(self):
        """ test selection only rows with samples information """

        self.init_test(self.current_func_name)
        anno_cols = FUNC_REFGENE_COL_NAME
        anno_cols += "," + EXONICFUNC_REFGENE_COL_NAME
        anno_cols += "," + GENE_REFGENE_COL_NAME
        anno_cols += "," + GENEDETAIL_REFGENE_COL_NAME
        anno_cols += "," + SWEGEN_AF_COL_NAME
        anno_cols += "," + WES294_OAF_BRCS_AF_COL_NAME
        jobs_setup_file = self.__create_jobs_setup_file(anno_cols=anno_cols,
                                                        sample_info="8:Co-35:Co-37,13:Co-95,275:Co-1262:Co-618,296:Co-793:Co-876",
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report()
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         42,
                         "Incorrect number of rows"
                         )

        jobs_setup_file = self.__create_jobs_setup_file(anno_cols=anno_cols,
                                                        sample_info="119:Co-309",
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report()
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         2,
                         "Incorrect number of rows"
                         )

        jobs_setup_file = self.__create_jobs_setup_file(anno_cols=anno_cols,
                                                        sample_info="119:Co-309,Co-1209,2016-00043",
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report()
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         49,
                         "Incorrect number of rows"
                         )

        jobs_setup_file = self.__create_jobs_setup_file(anno_cols=anno_cols,
                                                        sample_info="119:Co-309,2016-00043",
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report()
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         26,
                         "Incorrect number of rows"
                         )

        jobs_setup_file = self.__create_jobs_setup_file(anno_cols=anno_cols,
                                                        sample_info="2016-00043",
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report()
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         25,
                         "Incorrect number of rows"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_multiple_regions_1(self):
        """ test with multiple report_regions """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(report_regions="10:89683800-89685413,11")
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report()
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         14,
                         "Incorrect number of rows")

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_unicode_1(self):
        """ test if unicode character 'ö' is allowed in the report """

        self.init_test(self.current_func_name)
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_NON_INTRONIC
        jobs_setup_file = self.__create_jobs_setup_file(anno_excl_tags=DFLT_TEST_ANNO_EXCL_TAGS,
                                                        rows_filter_actions=rows_filter_actions,
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="17")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(),
                         2,
                         "report with swedish unicode character cannot be generated correctly"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_expression_action_del_row_1(self):
        """
        test if expreesion patterns and expression actions can be used
        for deleting rows with a given pattern
        """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(expression_patterns='exp_ns_snv:"ExonicFunc.refGene"==\'synonymous_SNV\'',
                                                        expression_usages="exp_ns_snv:DELETE_ROW",
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="6:78171940-78172992,18:28610987-28611790")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         9,
                         "Incorrect number of rows"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_expression_action_del_row_2(self):
        """
        test if expreesion patterns and expression actions can be used
        for deleting rows by numeric condition(s)
        """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(expression_patterns='exp_OAF1:"AXEQ_CHR3_6_14_18_PF"==1.0000',
                                                        expression_usages="exp_OAF1:DELETE_ROW",
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="6:78171940-78172992,18:28610987-28611790")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         10,
                         "Incorrect number of rows"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_expression_action_del_row_3(self):
        """
        test if expreesion patterns and expression actions can be used
        for deleting rows with multiallelic
        """

        self.init_test(self.current_func_name)
        anno_cols = list(DFLT_TEST_MUTREP_COLS)
        anno_cols.append("OAF_EARLYONSET_AF")
        jobs_setup_file = self.__create_jobs_setup_file(report_regions="22",
                                                        expression_patterns='exp_OAF0:"OAF_EARLYONSET_AF"==0.0000,exp_OAF1:"OAF_EARLYONSET_AF"==1.0000,exp_OAF_empty:"OAF_EARLYONSET_AF"==\'\',exp_OAF_NA:"OAF_EARLYONSET_AF"==\'NA\'',
                                                        expression_usages="exp_OAF0:DELETE_ROW,exp_OAF1:DELETE_ROW,exp_OAF_empty:DELETE_ROW,exp_OAF_NA:DELETE_ROW",
                                                        anno_cols=anno_cols,
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(mc.report_layout.report_regions)
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         4,
                         "Incorrect number of rows"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_expression_action_del_row_4(self):
        """
        test if expreesion patterns and expression actions can be used
        for deleting rows with "true synonymous" variants
        """

        self.init_test(self.current_func_name)
        anno_cols = list(DFLT_TEST_MUTREP_COLS)
        anno_cols.append("OAF_EARLYONSET_AF")
        anno_cols.append("dpsi_zscore")
        jobs_setup_file = self.__create_jobs_setup_file(expression_patterns='true_synonymous:("dpsi_zscore"!=\'\')and(float("dpsi_zscore")>-2)and("ExonicFunc.refGene"==\'synonymous_SNV\')and(float("dpsi_zscore")<2)',
                                                        expression_usages="true_synonymous:DELETE_ROW",
                                                        anno_cols=anno_cols,
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="1")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         2,
                         "Incorrect number of rows"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_expression_action_del_row_5(self):
        """
        test deleting rows using max_ref_maf column
        """

        self.init_test(self.current_func_name)
        anno_cols = list(DFLT_TEST_MUTREP_COLS)
        anno_cols.append(MAX_REF_MAF_COL_NAME)
        anno_cols += REF_MAF_COL_NAMES
        jobs_setup_file = self.__create_jobs_setup_file(expression_patterns='max_ref_maf_10percent:float("MAX_REF_MAF")>0.1',
                                                        expression_usages="max_ref_maf_10percent:DELETE_ROW",
                                                        anno_cols=anno_cols,
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report()
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         26,
                         "Incorrect number of rows"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_expression_action_color_row_1(self):
        """
        test if expreesion patterns and expression actions can be used
        for coloring rows with a given pattern
        """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(report_regions="6:78171940-78172992,18:28610987-28611790",
                                                        expression_patterns='exp_splicing:"Func.refGene"==\'splicing\'',
                                                        expression_usages="exp_splicing:COLOR_ROW:ROSY_BROWN",
#                                                        call_detail="YES",
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(mc.report_layout.report_regions)
        xu = XlsUtils(xls_file)
        col_idx = xu.get_col_idx(col_name="Func.refGene", sheet_idx=0)
        self.assertEqual(xu.get_cell_value(8, col_idx, sheet_idx=0),
                         "splicing",
                         "Incorrect cell value"
                         )
        exp_rgb = "FF" + COLORS_RGB["ROSY_BROWN"][-6:]
        self.assertEqual(xu.get_cell_rgb(8, 2, sheet_idx=0),
                         exp_rgb,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_rgb(8, col_idx, sheet_idx=0),
                         exp_rgb,
                         "Incorrect color"
                         )
        self.assertNotEqual(xu.get_cell_rgb(7, col_idx, sheet_idx=0),
                            exp_rgb,
                            "Incorrect color"
                            )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_expression_action_color_col_1(self):
        """
        test if expreesion patterns and expression actions can be used
        for coloring rows with a given pattern
        """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(expression_patterns='exp_exonic:"Func.refGene"==\'exonic\',exp_ns_snv:"ExonicFunc.refGene"==\'nonsynonymous_SNV\'',
                                                        expression_usages="exp_exonic:COLOR_ROW:ROSY_BROWN,exp_ns_snv:COLOR_COLUMN:cytoBand:TEAL,exp_exonic:COLOR_COLUMN:Gene.refGene:YELLOW",
#                                                        call_detail="YES",
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="6:78171940-78172992,18:28610987-28611790")
        xu = XlsUtils(xls_file)
        col_idx = xu.get_col_idx(col_name="Func.refGene", sheet_idx=0)
        self.assertEqual(xu.get_cell_value(4, col_idx, sheet_idx=0),
                         "exonic",
                         "Incorrect cell value"
                         )
        exp_rgb = "FF" + COLORS_RGB["ROSY_BROWN"][-6:]
        self.assertEqual(xu.get_cell_rgb(5, col_idx, sheet_idx=0),
                         exp_rgb,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_rgb(5, col_idx, sheet_idx=0),
                         exp_rgb,
                         "Incorrect color"
                         )
        self.assertNotEqual(xu.get_cell_rgb(2, col_idx, sheet_idx=0),
                            exp_rgb,
                            "Incorrect color"
                            )
        col_idx = xu.get_col_idx(col_name="cytoBand", sheet_idx=0)
        exp_rgb = "FF" + COLORS_RGB["TEAL"][-6:]
        self.assertEqual(xu.get_cell_rgb(4, col_idx, sheet_idx=0),
                         exp_rgb,
                         "Incorrect color"
                         )
        col_idx = xu.get_col_idx(col_name="Gene.refGene", sheet_idx=0)
        exp_rgb = "FF" + COLORS_RGB["YELLOW"][-6:]
        self.assertEqual(xu.get_cell_rgb(4, col_idx, sheet_idx=0),
                         exp_rgb,
                         "Incorrect color"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_expression_action_color_col_2(self):
        """
        test multiple coloring in one column with different expression
        """

        self.init_test(self.current_func_name)
        anno_cols = list(DFLT_TEST_MUTREP_COLS)
        anno_cols.append(SWEGEN_AF_COL_NAME)
        expression_patterns='swegen_rare:("SWEGEN_AF"!=\'\')'
        expression_patterns+='and(float("SWEGEN_AF")<0.03)'
        expression_usages='swegen_rare:COLOR_COLUMN:SWEGEN_AF:YELLOW'
        expression_patterns+=',swegen_very_rare:("SWEGEN_AF"!=\'\')'
        expression_patterns+='and(float("SWEGEN_AF")<0.01)'
        expression_usages+=',swegen_very_rare:COLOR_COLUMN:SWEGEN_AF:XLS_ORANGE'
        expression_patterns+=',swegen_empty:("SWEGEN_AF"==\'\')'
        expression_usages+=',swegen_empty:COLOR_COLUMN:SWEGEN_AF:XLS_ORANGE'
        expression_patterns+=',kg_all_rare:("1000g2014oct_all"!=\'\')'
        expression_patterns+='and(float("1000g2014oct_all")<0.03)'
        expression_usages+=',kg_all_rare:COLOR_COLUMN:1000g2014oct_all:YELLOW'
        jobs_setup_file = self.__create_jobs_setup_file(expression_patterns=expression_patterns,
                                                        expression_usages=expression_usages,
                                                        anno_cols=anno_cols,
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report()
        xu = XlsUtils(xls_file)
        col_idx = xu.get_col_idx(col_name="SWEGEN_AF", sheet_idx=0)
        exp_yellow = "FF" + COLORS_RGB["YELLOW"][-6:]
        exp_orange = "FFFDB50A"
        self.assertEqual(xu.get_cell_rgb(2, col_idx, sheet_idx=0),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_rgb(3, col_idx, sheet_idx=0),
                         exp_yellow,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_rgb(4, col_idx, sheet_idx=0),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_rgb(5, col_idx, sheet_idx=0),
                         exp_orange,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_rgb(6, col_idx, sheet_idx=0),
                         exp_orange,
                         "Incorrect color"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_coloring_shared_1(self):
        """ test coloring shared samples when the option is on """

        self.init_test(self.current_func_name)
        sample_info = join_path(self.data_dir,
                                "sample.info")
        jobs_setup_file = self.__create_jobs_setup_file(sample_info=sample_info,
                                                        coloring_shared=True,
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(mc.report_layout.report_regions)
        xu = XlsUtils(xls_file)
        sample_col_idx = xu.get_col_idx("256-Co-388")
        exp_rgb = "FF" + COLORS_RGB["LIGHT_BLUE"][-6:]
        self.assertEqual(xu.get_cell_value(2, sample_col_idx),
                         "wt",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(2, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_value(3, sample_col_idx),
                         "het",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(3, sample_col_idx),
                         exp_rgb,
                         "Incorrect color"
                         )
        sample_col_idx = xu.get_col_idx("13-Co-95")
        exp_rgb = "FF" + COLORS_RGB["ICEBLUE"][-6:]
        self.assertEqual(xu.get_cell_value(2, sample_col_idx),
                         ".",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(2, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_value(3, sample_col_idx),
                         "hom",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(3, sample_col_idx),
                         exp_rgb,
                         "Incorrect color"
                         )
        sample_col_idx = xu.get_col_idx("110-Co-1313")
        self.assertEqual(xu.get_cell_value(2, sample_col_idx),
                         "wt",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(2, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_value(17, sample_col_idx),
                         "het",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(17, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        sample_col_idx = xu.get_col_idx("771-08F")
        self.assertEqual(xu.get_cell_value(2, sample_col_idx),
                         ".",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(2, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_value(3, sample_col_idx),
                         "hom",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(3, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_value(15, sample_col_idx),
                         "het",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(15, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )

#    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
#    def test_cal_est_kvot_1(self):
#        """
#        test variable type handling in test_cal_ors
#        """
#
#        self.init_test(self.current_func_name)
#        db_file = join_path(self.data_dir,
#                                        "input.vcf.gz")
#        project_name = self.test_function
#        anno_cols = list(DFLT_TEST_MUTREP_COLS)
#        anno_cols.append(KG2014OCT_EUR_COL_NAME)
#        anno_cols.append(EST_KVOT_EARLYONSET_VS_BRC_COL_NAME)
#        anno_cols.append(WES294_OAF_EARLYONSET_AF_COL_NAME)
#        anno_cols.append(WES294_OAF_BRCS_AF_COL_NAME)
#        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
#                                                        db_file=db_file,
#                                                        report_regions="22",
#                                                        anno_cols=anno_cols,
#                                                        )
#        pl = MutRepController(jobs_setup_file=jobs_setup_file)
#        pl.gen_summary_report(pl.report_layout.report_regions)
#        xls_file = join_path(self.working_dir,
#                             "rpts",
#                             project_name+"_summary.xlsx")
#        xu = XlsUtils(xls_file)
#        self.assertEqual(xu.count_rows(sheet_idx=0),
#                         23,
#                         "Incorrect number of rows"
#                         )
#
#    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
#    def test_cal_est_kvot_2(self):
#        """
#        test if there is mutation in reference for cal_est_ors
#        """
#
#        self.init_test(self.current_func_name)
#        db_file = join_path(self.data_dir,
#                                        "input.vcf.gz")
#        sample_info = join_path(self.data_dir,
#                                "sample.info")
#        project_name = self.test_function
#        anno_cols = list(DFLT_TEST_MUTREP_COLS)
#        anno_cols.append(KG2014OCT_EUR_COL_NAME)
#        anno_cols.append(EST_KVOT_EARLYONSET_VS_BRC_COL_NAME)
#        anno_cols.append(EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME)
#        anno_cols.append(EST_KVOT_EARLYONSET_VS_KG_EUR_COL_NAME)
#        anno_cols.append(WES294_OAF_EARLYONSET_AF_COL_NAME)
#        anno_cols.append(WES294_OAF_EARLYONSET_GF_COL_NAME)
#        anno_cols.append(WES294_OAF_BRCS_AF_COL_NAME)
#        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
#                                                        db_file=db_file,
#                                                        report_regions="22",
#                                                        sample_info=sample_info,
#                                                        anno_cols=anno_cols,
#                                                        )
#        pl = MutRepController(jobs_setup_file=jobs_setup_file)
#        pl.gen_summary_report(pl.report_layout.report_regions)
#        xls_file = join_path(self.working_dir,
#                             "rpts",
#                             project_name+"_summary.xlsx")
#        xu = XlsUtils(xls_file)
#        ors_col_idx = xu.get_col_idx(EST_KVOT_EARLYONSET_VS_BRC_COL_NAME)
#        self.assertEqual(xu.get_cell_value(4, ors_col_idx),
#                         "NA",
#                         "Incorect ORS estimation"
#                         )
#        self.assertEqual(xu.get_cell_value(5, ors_col_idx),
#                         1.5273,
#                         "Incorect ORS estimation"
#                         )
#        self.assertEqual(xu.get_cell_value(6, ors_col_idx),
#                         "INF",
#                         "Incorect ORS estimation"
#                         )
#        ors_col_idx = xu.get_col_idx(EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME)
#        self.assertEqual(xu.get_cell_value(4, ors_col_idx),
#                         "NA",
#                         "Incorect ORS estimation"
#                         )
#        self.assertEqual(xu.get_cell_value(5, ors_col_idx),
#                         1.3446,
#                         "Incorect ORS estimation"
#                         )
#        self.assertEqual(xu.get_cell_value(6, ors_col_idx),
#                         0.5681,
#                         "Incorect ORS estimation"
#                         )
#        ors_col_idx = xu.get_col_idx(EST_KVOT_EARLYONSET_VS_KG_EUR_COL_NAME)
#        self.assertEqual(xu.get_cell_value(4, ors_col_idx),
#                         "NA",
#                         "Incorect ORS estimation"
#                         )
#        self.assertEqual(xu.get_cell_value(5, ors_col_idx),
#                         1.2375,
#                         "Incorect ORS estimation"
#                         )
#        self.assertEqual(xu.get_cell_value(6, ors_col_idx),
#                         0.9378,
#                         "Incorect ORS estimation"
#                         )
#
#    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
#    def test_criteria_sheet_1(self):
#        """ test summary with multiple report_regions and many sample infos """
#
#        self.init_test(self.current_func_name)
#        db_file = join_path(self.data_dir,
#                                        "chr6_18.vcf.gz")
#        project_name = self.test_function
#        rows_filter_actions = JOBS_SETUP_RPT_FILTER_NON_INTRONIC
#        jobs_setup_file = self.__create_jobs_setup_file(project_name=project_name,
#                                                        db_file=db_file,
#                                                        report_regions="6:78171940-78172992,18:28610987-28611790",
#                                                        sample_info="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
#                                                        call_detail="YES",
#                                                        rows_filter_actions=rows_filter_actions,
#                                                        )
#        pl = MutRepController(jobs_setup_file=jobs_setup_file)
#        pl.gen_summary_report(pl.report_layout.report_regions)
#        xls_file = join_path(self.working_dir,
#                             "rpts",
#                             project_name+"_summary.xlsx")
#        xu = XlsUtils(xls_file)
#        self.assertEqual(xu.count_rows(sheet_idx=0),
#                         8,
#                         "Incorrect number of rows in the variants sheet"
#                         )
#        self.assertEqual(xu.count_cols(col_name1="1234-Alb-31",
#                                       col_name2="6789-Al-65",
#                                       sheet_idx=0),
#                         11,
#                         "Incorrect number of columns"
#                         )
#        self.assertEqual(xu.nsheets,
#                         2,
#                         "Incorrect number of sheets"
#                         )
#        self.assertEqual(xu.count_rows(sheet_idx=1),
#                         6,
#                         "Incorrect number of rows the criteria sheet"
#                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_header_corrections_1(self):
        """ test if column headers can be replaced with better names  """

        self.init_test(self.current_func_name)
        anno_cols = list(DFLT_TEST_MUTREP_COLS)
        anno_cols.append("OAF_EARLYONSET_AF")
        anno_cols.append("OAF_BRC_CRC_PROSTATE_AF")
        header_corrections = "OAF_EARLYONSET_AF:EARLYONSET_AF"
        header_corrections += ",OAF_BRC_CRC_PROSTATE_AF:ALL_EXOME_AF"
        jobs_setup_file = self.__create_jobs_setup_file(header_corrections=header_corrections,
                                                        anno_cols=anno_cols,
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="22")
        xu = XlsUtils(xls_file)
        oaf_col_idx = xu.get_col_idx("ALL_EXOME_AF")
        self.assertEqual(xu.get_cell_value(4, oaf_col_idx),
                         0.1250,
                         "Incorect allele frequency"
                         )
        self.assertEqual(xu.get_cell_value(6, oaf_col_idx),
                         None,
                         "Incorect allele frequency"
                         )
        self.assertEqual(xu.get_cell_value(7, oaf_col_idx),
                         None,
                         "Incorect allele frequency"
                         )
        self.assertEqual(xu.get_cell_value(8, oaf_col_idx),
                         None,
                         "Incorect allele frequency"
                         )
        self.assertEqual(xu.get_cell_value(9, oaf_col_idx),
                         0.0571,
                         "Incorect allele frequency"
                         )
        self.assertEqual(xu.get_cell_value(24, oaf_col_idx),
                         0,
                         "Incorect allele frequency"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_filter_genes_1(self):
        """
        test if gene search can be done correctly
        both inside the genes and surrounding.
        """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(sample_info="8:Co-35:Co-37,13:Co-95,275:Co-1262:Co-618,296:Co-793:Co-876")
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="9")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         42,
                         "Incorrect number of rows in the variants sheet"
                         )

        filter_genes = "ANKRD19P"
        jobs_setup_file = self.__create_jobs_setup_file(sample_info="8:Co-35:Co-37,13:Co-95,275:Co-1262:Co-618,296:Co-793:Co-876",
                                                        filter_genes=filter_genes,
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="9")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         19,
                         "Incorrect number of rows in the variants sheet"
                         )

        filter_genes = "ANKRD19P,IPPK"
        jobs_setup_file = self.__create_jobs_setup_file(sample_info="8:Co-35:Co-37,13:Co-95,275:Co-1262:Co-618,296:Co-793:Co-876",
                                                        filter_genes=filter_genes,
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="9")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         25,
                         "Incorrect number of rows in the variants sheet"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_color_genes_1(self):
        """
        test if basic gene search (one gene, full name) can be done correctly
        both inside the genes and surrounding.
        """

        self.init_test(self.current_func_name)
        color_genes = "ANKRD19P"
        jobs_setup_file = self.__create_jobs_setup_file(sample_info="8:Co-35:Co-37,13:Co-95,275:Co-1262:Co-618,296:Co-793:Co-876",
                                                        color_genes=color_genes,
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="9")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         42,
                         "Incorrect number of rows"
                         )
        refgene_col_idx = xu.get_col_idx(GENE_REFGENE_COL_NAME)
        exp_rgb = "FF" + COLORS_RGB["XLS_GREEN"][-6:]
        self.assertEqual(xu.get_cell_rgb(2, refgene_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_rgb(20, refgene_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_rgb(30, refgene_col_idx),
                         exp_rgb,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_rgb(40, refgene_col_idx),
                         exp_rgb,
                         "Incorrect color"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_color_genes_2(self):
        """
        test if a little advance gene search (two genes, full name) can be done correctly
        both inside the genes and surrounding.
        """

        self.init_test(self.current_func_name)
        color_genes = "ANKRD19P,IPPK"
        jobs_setup_file = self.__create_jobs_setup_file(sample_info="8:Co-35:Co-37,13:Co-95,275:Co-1262:Co-618,296:Co-793:Co-876",
                                                        color_genes=color_genes,
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="9")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         42,
                         "Incorrect number of rows"
                         )
        refgene_col_idx = xu.get_col_idx(GENE_REFGENE_COL_NAME)
        exp_rgb = "FF" + COLORS_RGB["XLS_GREEN"][-6:]
        self.assertEqual(xu.get_cell_rgb(2, refgene_col_idx),
                         exp_rgb,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_rgb(20, refgene_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_rgb(30, refgene_col_idx),
                         exp_rgb,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_rgb(40, refgene_col_idx),
                         exp_rgb,
                         "Incorrect color"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_coloring_zygosity_1(self):
        """ test coloring shared samples when the option is off """

        self.init_test(self.current_func_name)
        custom_excl_tags = DFLT_TEST_ANNO_EXCL_TAGS
        custom_excl_tags += "," + AXEQ_CHR3_6_14_18_COLS_TAG
        custom_excl_tags += "," + AXEQ_CHR5_19_COLS_TAG
        custom_excl_tags += "," + LJB_SCORE_COLS_TAG
        sample_info = join_path(self.data_dir,
                                "sample.info")
        jobs_setup_file = self.__create_jobs_setup_file(anno_cols=ALL_MUTREP_ANNO_COLS,
                                                        anno_excl_tags=custom_excl_tags,
                                                        sample_info=sample_info,
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="9")
        xu = XlsUtils(xls_file)
        sample_col_idx = xu.get_col_idx("256-Co-388")
        self.assertEqual(xu.get_cell_value(2, sample_col_idx),
                         "wt",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(2, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_value(3, sample_col_idx),
                         "het",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(3, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        sample_col_idx = xu.get_col_idx("13-Co-95")
        self.assertEqual(xu.get_cell_value(2, sample_col_idx),
                         ".",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(2, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_value(3, sample_col_idx),
                         "hom",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(3, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        sample_col_idx = xu.get_col_idx("110-Co-1313")
        self.assertEqual(xu.get_cell_value(2, sample_col_idx),
                         "wt",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(2, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_value(18, sample_col_idx),
                         "het",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(18, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        sample_col_idx = xu.get_col_idx("771-08F")
        self.assertEqual(xu.get_cell_value(2, sample_col_idx),
                         ".",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(2, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_value(3, sample_col_idx),
                         "hom",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(3, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_value(26, sample_col_idx),
                         "het",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(26, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_coloring_zygosity_2(self):
        """ test coloring shared samples when the option is on """

        self.init_test(self.current_func_name)
        custom_excl_tags = DFLT_TEST_ANNO_EXCL_TAGS
        custom_excl_tags += "," + AXEQ_CHR3_6_14_18_COLS_TAG
        custom_excl_tags += "," + AXEQ_CHR5_19_COLS_TAG
        custom_excl_tags += "," + LJB_SCORE_COLS_TAG
        sample_info = join_path(self.data_dir,
                                "sample.info")
        jobs_setup_file = self.__create_jobs_setup_file(anno_cols=ALL_MUTREP_ANNO_COLS,
                                                        anno_excl_tags=custom_excl_tags,
                                                        sample_info=sample_info,
                                                        coloring_zygosity=True,
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="9")
        xu = XlsUtils(xls_file)
        sample_col_idx = xu.get_col_idx("256-Co-388")
        exp_rgb = "FF" + COLORS_RGB["XLS_LIGHT_GREEN"][-6:]
        self.assertEqual(xu.get_cell_value(2, sample_col_idx),
                         "wt",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(2, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_value(3, sample_col_idx),
                         "het",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(3, sample_col_idx),
                         exp_rgb,
                         "Incorrect color"
                         )
        sample_col_idx = xu.get_col_idx("13-Co-95")
        exp_rgb = "FF" + COLORS_RGB["XLS_DARK_GREEN"][-6:]
        self.assertEqual(xu.get_cell_value(2, sample_col_idx),
                         ".",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(2, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_value(3, sample_col_idx),
                         "hom",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(3, sample_col_idx),
                         exp_rgb,
                         "Incorrect color"
                         )
        sample_col_idx = xu.get_col_idx("1207-2818-07D")
        exp_rgb = "FF" + COLORS_RGB["XLS_LIGHT_GREEN"][-6:]
        self.assertEqual(xu.get_cell_value(25, sample_col_idx),
                         "het",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(25, sample_col_idx),
                         exp_rgb,
                         "Incorrect color"
                         )
        sample_col_idx = xu.get_col_idx("771-08F")
        self.assertEqual(xu.get_cell_value(2, sample_col_idx),
                         ".",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(2, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )
        exp_rgb = "FF" + COLORS_RGB["XLS_DARK_GREEN"][-6:]
        self.assertEqual(xu.get_cell_value(3, sample_col_idx),
                         "hom",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(3, sample_col_idx),
                         exp_rgb,
                         "Incorrect color"
                         )
        exp_rgb = "FF" + COLORS_RGB["XLS_LIGHT_GREEN"][-6:]
        self.assertEqual(xu.get_cell_value(11, sample_col_idx),
                         "het",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(11, sample_col_idx),
                         exp_rgb,
                         "Incorrect color"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_datasets_1(self):
        """
        test if samples in mutation report can be displayed in groups
          - no group
        """

        self.init_test(self.current_func_name)
        custom_excl_tags = DFLT_TEST_ANNO_EXCL_TAGS
        custom_excl_tags += "," + AXEQ_CHR3_6_14_18_COLS_TAG
        custom_excl_tags += "," + AXEQ_CHR5_19_COLS_TAG
        custom_excl_tags += "," + LJB_SCORE_COLS_TAG
#        rows_filter_actions = JOBS_SETUP_RPT_FILTER_HAS_MUTATION
#        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_PASS_VQSR
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_NON_INTERGENIC
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_INTRONIC
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_UPSTREAM
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_UTR
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS
        expression_patterns = OrderedDict()
        expression_patterns['filter_NA_young'] = '"' + WES294_OAF_EARLYONSET_AF_COL_NAME + '"==\'NA\''
        expression_patterns['filter_empty_young'] = '"' + WES294_OAF_EARLYONSET_AF_COL_NAME + '"==\'\''
        expression_patterns['filter_unique_young'] = 'float("' + WES294_OAF_EARLYONSET_AF_COL_NAME + '")<0.0196'
#        expression_patterns['filter_kvot_young_vs_exac'] = '("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '"!=\'INF\')'
#        expression_patterns['filter_kvot_young_vs_exac'] += 'and("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '"!=\'\')'
#        expression_patterns['filter_kvot_young_vs_exac'] += 'and("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '"!=\'NA\')'
#        expression_patterns['filter_kvot_young_vs_exac'] += 'and(float("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '")<1.3)'
#        expression_patterns['filter_kvot_young_vs_brc'] = '("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '"!=\'INF\')'
#        expression_patterns['filter_kvot_young_vs_brc'] += 'and("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '"!=\'\')'
#        expression_patterns['filter_kvot_young_vs_brc'] += 'and(float("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '")<1.3)'
        expression_patterns['filter_ncRNA_exonic'] = '"' + FUNC_REFGENE_COL_NAME + '"==\'ncRNA_exonic\''
        sample_info = join_path(self.data_dir,
                                "sample.info")
        jobs_setup_file = self.__create_jobs_setup_file(anno_cols=ALL_MUTREP_ANNO_COLS,
                                                        anno_excl_tags=custom_excl_tags,
                                                        sample_info=sample_info,
                                                        rows_filter_actions=rows_filter_actions,
                                                        expression_patterns=",".join(map(lambda x: x+":"+expression_patterns[x], expression_patterns)),
                                                        expression_usages=",".join(map(lambda x: x+":DELETE_ROW", expression_patterns)),
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="2,6,7,19")
        xu = XlsUtils(xls_file)
        last_anno_col_idx = xu.get_col_idx("CADD_phred")
        first_sample_col_idx = xu.get_col_idx("1003-06o")
        self.assertEqual(first_sample_col_idx-last_anno_col_idx,
                         2,
                         "Incorrect report layout"
                         )
        self.assertEqual(xu.get_cell_value(5, first_sample_col_idx),
                         "het",
                         "Incorrect cell value"
                         )
        last_sample_col_idx = xu.get_col_idx("L862-Br-716")
        self.assertEqual(last_sample_col_idx-first_sample_col_idx,
                         290,
                         "Incorrect report layout"
                         )
        self.assertEqual(xu.get_cell_value(5, last_sample_col_idx),
                         "het",
                         "Incorrect cell value"
                         )
        sample_col_idx = xu.get_col_idx("L862-Br-710")
        self.assertEqual(xu.get_cell_value(6, sample_col_idx),
                         ".",
                         "Incorrect cell value"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_datasets_2(self):
        """
        test if samples in mutation report can be displayed in groups
          - all the samples have group information
        """

        self.init_test(self.current_func_name)
        custom_excl_tags = DFLT_TEST_ANNO_EXCL_TAGS
        custom_excl_tags += "," + AXEQ_CHR3_6_14_18_COLS_TAG
        custom_excl_tags += "," + AXEQ_CHR5_19_COLS_TAG
        custom_excl_tags += "," + LJB_SCORE_COLS_TAG
#        rows_filter_actions = JOBS_SETUP_RPT_FILTER_HAS_MUTATION
#        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_PASS_VQSR
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_NON_INTERGENIC
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_INTRONIC
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_UPSTREAM
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_UTR
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS
        expression_patterns = OrderedDict()
        expression_patterns['filter_NA_young'] = '"' + WES294_OAF_EARLYONSET_AF_COL_NAME + '"==\'NA\''
        expression_patterns['filter_empty_young'] = '"' + WES294_OAF_EARLYONSET_AF_COL_NAME + '"==\'\''
        expression_patterns['filter_unique_young'] = 'float("' + WES294_OAF_EARLYONSET_AF_COL_NAME + '")<0.0196'
#        expression_patterns['filter_kvot_young_vs_exac'] = '("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '"!=\'INF\')'
#        expression_patterns['filter_kvot_young_vs_exac'] += 'and("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '"!=\'\')'
#        expression_patterns['filter_kvot_young_vs_exac'] += 'and("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '"!=\'NA\')'
#        expression_patterns['filter_kvot_young_vs_exac'] += 'and(float("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '")<1.3)'
#        expression_patterns['filter_kvot_young_vs_brc'] = '("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '"!=\'INF\')'
#        expression_patterns['filter_kvot_young_vs_brc'] += 'and("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '"!=\'\')'
#        expression_patterns['filter_kvot_young_vs_brc'] += 'and(float("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '")<1.3)'
        expression_patterns['filter_ncRNA_exonic'] = '"' + FUNC_REFGENE_COL_NAME + '"==\'ncRNA_exonic\''
        sample_info = join_path(self.data_dir,
                                "sample.info")
        jobs_setup_file = self.__create_jobs_setup_file(anno_cols=ALL_MUTREP_ANNO_COLS,
                                                        anno_excl_tags=custom_excl_tags,
                                                        sample_info=sample_info,
                                                        rows_filter_actions=rows_filter_actions,
                                                        expression_patterns=",".join(map(lambda x: x+":"+expression_patterns[x], expression_patterns)),
                                                        expression_usages=",".join(map(lambda x: x+":DELETE_ROW", expression_patterns)),
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="2,6,7,19")
        xu = XlsUtils(xls_file)
        last_anno_col_idx = xu.get_col_idx("CADD_phred")
        first_sample_col_idx = xu.get_col_idx("1003-06o")
        self.assertEqual(first_sample_col_idx-last_anno_col_idx,
                         2,
                         "Incorrect report layout"
                         )
        self.assertEqual(xu.get_cell_value(5, first_sample_col_idx),
                         "het",
                         "Incorrect cell value"
                         )
        last_sample_col_idx = xu.get_col_idx("L862-Br-716")
        self.assertEqual(last_sample_col_idx-first_sample_col_idx,
                         297,
                         "Incorrect report layout"
                         )
        self.assertEqual(xu.get_cell_value(5, last_sample_col_idx),
                         "het",
                         "Incorrect cell value"
                         )
        sample_col_idx = xu.get_col_idx("L862-Br-710")
        self.assertEqual(xu.get_cell_value(6, sample_col_idx),
                         ".",
                         "Incorrect cell value"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_show_shared_variants_1(self):
        """ test if shared mutations can be shown """

        self.init_test(self.current_func_name)
        custom_excl_tags = DFLT_TEST_ANNO_EXCL_TAGS
        custom_excl_tags += "," + AXEQ_CHR3_6_14_18_COLS_TAG
        custom_excl_tags += "," + AXEQ_CHR5_19_COLS_TAG
        custom_excl_tags += "," + LJB_SCORE_COLS_TAG
        sample_info = join_path(self.data_dir,
                                "sample.info")
        jobs_setup_file = self.__create_jobs_setup_file(anno_excl_tags=custom_excl_tags,
                                                        sample_info=sample_info,
                                                        show_shared_variants=True,
                                                        coloring_shared=True,
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="3")
        xu = XlsUtils(xls_file)
        fam_col_idx = xu.get_col_idx("F0007826")
        self.assertEqual(xu.get_cell_value(2, fam_col_idx),
                         "not shared",
                         "Incorrect cell value"
                         )
        exp_rgb = "FF" + COLORS_RGB["XLS_GRAY1_2"][-6:]
        self.assertEqual(xu.get_cell_rgb(2, fam_col_idx),
                         exp_rgb,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_value(3, fam_col_idx),
                         "shared",
                         "Incorrect cell value"
                         )
        exp_rgb = "FF" + COLORS_RGB["LIGHT_BLUE"][-6:]
        self.assertEqual(xu.get_cell_rgb(3, fam_col_idx),
                         exp_rgb,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_value(10, fam_col_idx),
                         "not shared",
                         "Incorrect cell value"
                         )
        exp_rgb = "FF" + COLORS_RGB["XLS_GRAY1_2"][-6:]
        self.assertEqual(xu.get_cell_rgb(10, fam_col_idx),
                         exp_rgb,
                         "Incorrect color"
                         )

    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_filter_non_recessive_gene_1(self):
        """ test basic filtering non-recessive gene """

        self.init_test(self.current_func_name)
        custom_excl_tags = DFLT_TEST_ANNO_EXCL_TAGS
        custom_excl_tags += "," + AXEQ_CHR3_6_14_18_COLS_TAG
        custom_excl_tags += "," + AXEQ_CHR5_19_COLS_TAG
        custom_excl_tags += "," + LJB_SCORE_COLS_TAG
#        rows_filter_actions = JOBS_SETUP_RPT_FILTER_PASS_VQSR
#        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_INTERGENIC
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_NON_INTRONIC
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_UPSTREAM
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_UTR
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_RECESSIVE_GENE
        expression_patterns = OrderedDict()
        expression_patterns['filter_NA_young'] = '"' + WES294_OAF_EARLYONSET_AF_COL_NAME + '"==\'NA\''
        expression_patterns['filter_empty_young'] = '"' + WES294_OAF_EARLYONSET_AF_COL_NAME + '"==\'\''
        expression_patterns['filter_unique_young'] = 'float("' + WES294_OAF_EARLYONSET_AF_COL_NAME + '")<0.0196'
#        expression_patterns['filter_kvot_young_vs_exac'] = '("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '"!=\'INF\')'
#        expression_patterns['filter_kvot_young_vs_exac'] += 'and("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '"!=\'\')'
#        expression_patterns['filter_kvot_young_vs_exac'] += 'and("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '"!=\'NA\')'
#        expression_patterns['filter_kvot_young_vs_exac'] += 'and(float("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '")<1.3)'
#        expression_patterns['filter_kvot_young_vs_brc'] = '("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '"!=\'INF\')'
#        expression_patterns['filter_kvot_young_vs_brc'] += 'and("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '"!=\'\')'
#        expression_patterns['filter_kvot_young_vs_brc'] += 'and(float("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '")<1.3)'
        expression_patterns['filter_ncRNA_exonic'] = '"' + FUNC_REFGENE_COL_NAME + '"==\'ncRNA_exonic\''
        sample_info = join_path(self.data_dir,
                                "sample.info")
        jobs_setup_file = self.__create_jobs_setup_file(anno_excl_tags=custom_excl_tags,
                                                        sample_info=sample_info,
                                                        rows_filter_actions=rows_filter_actions,
                                                        expression_patterns=",".join(map(lambda x: x+":"+expression_patterns[x], expression_patterns)),
                                                        expression_usages=",".join(map(lambda x: x+":DELETE_ROW", expression_patterns)),
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="2,6,7,19")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         11,
                         "Incorrect number of rows in the variants sheet"
                         )
        info_col_idx = xu.get_col_idx(COMPOUND_HETEROZYGOTE_AFFECTED_COUNT_COL_NAME)
        self.assertEqual(xu.get_cell_value(2, info_col_idx),
                         46,
                         "Incorrect cell value",
                         )
        self.assertEqual(xu.get_cell_value(9, info_col_idx),
                         0,
                         "Incorrect cell value",
                         )
        self.assertEqual(xu.get_cell_value(10, info_col_idx),
                         7,
                         "Incorrect cell value",
                         )
        info_col_idx = xu.get_col_idx(COMPOUND_HETEROZYGOTE_FREQ_RATIO_COL_NAME)
        self.assertEqual(xu.get_cell_value(8, info_col_idx),
                         "0.90 vs 0.78",
                         "Incorrect cell value",
                         )
        info_col_idx = xu.get_col_idx(HOMOZYGOTE_AFFECTED_COUNT_COL_NAME)
        self.assertEqual(xu.get_cell_value(10, info_col_idx),
                         0,
                         "Incorrect cell value",
                         )
        info_col_idx = xu.get_col_idx(HOMOZYGOTE_FREQ_RATIO_COL_NAME)
        self.assertEqual(xu.get_cell_value(9, info_col_idx),
                         "0.18 vs 0.05",
                         "Incorrect cell value",
                         )
        sample_col_idx = xu.get_col_idx("1199-05o")
        exp_rgb = "FF" + COLORS_RGB["XLS_GREEN"][-6:]
        self.assertEqual(xu.get_cell_value(7, sample_col_idx),
                         "het",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(7, sample_col_idx),
                         exp_rgb,
                         "Incorrect color"
                         )
        self.assertEqual(xu.get_cell_value(9, sample_col_idx),
                         "het",
                         "Incorrect cell value"
                         )
        self.assertEqual(xu.get_cell_rgb(9, sample_col_idx),
                         RGB_NO_FILL,
                         "Incorrect color"
                         )

#    @unittest.skipUnless(FULL_SYSTEM_TEST or MUTREP_TEST or XLS_TEST or DB_TEST, "taking too long time to test")
    def test_filter_non_recessive_gene_2(self):
        """ test filtering non-recessive gene must include variants found in the controls"""

        self.init_test(self.current_func_name)
        anno_cols = []
        anno_cols.append(FUNC_REFGENE_COL_NAME)
        anno_cols.append(EXONICFUNC_REFGENE_COL_NAME)
        anno_cols.append(GENE_REFGENE_COL_NAME)
        anno_cols.append(GENEDETAIL_REFGENE_COL_NAME)
        anno_cols.append(CYTOBAND_COL_NAME)
        anno_cols.append(MAX_REF_MAF_COL_NAME)
        anno_cols += REF_MAF_COL_NAMES
        anno_cols.append(WES294_OAF_EARLYONSET_AF_COL_NAME)
        anno_cols.append(WES294_OAF_EARLYONSET_GF_COL_NAME)
        anno_cols.append(SPIDEX_DPSI_ZSCORE_COL_NAME)
        anno_cols.append(PATHOGENIC_COUNT_COL_NAME)
        anno_cols += PREDICTION_COLS
        anno_cols.append(LJB_METASVM_PREDICTION_COL_NAME)
        anno_cols.append(LJB_METALR_PREDICTION_COL_NAME)
        rows_filter_actions = JOBS_SETUP_RPT_FILTER_NON_INTRONIC
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_UPSTREAM
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_UTR
        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_RECESSIVE_GENE

        expression_patterns = OrderedDict()
        expression_patterns['true_synonymous'] = '("dpsi_zscore"!=\'\')and("ExonicFunc.refGene"==\'synonymous_SNV\')and(float("dpsi_zscore")>-2)and(float("dpsi_zscore")<2)'
        expression_patterns['likely_synonymous'] = '("dpsi_zscore"==\'\')and("ExonicFunc.refGene"==\'synonymous_SNV\')'
        expression_patterns['MAX_REF_MAF_COMMON'] = '("MAX_REF_MAF"!=\'\')and(float("MAX_REF_MAF")<=0.80)and(float("MAX_REF_MAF")>=0.20)'
        expression_patterns['Pathogenic5'] = 'int("Pathogenic_count")<5'
        expression_patterns['EARLYONSET_lowQC'] = '(float("OAF_EARLYONSET_GF")<0.5)'
        expression_patterns['filter_ncRNA_exonic'] = '"' + FUNC_REFGENE_COL_NAME + '"==\'ncRNA_exonic\''
        sample_info = join_path(self.data_dir,
                                "sample.info")
        jobs_setup_file = self.__create_jobs_setup_file(anno_cols=anno_cols,
                                                        sample_info=sample_info,
                                                        rows_filter_actions=rows_filter_actions,
                                                        expression_patterns=",".join(map(lambda x: x+":"+expression_patterns[x], expression_patterns)),
                                                        expression_usages=",".join(map(lambda x: x+":DELETE_ROW", expression_patterns)),
                                                        )
        mc = MutRepController(jobs_setup_file, verbose=False)
        xls_file = mc.gen_report(report_regions="8")
        xu = XlsUtils(xls_file)
        self.assertEqual(xu.count_rows(sheet_idx=0),
                         9,
                         "Incorrect number of rows in the variants sheet"
                         )
