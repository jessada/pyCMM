# This Python file uses the following encoding: utf-8
import unittest
from os.path import join as join_path
from os.path import dirname
from collections import OrderedDict
from pycmm.template import SafeTester
from pycmm.utils.xlsutils import XlsUtils
from pycmm.settings import SLURM_MUTREP_TEST
from pycmm.settings import FAST_PROJECT_CODE
from pycmm.settings import SLOW_PROJECT_CODE
from pycmm.settings import ALL_MUTREP_ANNO_COLS
from pycmm.settings import DFLT_MUTREP_ALLOC_TIME
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
from pycmm.flow.mutrep import MutRepPipeline
from pycmm.flow.mutrep import create_jobs_setup_file
from pycmm.flow.mutrep import JOBS_SETUP_RPT_FILTER_RARE
from pycmm.flow.mutrep import JOBS_SETUP_RPT_FILTER_NON_INTERGENIC
from pycmm.flow.mutrep import JOBS_SETUP_RPT_FILTER_NON_INTRONIC
from pycmm.flow.mutrep import JOBS_SETUP_RPT_FILTER_HAS_MUTATION
from pycmm.flow.mutrep import JOBS_SETUP_RPT_FILTER_HAS_SHARED

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

#class TestMutRepPipelineSlurm(SafeTester):
#
#    def __init__(self, methodName):
#        super(TestMutRepPipelineSlurm, self).__init__(methodName=methodName,
#                                                      test_module_name=__name__,
#                                                      )
#
#    def setUp(self):
#        pass
#
#    def __create_jobs_setup_file(self,
#                                 dataset_name=None,
#                                 project_code=SLOW_PROJECT_CODE,
#                                 rpt_alloc_time=DFLT_MUTREP_ALLOC_TIME,
#                                 sample_infos=None,
#                                 anno_cols=DFLT_TEST_MUTREP_COLS,
#                                 anno_excl_tags=None,
#                                 annotated_vcf_tabix=None,
#                                 report_regions=DFLT_TEST_REPORT_REGIONS,
#                                 frequency_ratios=None,
#                                 split_chrom=False,
#                                 summary_families_sheet=False,
#                                 call_detail=False,
#                                 rows_filter_actions=None,
#                                 only_summary=False,
#                                 only_families=False,
#                                 ):
#        jobs_setup_file = join_path(self.working_dir,
#                                    self.test_function+'_jobs_setup.txt')
#        if dataset_name is None:
#            dataset_name = self.test_function
#        create_jobs_setup_file(dataset_name=dataset_name,
#                               project_out_dir=self.working_dir,
#                               sample_infos=sample_infos,
#                               project_code=project_code,
#                               rpt_alloc_time=rpt_alloc_time,
#                               anno_cols=",".join(anno_cols),
#                               anno_excl_tags=anno_excl_tags,
#                               annotated_vcf_tabix=annotated_vcf_tabix,
#                               report_regions=report_regions,
#                               frequency_ratios=frequency_ratios,
#                               split_chrom=split_chrom,
#                               summary_families_sheet=summary_families_sheet,
#                               call_detail=call_detail,
#                               rows_filter_actions=rows_filter_actions,
#                               only_summary=only_summary,
#                               only_families=only_families,
#                               out_jobs_setup_file=jobs_setup_file,
#                               )
#        return jobs_setup_file
#
#    def test_load_jobs_info_1(self):
#        """ test if layout configurations are loaded correctly """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file()
#        pl = MutRepPipeline(jobs_setup_file)
#        self.assertEqual(pl.rpt_alloc_time,
#                         DFLT_MUTREP_ALLOC_TIME,
#                         "MutRepPipeline cannot correctly read meta info 'report allocation time' from jobs setup file")
#
#    def test_load_jobs_info_2(self):
#        """ test if non-default layout configurations are loaded correctly """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        rows_filter_actions = JOBS_SETUP_RPT_FILTER_RARE
#        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_INTERGENIC
#        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_NON_INTRONIC
#        rows_filter_actions += ',' + JOBS_SETUP_RPT_FILTER_HAS_MUTATION
#        jobs_setup_file = self.__create_jobs_setup_file(rpt_alloc_time="2-00:00:00",
#                                                        )
#        pl = MutRepPipeline(jobs_setup_file)
#        self.assertEqual(pl.rpt_alloc_time,
#                         "2-00:00:00",
#                         "MutRepPipeline cannot correctly read meta info 'report allocation time' from jobs setup file")
#
#    @unittest.skipUnless(SLURM_MUTREP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_summary_reports_1(self):
#        """ test generating slurm summary reports (3 families) w/o splite_chrom """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        annotated_vcf_tabix = join_path(self.data_dir,
#                                        "chr6_18.vcf.gz")
#        rpt_out_file = join_path(self.working_dir,
#                                 self.current_func_name + ".xlsx")
#        jobs_setup_file = self.__create_jobs_setup_file(rpt_alloc_time="10:00:00",
#                                                        annotated_vcf_tabix=annotated_vcf_tabix,
#                                                        report_regions="6:78171941-78172992,18:28610988-28611790",
#                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
#                                                        call_detail="YES",
#                                                        )
#        self.dbg(annotated_vcf_tabix)
#        pl = MutRepPipeline(jobs_setup_file)
#        pl.gen_summary_reports()
#
#    @unittest.skipUnless(SLURM_MUTREP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_summary_reports_2(self):
#        """ test generating slurm summary reports (3 families) with splite_chrom w/o report regions """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        annotated_vcf_tabix = join_path(self.data_dir,
#                                        "chr6_18.vcf.gz")
#        rpt_out_file = join_path(self.working_dir,
#                                 self.current_func_name + ".xlsx")
#        jobs_setup_file = self.__create_jobs_setup_file(annotated_vcf_tabix=annotated_vcf_tabix,
#                                                        report_regions=None,
#                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
#                                                        call_detail="YES",
#                                                        split_chrom=True,
#                                                        )
#        pl = MutRepPipeline(jobs_setup_file)
#        pl.gen_summary_reports()
#
#    @unittest.skipUnless(SLURM_MUTREP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_summary_reports_3(self):
#        """ test generating slurm summary reports (3 families) with splite_chrom with report regions """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        annotated_vcf_tabix = join_path(self.data_dir,
#                                        "chr6_18.vcf.gz")
#        rpt_out_file = join_path(self.working_dir,
#                                 self.current_func_name + ".xlsx")
#        jobs_setup_file = self.__create_jobs_setup_file(annotated_vcf_tabix=annotated_vcf_tabix,
#                                                        report_regions="6:78171941-78172992,18:28610988-28611790",
#                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
#                                                        call_detail="YES",
#                                                        split_chrom=True,
#                                                        )
#        pl = MutRepPipeline(jobs_setup_file)
#        pl.gen_summary_reports()
#
#    @unittest.skipUnless(SLURM_MUTREP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_summary_reports_4(self):
#        """ test generating slurm summary reports for sporadic samples """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        annotated_vcf_tabix = join_path(self.data_dir,
#                                        "input.vcf.gz")
#        rpt_out_file = join_path(self.working_dir,
#                                 self.current_func_name + ".xlsx")
#        anno_cols = list(DFLT_TEST_MUTREP_COLS)
#        anno_cols.append("LRT_pred")
#        anno_cols.append("LRT_score")
#        jobs_setup_file = self.__create_jobs_setup_file(annotated_vcf_tabix=annotated_vcf_tabix,
#                                                        anno_cols=anno_cols,
#                                                        project_code=None,
#                                                        report_regions=None,
#                                                        call_detail=False,
#                                                        sample_infos="1003-06o,1068-05o,1094-08F,1199-05o,1205-05o,1307-09D,142-09F,153-09F,2014-04018,230-06o,258-06o,264-08F,2750-13D,294-07LF,295-08F,296-09F,323-09F,329-09F,344-07LF,376-05o,386-04o,4270-11D,4283-13D,4735-11D,4773-11D,48-07LF,49-06o,5053-10D,566-04o,569-04o,59-04o,61-09F,664-08F,679-05o,711-14D,728-05o,749-08F,760-04o,771-08F,784-05o,813-06o,91-04o,95-04o,Co-1161,Co-1190,Co-121,Co-1524,Co-245,Co-278,Co-286,Co-489",
#                                                        )
#        pl = MutRepPipeline(jobs_setup_file)
#        pl.gen_summary_reports()
#
#    @unittest.skipUnless(SLURM_MUTREP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_families_reports_1(self):
#        """ test generating slurm families reports (3 families) w/o splite_chrom """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        annotated_vcf_tabix = join_path(self.data_dir,
#                                        "input.vcf.gz")
#        jobs_setup_file = self.__create_jobs_setup_file(annotated_vcf_tabix=annotated_vcf_tabix,
#                                                        report_regions="6",
#                                                        call_detail="YES",
#                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
#                                                        )
#        pl = MutRepPipeline(jobs_setup_file)
#        pl.gen_families_reports()
#
#    @unittest.skipUnless(SLURM_MUTREP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_families_reports_2(self):
#        """ test generating slurm families reports (3 families) with splite_chrom w/o report_regions """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        annotated_vcf_tabix = join_path(self.data_dir,
#                                        "chr6_18.vcf.gz")
#        jobs_setup_file = self.__create_jobs_setup_file(annotated_vcf_tabix=annotated_vcf_tabix,
#                                                        report_regions=None,
#                                                        call_detail="YES",
#                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
#                                                        split_chrom=True,
#                                                        )
#        pl = MutRepPipeline(jobs_setup_file)
#        pl.gen_families_reports()
#
#    @unittest.skipUnless(SLURM_MUTREP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_families_reports_3(self):
#        """ test generating slurm families reports (3 families) with splite_chrom and with report_regions """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        annotated_vcf_tabix = join_path(self.data_dir,
#                                        "chr6_18.vcf.gz")
#        jobs_setup_file = self.__create_jobs_setup_file(annotated_vcf_tabix=annotated_vcf_tabix,
#                                                        report_regions="6:78171941-78172992,18:28610988-28611790",
#                                                        call_detail="YES",
#                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
#                                                        split_chrom=True,
#                                                        )
#        pl = MutRepPipeline(jobs_setup_file)
#        pl.gen_families_reports()
#
#    def tearDown(self):
#        self.remove_working_dir()
