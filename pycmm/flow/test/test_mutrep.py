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
DFLT_TEST_MUTREP_COLS.append("ExonicFunc.refGene")
DFLT_TEST_MUTREP_COLS.append("Gene.refGene")
DFLT_TEST_MUTREP_COLS.append("GeneDetail.refGene")
DFLT_TEST_MUTREP_COLS.append("cytoBand")
DFLT_TEST_MUTREP_COLS.append("1000g2014oct_all")

DFLT_TEST_REPORT_REGIONS = "18:12512255-14542551"
DFLT_TEST_FREQ_RATIOS = DFLT_MUTREP_FREQ_RATIOS

class TestMutRepPipeline(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            test_module_name=__name__,
                            )

    def setUp(self):
        mylogger.getLogger(__name__)

    def __create_jobs_setup_file(self,
                                 vcf_tabix_file,
                                 dataset_name=None,
                                 project_code=SLOW_PROJECT_CODE,
                                 vcf_region="18",
                                 sample_infos=None,
                                 anno_cols=DFLT_TEST_MUTREP_COLS,
                                 annotated_vcf_tabix=None,
                                 report_regions=DFLT_TEST_REPORT_REGIONS,
                                 call_info="NO",
                                 frequency_ratios=DFLT_TEST_FREQ_RATIOS,
                                 split_chrom=False,
                                 exclude_common=False,
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        if dataset_name is None:
            dataset_name = self.test_function
        create_jobs_setup_file(dataset_name=dataset_name,
                               project_out_dir=self.working_dir,
                               vcf_tabix_file=vcf_tabix_file,
                               vcf_region=vcf_region,
                               sample_infos=sample_infos,
                               project_code=project_code,
                               anno_cols=",".join(anno_cols),
                               annotated_vcf_tabix=annotated_vcf_tabix,
                               report_regions=report_regions,
                               call_info=call_info,
                               frequency_ratios=frequency_ratios,
                               split_chrom=split_chrom,
                               exclude_common=exclude_common,
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
        self.assertEqual(pl.report_layout.anno_cols[3],
                         "GeneDetail.refGene",
                         "MutRepPipeline cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(pl.report_layout.anno_cols[3],
                         "GeneDetail.refGene",
                         "MutRepPipeline cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(pl.report_layout.anno_cols[4],
                         "cytoBand",
                         "MutRepPipeline cannot correctly read report layout info 'layout columns' from jobs setup file")
        self.assertEqual(pl.report_layout.report_regions[0].chrom,
                         "18",
                         "MutRepPipeline cannot correctly read report layout info 'report regions' from jobs setup file")
        self.assertEqual(pl.annotated_vcf_tabix,
                         pl.annovar_config.annotated_vcf + ".gz",
                         "MutRepPipeline cannot correctly determine report layout info 'annotated vcf tabix' file")
        self.assertEqual(pl.report_layout.call_info,
                         False,
                         "MutRepPipeline cannot correctly read report layout info 'call info' from jobs setup file")
        self.assertEqual(pl.report_layout.split_chrom,
                         False,
                         "MutRepPipeline cannot correctly read report layout info 'split chrom' from jobs setup file")
        self.assertEqual(pl.report_layout.exclude_common,
                         False,
                         "MutRepPipeline cannot correctly read report layout info 'exclude common' from jobs setup file")

    def test_load_jobs_info_2(self):
        """ test if non-default layout configurations are loaded correctly """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        dummy_vcf_tabix_file = "/path/to/vcf_tabix_file"
        dummy_annotated_vcf_tabix = "/path/to/annotated_vcf_tabix.vcf.gz"
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=dummy_vcf_tabix_file,
                                                        annotated_vcf_tabix=dummy_annotated_vcf_tabix,
                                                        report_regions="6:78161823-78164117,"+DFLT_TEST_REPORT_REGIONS+",22",
                                                        call_info="YES",
                                                        frequency_ratios=None,
                                                        split_chrom=True,
                                                        exclude_common=True,
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
        self.assertEqual(pl.annotated_vcf_tabix,
                         dummy_annotated_vcf_tabix,
                         "MutRepPipeline cannot correctly determine report layout info 'annotated vcf tabix' file")
        self.assertEqual(pl.report_layout.call_info,
                         True,
                         "MutRepPipeline cannot correctly read report layout info 'call info' from jobs setup file")
        self.assertEqual(pl.report_layout.split_chrom,
                         True,
                         "MutRepPipeline cannot correctly read report layout info 'split chrom' from jobs setup file")
        self.assertEqual(pl.report_layout.exclude_common,
                         True,
                         "MutRepPipeline cannot correctly read report layout info 'exclude common' from jobs setup file")

    def test_load_jobs_info_3(self):
        """ test if non-default (None) layout configurations are loaded correctly """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        dummy_vcf_tabix_file = "/path/to/vcf_tabix_file"
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=dummy_vcf_tabix_file,
                                                        report_regions=None,
                                                        frequency_ratios="ExAC:0.5",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        self.assertEqual(pl.report_layout.report_regions,
                         None,
                         "MutRepPipeline cannot correctly read report layout info 'report regions' from jobs setup file")
        self.assertEqual(pl.report_layout.freq_ratios["ExAC"],
                         0.5,
                         "MutRepPipeline cannot correctly read report layout info 'frequency ratios' from jobs setup file")

#    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_summary_report_1(self):
        """ test if summary report with default configuration can be correctly generated """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "chr6_18.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "chr6_18.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        project_code=None,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions=None,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)

    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_summary_report_2(self):
        """ test if summary reprot with custom configuration can be correctly generated """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "chr18.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "chr18.vcf.gz")
        anno_cols=DFLT_TEST_MUTREP_COLS
        anno_cols.remove("cytoBand")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        project_code=None,
                                                        anno_cols=anno_cols,
                                                        call_info="YES",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)

#    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_summary_report_3(self):
        """ test summary with multiple report_regions """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "chr6_18.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "chr6_18.vcf.gz")
        rpt_out_file = join_path(self.working_dir,
                                 self.current_func_name + ".xlsx")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        project_code=None,
                                                        report_regions="6:78171941-78172992,18:28610988-28611790",
                                                        call_info="YES",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions, out_file=rpt_out_file)

#    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_summary_report_4(self):
        """
        test if mutations with more than one alternate alleles in summary
        is correctly generated
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        project_code=None,
                                                        report_regions="6",
                                                        call_info="NO",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)

#    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_summary_report_5(self):
        """ test summary with multiple report_regions and many sample infos """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "chr6_18.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "chr6_18.vcf.gz")
        rpt_out_file = join_path(self.working_dir,
                                 self.current_func_name + ".xlsx")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        project_code=None,
                                                        report_regions="6:78171941-78172992,18:28610988-28611790",
                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
                                                        call_info="YES",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)

#    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_summary_report_6(self):
        """ test summary report that will exclude common mutation MAF > 0.2  """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        project_code=None,
                                                        report_regions="6",
                                                        call_info="NO",
                                                        exclude_common=True,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_report(pl.report_layout.report_regions)

#    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_summary_reports_1(self):
        """ test generating offline summary report (3 families)"""

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "chr6_18.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "chr6_18.vcf.gz")
        rpt_out_file = join_path(self.working_dir,
                                 self.current_func_name + ".xlsx")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        project_code=None,
                                                        report_regions="6:78171941-78172992,18:28610988-28611790",
                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
                                                        call_info="YES",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_reports()

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_summary_reports_2(self):
        """ test generating slurm summary reports (3 families) w/o splite_chrom """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "chr6_18.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "chr6_18.vcf.gz")
        rpt_out_file = join_path(self.working_dir,
                                 self.current_func_name + ".xlsx")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions="6:78171941-78172992,18:28610988-28611790",
                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
                                                        call_info="YES",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_reports()

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_summary_reports_3(self):
        """ test generating slurm summary reports (3 families) with splite_chrom w/o report regions """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "chr6_18.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "chr6_18.vcf.gz")
        rpt_out_file = join_path(self.working_dir,
                                 self.current_func_name + ".xlsx")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions=None,
                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
                                                        call_info="YES",
                                                        split_chrom=True,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_reports()

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_summary_reports_4(self):
        """ test generating slurm summary reports (3 families) with splite_chrom with report regions """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "chr6_18.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "chr6_18.vcf.gz")
        rpt_out_file = join_path(self.working_dir,
                                 self.current_func_name + ".xlsx")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions="6:78171941-78172992,18:28610988-28611790",
                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
                                                        call_info="YES",
                                                        split_chrom=True,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_summary_reports()

#    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_family_report_1(self):
        """ test with only one family which has only one members """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        project_code=None,
                                                        call_info="YES",
                                                        report_regions="6",
                                                        sample_infos="6789:Al-65",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_family_report('6789', pl.report_layout.report_regions)

#    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_family_report_2(self):
        """ test with only one family which has two members """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        project_code=None,
                                                        call_info="YES",
                                                        report_regions="6",
                                                        sample_infos="1234:Alb-31:Br-466",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_family_report('1234', pl.report_layout.report_regions)

#    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_family_report_3(self):
        """ test with only one family which has three members """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        project_code=None,
                                                        call_info="YES",
                                                        report_regions="6",
                                                        sample_infos="6067:Br-432:Al-161:Br-504",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_family_report('6067', pl.report_layout.report_regions)

#    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_families_reports_1(self):
        """ test generating offline families reports (1 family)"""

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        project_code=None,
                                                        call_info="YES",
                                                        report_regions="6",
                                                        sample_infos="6067:Br-432:Al-161:Br-504",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_families_reports()

#    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_families_reports_2(self):
        """ test generating offline families reports (3 families)"""

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        project_code=None,
                                                        report_regions="6",
                                                        call_info="YES",
                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_families_reports()

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_families_reports_3(self):
        """ test generating slurm families reports (3 families) w/o splite_chrom """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions="6",
                                                        call_info="YES",
                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_families_reports()

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_families_reports_4(self):
        """ test generating slurm families reports (3 families) with splite_chrom w/o report_regions """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "chr6_18.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "chr6_18.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions=None,
                                                        call_info="YES",
                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
                                                        split_chrom=True,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_families_reports()

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_families_reports_5(self):
        """ test generating slurm families reports (3 families) with splite_chrom and with report_regions """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "chr6_18.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "chr6_18.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        report_regions="6:78171941-78172992,18:28610988-28611790",
                                                        call_info="YES",
                                                        sample_infos="1234:Alb-31:Br-466,6067:Br-432:Al-161:Br-504,6789:Al-65",
                                                        split_chrom=True,
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_families_reports()

    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_families_reports_6(self):
        """ test generating 101 CRC families reports """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        vcf_tabix_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        annotated_vcf_tabix = join_path(self.data_dir,
                                        "input.vcf.gz")
        jobs_setup_file = self.__create_jobs_setup_file(vcf_tabix_file=vcf_tabix_file,
                                                        annotated_vcf_tabix=annotated_vcf_tabix,
                                                        project_code=None,
                                                        report_regions=None,
                                                        call_info="NO",
                                                        sample_infos="8:Co-35:Co-37,12:Co-89:Co-90,13:Co-95,26:Co-131:Co-135,31:1793-11D:1322-11D,87:Co-218:Co-258,91:Co-454:Co-700,94:Co-238,110:1526-02D:Co-1301,134:Co-460:Co-553,141:Co-305:Co-785,185:Co-603:Co-669,191:Co-384,214:Co-484,216:Co-367:Co-446,221:Co-358,227:Co-364,231:Co-555:Co-572,254:Co-616:Co-1156,275:Co-618:Co-1262,288:Co-1141,296:Co-793:Co-876,301:Co-837:Co-840:Co-1053,306:Co-779,309:Co-783,312:Co-1116,315:1462-01D,325:Co-851:Co-859,348:Co-846:Co-857,350:1104-03D:Co-866,409:Co-1254,415:Co-1031:Co-1037,425:Co-1458:Co-1595,434:Co-1051:Co-1534,445:Co-1157:Co-1158,478:Co-1207:Co-1274,485:Co-1302:Co-1322,532:Co-1583:Co-1584,574:468-04:474-05,578:531-04o:Co-1349,650:398-05o:729-05o,695:Co-1354:Co-1359:Co-1368,739:529-05:Co-1467,740:602-05o:Co-1373:Co-1383,849:Co-1764:Co-1765,869:Co-1685,871:Co-1618:Co-1661,918:134-06:354-06,975:Co-1591:Co-1600,1025:Co-1529,1085:Co-1518,1113:642-06:Co-1538,1206:1052-05D:Co-1552,1207:2818-07D,1213:Co-1666,1252:Co-1719,1290:Co-1723,prostate:P001:P002:P003",
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        pl.gen_families_reports()
        mylogger.debug(len(pl.samples_list))
#        pl.samples_list.sort
#        mylogger.debug(",".join(pl.samples_list))

    def tearDown(self):
        self.remove_working_dir()
