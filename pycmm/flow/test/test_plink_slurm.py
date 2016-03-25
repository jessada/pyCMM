# This Python file uses the following encoding: utf-8
import unittest
from os.path import join as join_path
from os.path import dirname
from pycmm.settings import FULL_SYSTEM_TEST
from pycmm.settings import FAST_PROJECT_CODE
from pycmm.settings import SLOW_PROJECT_CODE
from pycmm.template import SafeTester
from pycmm.utils.dnalib import ALL_CHROMS
from pycmm.utils.plinkutils import HapAssocUtils
from pycmm.flow.plink import PlinkPipeline
from pycmm.flow.plink import create_jobs_setup_file
from pycmm.flow.plink import DFLT_CUTOFF_PVALUE
from pycmm.flow.plink import DFLT_HAP_WINDOW_SIZES

PLINK_TEST = False

class TestPlinkPipelineSlurm(SafeTester):

    def __init__(self, methodName):
        super(TestPlinkPipelineSlurm, self).__init__(methodName=methodName,
                                                     test_module_name=__name__,
                                                     )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self,
                                 project_name=None,
                                 project_out_dir=None,
                                 input_file_prefix=None,
                                 input_binary=None,
                                 dna_regions=None,
                                 cutoff_pvalue=None,
                                 hap_window_sizes=None,
                                 project_code=None,
#                                 sample_infos=None,
#                                 report_regions=DFLT_TEST_REPORT_REGIONS,
#                                 split_chrom=False,
#                                 summary_families_sheet=False,
#                                 only_summary=False,
#                                 only_families=False,
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        if project_name is None:
            project_name = self.test_function
        if project_out_dir is None:
            project_out_dir = self.working_dir
        if input_file_prefix is None:
            input_file_prefix = join_path(self.data_dir,
                                          "input")
        if input_binary is None:
            input_binary = True
        create_jobs_setup_file(project_name=project_name,
                               project_out_dir=project_out_dir,
                               input_file_prefix=input_file_prefix,
                               input_binary=input_binary,
                               dna_regions=dna_regions,
                               cutoff_pvalue=cutoff_pvalue,
                               hap_window_sizes=hap_window_sizes,
                               project_code=project_code,
#                               sample_infos=sample_infos,
#                               report_regions=report_regions,
#                               split_chrom=split_chrom,
#                               summary_families_sheet=summary_families_sheet,
#                               only_summary=only_summary,
#                               only_families=only_families,
#                               out_jobs_setup_file=jobs_setup_file,
                               )
        return jobs_setup_file

    @unittest.skipUnless(FULL_SYSTEM_TEST or PLINK_TEST, "taking too long time to test")
    def test_hap_assoc_n_wins_1_slurm(self):
        """
        test basic plink haplotype association with more than
        2 haplotype windows in SLURM
        - one region with both chromosome and positions
        - w/o phenotype file
        - w/o sample info
        - with all other default features
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(dna_regions="9:100911000-100946000",
                                                        project_code=SLOW_PROJECT_CODE,
                                                        hap_window_sizes="1,3,5",
                                                        )
        pl = PlinkPipeline(jobs_setup_file)
        hap_assoc_out = pl.run_hap_assoc_offline(pl.plink_params.dna_regions[0])
#        hat = HapAssocUtils(hap_assoc_out[0])
#        self.assertEqual(hat.nlines,
#                         20,
#                         "PlinkPipeline cannot correctly perform 'haplotype association study'")
#        self.assertEqual(hat.ncols,
#                         8,
#                         "PlinkPipeline cannot correctly perform 'haplotype association study'")

    def tearDown(self):
        self.remove_working_dir()
