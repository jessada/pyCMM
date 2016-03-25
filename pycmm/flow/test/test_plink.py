# This Python file uses the following encoding: utf-8
import unittest
from os.path import join as join_path
from os.path import dirname
from pycmm.settings import FULL_SYSTEM_TEST
from pycmm.template import SafeTester
from pycmm.utils.dnalib import ALL_CHROMS
from pycmm.utils.plinkutils import HapAssocUtils
from pycmm.flow.plink import PlinkPipeline
from pycmm.flow.plink import create_jobs_setup_file
from pycmm.flow.plink import DFLT_CUTOFF_PVALUE
from pycmm.flow.plink import DFLT_HAP_WINDOW_SIZES

PLINK_TEST = True

class TestPlinkPipeline(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            test_module_name=__name__,
                            )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self,
                                 project_name=None,
                                 project_out_dir=None,
                                 input_file_prefix=None,
                                 input_binary=None,
                                 input_dna_regions=None,
                                 cutoff_pvalue=None,
                                 hap_window_sizes=None,
                                 project_code=None,
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
                               input_dna_regions=input_dna_regions,
                               cutoff_pvalue=cutoff_pvalue,
                               hap_window_sizes=hap_window_sizes,
                               project_code=project_code,
                               )
        return jobs_setup_file

    def test_load_jobs_info_1(self):
        """ test if default PLINK pipeline configurations are loaded correctly """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        pl = PlinkPipeline(jobs_setup_file)
        self.assertEqual(pl.project_out_dir,
                         self.working_dir,
                         "PlinkPipeline cannot correctly identify 'output dir' from jobs setup file")
        self.assertEqual(pl.project_name,
                         self.current_func_name,
                         "PlinkPipeline cannot correctly identify 'dataset name' from jobs setup file")
        self.assertEqual(pl.project_code,
                         None,
                         "PlinkPipeline cannot correctly identify 'project code' from jobs setup file")
        exp_input_file_prefix = join_path(self.data_dir,
                                          "input")
        self.assertEqual(pl.plink_params.input_file_prefix,
                         exp_input_file_prefix,
                         "PlinkPipeline cannot correctly identify 'input file prefix' from jobs setup file")
        self.assertEqual(pl.plink_params.input_binary,
                         True,
                         "PlinkPipeline cannot correctly identify 'input binary' from jobs setup file")
        self.assertEqual(pl.plink_params.input_dna_regions[0].chrom,
                         "1",
                         "PlinkPipeline cannot correctly identify 'dna regions' from jobs setup file")
        self.assertEqual(pl.plink_params.input_dna_regions[5].chrom,
                         "6",
                         "PlinkPipeline cannot correctly identify 'dna regions' from jobs setup file")
        self.assertEqual(pl.plink_params.input_dna_regions[22].chrom,
                         "MT",
                         "PlinkPipeline cannot correctly identify 'dna regions' from jobs setup file")
        self.assertEqual(pl.plink_params.input_dna_regions[23].chrom,
                         "X",
                         "PlinkPipeline cannot correctly identify 'dna regions' from jobs setup file")
        self.assertEqual(pl.plink_params.input_dna_regions[24].chrom,
                         "Y",
                         "PlinkPipeline cannot correctly identify 'dna regions' from jobs setup file")
        self.assertEqual(pl.plink_params.cutoff_pvalue,
                         DFLT_CUTOFF_PVALUE,
                         "PlinkPipeline cannot correctly identify 'cutoff p-value' from jobs setup file")
        self.assertEqual(pl.plink_params.hap_window_sizes,
                         DFLT_HAP_WINDOW_SIZES,
                         "PlinkPipeline cannot correctly identify 'haplotype window size' from jobs setup file")
        
    def test_load_jobs_info_2(self):
        """ test if non-default PLINK pipeline configurations are loaded correctly """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(project_code="b2011097",
                                                        input_file_prefix="/tmp/file",
                                                        input_binary=False,
                                                        input_dna_regions="14:3456-7890",
                                                        cutoff_pvalue="0.001",
                                                        hap_window_sizes="3",
                                                        )
        pl = PlinkPipeline(jobs_setup_file)
        self.assertEqual(pl.project_code,
                         "b2011097",
                         "PlinkPipeline cannot correctly identify 'project code' from jobs setup file")
        exp_input_file_prefix = "/tmp/file"
        self.assertEqual(pl.plink_params.input_file_prefix,
                         exp_input_file_prefix,
                         "PlinkPipeline cannot correctly identify 'input file prefix' from jobs setup file")
        self.assertEqual(pl.plink_params.input_binary,
                         False,
                         "PlinkPipeline cannot correctly identify 'input binary' from jobs setup file")
        self.assertEqual(pl.plink_params.input_dna_regions[0].chrom,
                         "14",
                         "PlinkPipeline cannot correctly identify 'dna regions chromosome' from jobs setup file")
        self.assertEqual(pl.plink_params.input_dna_regions[0].start_pos,
                         "3456",
                         "PlinkPipeline cannot correctly identify 'dna regions start position' from jobs setup file")
        self.assertEqual(pl.plink_params.input_dna_regions[0].end_pos,
                         "7890",
                         "PlinkPipeline cannot correctly identify 'dna regions end position' from jobs setup file")
        self.assertEqual(pl.plink_params.cutoff_pvalue,
                         "0.001",
                         "PlinkPipeline cannot correctly identify 'cutoff p-value' from jobs setup file")
        self.assertEqual(pl.plink_params.hap_window_sizes,
                         [3],
                         "PlinkPipeline cannot correctly identify 'haplotype window size' from jobs setup file")
        
    def test_load_jobs_info_3(self):
        """ test if special PLINK pipeline configurations are loaded correctly """

        self.init_test(self.current_func_name)
        input_dna_regions = "5"
        input_dna_regions += ",X:1456-2334"
        input_dna_regions += ",Y:3456-7890"
        jobs_setup_file = self.__create_jobs_setup_file(hap_window_sizes="3,5,10,25",
                                                        input_dna_regions=input_dna_regions,
                                                        )
        pl = PlinkPipeline(jobs_setup_file)
        self.assertEqual(pl.plink_params.input_dna_regions[0].chrom,
                         "5",
                         "PlinkPipeline cannot correctly identify 'dna regions chromosome' from jobs setup file")
        self.assertEqual(pl.plink_params.input_dna_regions[0].start_pos,
                         None,
                         "PlinkPipeline cannot correctly identify 'dna regions start position' from jobs setup file")
        self.assertEqual(pl.plink_params.input_dna_regions[0].end_pos,
                         None,
                         "PlinkPipeline cannot correctly identify 'dna regions end position' from jobs setup file")
        self.assertEqual(pl.plink_params.input_dna_regions[1].chrom,
                         "X",
                         "PlinkPipeline cannot correctly identify 'dna regions chromosome' from jobs setup file")
        self.assertEqual(pl.plink_params.input_dna_regions[1].start_pos,
                         "1456",
                         "PlinkPipeline cannot correctly identify 'dna regions start position' from jobs setup file")
        self.assertEqual(pl.plink_params.input_dna_regions[1].end_pos,
                         "2334",
                         "PlinkPipeline cannot correctly identify 'dna regions end position' from jobs setup file")
        self.assertEqual(pl.plink_params.input_dna_regions[2].chrom,
                         "Y",
                         "PlinkPipeline cannot correctly identify 'dna regions chromosome' from jobs setup file")
        self.assertEqual(pl.plink_params.input_dna_regions[2].start_pos,
                         "3456",
                         "PlinkPipeline cannot correctly identify 'dna regions start position' from jobs setup file")
        self.assertEqual(pl.plink_params.input_dna_regions[2].end_pos,
                         "7890",
                         "PlinkPipeline cannot correctly identify 'dna regions end position' from jobs setup file")
        self.assertEqual(pl.plink_params.hap_window_sizes,
                         [3,5,10,25],
                         "PlinkPipeline cannot correctly identify 'haplotype window size' from jobs setup file")
        self.assertEqual(len(pl.plink_params.hap_window_sizes),
                         4,
                         "PlinkPipeline cannot correctly identify 'haplotype window size' from jobs setup file")

    @unittest.skipUnless(FULL_SYSTEM_TEST or PLINK_TEST, "taking too long time to test")
    def test_hap_assocs_offline_1(self):
        """
        test plink haplotype association with very basic paramiters
        - one region with both chromosome and positions
        - one haplotype window (1)
        - w/o phenotype file
        - w/o sample info
        - with all other default features
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(input_dna_regions="9:100911000-100946000",
                                                        )
        pl = PlinkPipeline(jobs_setup_file)
        hap_assoc_out = pl.run_hap_assocs_offline()
        self.dbg(hap_assoc_out)
        hat = HapAssocUtils(hap_assoc_out[0][0])
        self.assertEqual(hat.nlines,
                         20,
                         "PlinkPipeline cannot correctly perform 'haplotype association study'")
        self.assertEqual(hat.ncols,
                         8,
                         "PlinkPipeline cannot correctly perform 'haplotype association study'")

    @unittest.skipUnless(FULL_SYSTEM_TEST or PLINK_TEST, "taking too long time to test")
    def test_hap_assocs_offline_2(self):
        """
        test plink haplotype association with a little advanced paramiters
        - one region with both chromosome and positions
        - two haplotype window (1,3)
        - w/o phenotype file
        - w/o sample info
        - with all other default features
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(input_dna_regions="9:100911000-100946000",
                                                        hap_window_sizes="1,3",
                                                        )
        pl = PlinkPipeline(jobs_setup_file)
        hap_assoc_out = pl.run_hap_assocs_offline()
        hat = HapAssocUtils(hap_assoc_out[0][0])
        self.assertEqual(hat.nlines,
                         20,
                         "PlinkPipeline cannot correctly perform 'haplotype association study'")
        self.assertEqual(hat.ncols,
                         8,
                         "PlinkPipeline cannot correctly perform 'haplotype association study'")
        hat = HapAssocUtils(hap_assoc_out[0][1])
        self.assertEqual(hat.nlines,
                         46,
                         "PlinkPipeline cannot correctly perform 'haplotype association study'")
        self.assertEqual(hat.ncols,
                         8,
                         "PlinkPipeline cannot correctly perform 'haplotype association study'")

    def tearDown(self):
        self.remove_working_dir()
