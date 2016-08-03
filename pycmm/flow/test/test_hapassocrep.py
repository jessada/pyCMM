# This Python file uses the following encoding: utf-8
import unittest
from os.path import join as join_path
from pycmm.settings import FULL_SYSTEM_TEST
from pycmm.template import SafeTester
from pycmm.flow.plink import create_jobs_setup_file
from pycmm.flow.plink import DFLT_CUTOFF_PVALUE
from pycmm.flow.plink import DFLT_CUTOFF_ORS
from pycmm.flow.hapassocrep import HapAssocRepPipeline

class TestHapAssocRepPipeline(SafeTester):

    def __init__(self, methodName):
        super(TestHapAssocRepPipeline, self).__init__(methodName=methodName,
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
                                 phenotype_file=None,
                                 cutoff_pvalue=None,
                                 cutoff_ors=None,
                                 fam_hap_prefix=None,
                                 sample_info=None,
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
        if phenotype_file is None:
            phenotype_file = join_path(self.data_dir,
                                       "phenotype.txt")
        create_jobs_setup_file(project_name=project_name,
                               project_out_dir=project_out_dir,
                               input_file_prefix=input_file_prefix,
                               phenotype_file=phenotype_file,
                               input_binary=input_binary,
                               input_dna_regions=input_dna_regions,
                               cutoff_pvalue=cutoff_pvalue,
                               cutoff_ors=cutoff_ors,
                               fam_hap_prefix=fam_hap_prefix,
                               sample_info=sample_info,
                               )
        return jobs_setup_file

    def test_load_jobs_info_1(self):
        """
        test if default haplotype association study report configurations
        are loaded correctly
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        hap_assoc_file = join_path(self.data_dir,
                                   "input.assoc.hap")
        snp_info_file = join_path(self.data_dir,
                                  "input.snp.info")
        pl = HapAssocRepPipeline(hap_assoc_file=hap_assoc_file,
                                 snp_info_file=snp_info_file,
                                 jobs_setup_file=jobs_setup_file,
                                 )
        self.assertEqual(pl.rpt_params.cutoff_pvalue,
                         DFLT_CUTOFF_PVALUE,
                         "HapAssocRepPipeline cannot correctly identify 'cutoff p-value' from jobs setup file")
        self.assertEqual(pl.rpt_params.cutoff_ors,
                         DFLT_CUTOFF_ORS,
                         "HapAssocRepPipeline cannot correctly identify 'cutoff odds ratio' from jobs setup file")
        self.assertEqual(pl.rpt_params.fam_hap_prefix,
                         None,
                         "HapAssocRepPipeline cannot correctly identify 'familyt haplotype prefix' from jobs setup file")

    def test_load_jobs_info_2(self):
        """
        test if non-default haplotype association study report configurations
        are loaded correctly
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        fam_hap_prefix = join_path(self.data_dir,
                                   "input")
        jobs_setup_file = self.__create_jobs_setup_file(cutoff_pvalue=0.0005,
                                                        cutoff_ors=1.25,
                                                        fam_hap_prefix=fam_hap_prefix,
                                                        )
        hap_assoc_file = join_path(self.data_dir,
                                   "input.assoc.hap")
        snp_info_file = join_path(self.data_dir,
                                  "input.snp.info")
        pl = HapAssocRepPipeline(hap_assoc_file=hap_assoc_file,
                                 snp_info_file=snp_info_file,
                                 jobs_setup_file=jobs_setup_file,
                                 )
        self.assertEqual(pl.rpt_params.cutoff_pvalue,
                         0.0005,
                         "HapAssocRepPipeline cannot correctly identify 'cutoff p-value' from jobs setup file")
        self.assertEqual(pl.rpt_params.cutoff_ors,
                         1.25,
                         "HapAssocRepPipeline cannot correctly identify 'cutoff odds ratio' from jobs setup file")
        self.assertEqual(pl.rpt_params.fam_hap_prefix,
                         fam_hap_prefix,
                         "HapAssocRepPipeline cannot correctly identify 'familyt haplotype prefix' from jobs setup file")

    def test_gen_report_1(self):
        """ test generating basic haplotype associaiton study report """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        fam_hap_prefix = join_path(self.data_dir,
                                   "input")
        sample_info = "1:old_fam1_shared_only"
        sample_info += "," + "3:fam_3"
        sample_info += "," + "9:fam_9"
        sample_info += "," + "10:fam_10"
        sample_info += "," + "13:fam_13"
        sample_info += "," + "test:test_no_exist_sample"
        jobs_setup_file = self.__create_jobs_setup_file(cutoff_pvalue=0.005,
                                                        cutoff_ors=1.07,
                                                        fam_hap_prefix=fam_hap_prefix,
                                                        sample_info=sample_info,
                                                        )
        hap_assoc_file = join_path(self.data_dir,
                                   "input.assoc.hap")
        snp_info_file = join_path(self.data_dir,
                                  "input.snp.info")
        pl = HapAssocRepPipeline(hap_assoc_file=hap_assoc_file,
                                 snp_info_file=snp_info_file,
                                 jobs_setup_file=jobs_setup_file,
                                 )
        out_file = join_path(self.working_dir,
                             self.current_func_name + ".xlsx")
        self.assertEqual(pl.n_snps,
                         10,
                         "HapAssocRepPipeline cannot correctly identify 'number of snps'")
        self.assertEqual(pl.families_info["3"].members[0].gts["rs3780457"][0],
                         "G",
                         "HapAssocRepPipeline cannot correctly load family haplotype information")
        self.assertEqual(pl.families_info["3"].members[0].gts["rs1544205"][0],
                         "A",
                         "HapAssocRepPipeline cannot correctly load family haplotype information")
        self.assertEqual(pl.families_info["3"].members[0].gts["rs2417733"][0],
                         "A",
                         "HapAssocRepPipeline cannot correctly load family haplotype information")
        rpt_file = pl.gen_report(out_file)
        self.warning("************** report has not yet tested *****************")

    def tearDown(self):
        self.remove_working_dir()
