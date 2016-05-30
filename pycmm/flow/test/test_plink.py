# This Python file uses the following encoding: utf-8
import unittest
from os.path import join as join_path
from os.path import dirname
from pycmm.settings import FULL_SYSTEM_TEST
from pycmm.template import SafeTester
from pycmm.cmmlib.dnalib import ALL_CHROMS
from pycmm.cmmlib.plinklib import HapAssocUtils
from pycmm.cmmlib.plinklib import HapAssocReader
from pycmm.cmmlib.plinklib import LMissReader
from pycmm.cmmlib.plinklib import MapReader
from pycmm.flow.plink import PlinkPipeline
from pycmm.flow.plink import create_jobs_setup_file
from pycmm.flow.plink import DFLT_CUTOFF_PVALUE
from pycmm.flow.plink import DFLT_HAP_WINDOW_SIZES
from pycmm.flow.plink import JOBS_SETUP_HAP_ASSOC_FILTER_PVALUE005
from pycmm.flow.plink import JOBS_SETUP_HAP_ASSOC_FILTER_DISEASE_SNP

PLINK_TEST = False

class TestPlinkPipeline(SafeTester):

    def __init__(self, methodName):
        super(TestPlinkPipeline, self).__init__(methodName=methodName,
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
                                 hap_window_sizes=None,
                                 filter_criteria=None,
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
                               hap_window_sizes=hap_window_sizes,
                               filter_criteria=filter_criteria,
                               project_code=project_code,
                               out_jobs_setup_file=jobs_setup_file,
                               )
        return jobs_setup_file

    def test_load_jobs_info_1(self):
        """ test if default PLINK pipeline configurations are loaded correctly """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        pl = PlinkPipeline(jobs_setup_file=jobs_setup_file)
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
        exp_phenotype_file = join_path(self.data_dir,
                                       "phenotype.txt")
        self.assertEqual(pl.plink_params.phenotype_file,
                         exp_phenotype_file,
                         "PlinkPipeline cannot correctly identify 'phenotype file' from jobs setup file")
#        self.assertEqual(pl.plink_params.cutoff_pvalue,
#                         DFLT_CUTOFF_PVALUE,
#                         "PlinkPipeline cannot correctly identify 'cutoff p-value' from jobs setup file")
        self.assertEqual(pl.plink_params.hap_window_sizes,
                         DFLT_HAP_WINDOW_SIZES,
                         "PlinkPipeline cannot correctly identify 'haplotype window size' from jobs setup file")
        self.assertEqual(pl.plink_params.filter_criteria,
                         None,
                         "PlinkPipeline cannot correctly identify 'filter criteria' from jobs setup file")
        
    def test_load_jobs_info_2(self):
        """ test if non-default PLINK pipeline configurations are loaded correctly """

        self.init_test(self.current_func_name)
        filter_criteria = JOBS_SETUP_HAP_ASSOC_FILTER_PVALUE005
        filter_criteria += "," + JOBS_SETUP_HAP_ASSOC_FILTER_DISEASE_SNP
        jobs_setup_file = self.__create_jobs_setup_file(project_code="b2011097",
                                                        input_file_prefix="/tmp/file",
                                                        input_binary=False,
                                                        input_dna_regions="14:3456-7890",
                                                        cutoff_pvalue="0.001",
                                                        hap_window_sizes="3",
                                                        filter_criteria=filter_criteria,
                                                        )
        pl = PlinkPipeline(jobs_setup_file=jobs_setup_file)
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
#        self.assertEqual(pl.plink_params.cutoff_pvalue,
#                         "0.001",
#                         "PlinkPipeline cannot correctly identify 'cutoff p-value' from jobs setup file")
        self.assertEqual(pl.plink_params.hap_window_sizes,
                         [3],
                         "PlinkPipeline cannot correctly identify 'haplotype window size' from jobs setup file")
        self.assertEqual(pl.plink_params.filter_criteria.filter_pvalue005,
                         True,
                         "PlinkPipeline cannot correctly identify 'filter pvalue 0.05' from jobs setup file")
        self.assertEqual(pl.plink_params.filter_criteria.filter_disease_snp,
                         True,
                         "PlinkPipeline cannot correctly identify 'filter disease snp' from jobs setup file")
        
    def test_load_jobs_info_3(self):
        """ test if special PLINK pipeline configurations are loaded correctly """

        self.init_test(self.current_func_name)
        input_dna_regions = "5"
        input_dna_regions += ",X:1456-2334"
        input_dna_regions += ",Y:3456-7890"
        jobs_setup_file = self.__create_jobs_setup_file(hap_window_sizes="3,5,10,25",
                                                        input_dna_regions=input_dna_regions,
                                                        )
        pl = PlinkPipeline(jobs_setup_file=jobs_setup_file)
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

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(input_dna_regions="9:100911000-100946000",
                                                        )
        pl = PlinkPipeline(jobs_setup_file=jobs_setup_file)
        hap_assoc_out = pl.run_hap_assocs_offline()
        hat = HapAssocUtils(hap_assoc_out[pl.region_keys[0]][0])
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

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(input_dna_regions="9:100911000-100946000",
                                                        hap_window_sizes="1,3",
                                                        )
        pl = PlinkPipeline(jobs_setup_file=jobs_setup_file)
        hap_assoc_out = pl.run_hap_assocs_offline()
        hat = HapAssocUtils(hap_assoc_out[pl.region_keys[0]][0])
        self.assertEqual(hat.nlines,
                         20,
                         "PlinkPipeline cannot correctly perform 'haplotype association study'")
        self.assertEqual(hat.ncols,
                         8,
                         "PlinkPipeline cannot correctly perform 'haplotype association study'")
        hat = HapAssocUtils(hap_assoc_out[pl.region_keys[0]][1])
        self.assertEqual(hat.nlines,
                         46,
                         "PlinkPipeline cannot correctly perform 'haplotype association study'")
        self.assertEqual(hat.ncols,
                         8,
                         "PlinkPipeline cannot correctly perform 'haplotype association study'")

    @unittest.skipUnless(FULL_SYSTEM_TEST or PLINK_TEST, "taking too long time to test")
    def test_merge_hap_assocs_1(self):
        """ test basic merging assoc.hap files """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        pl = PlinkPipeline(jobs_setup_file=jobs_setup_file)
        input_files = []
        input_files.append(join_path(self.data_dir,
                                     "input1.assoc.hap"))
        input_files.append(join_path(self.data_dir,
                                     "input2.assoc.hap"))
        locus_prefixs = ["W1", "W3"]
        out_file = join_path(self.working_dir,
                             "merged_hap_assoc.assoc.hap"
                             )
        pl.merge_hap_assocs(input_files, out_file, locus_prefixs=locus_prefixs)
        hap_assoc_reader = HapAssocReader(file_name=out_file)
        self.assertEqual(len(list(hap_assoc_reader)),
                         58,
                         "merge_hap_assoc doesn't work correctly")
        hap_assoc_reader = HapAssocReader(file_name=out_file)
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.chisq,
                         "CHISQ",
                         "merge_hap_assoc doesn't work correctly")
        self.assertEqual(hap_assoc_rec.ors,
                         "ORS",
                         "merge_hap_assoc doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.locus,
                         "W1_WIN4",
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.ors,
                         1.1606,
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.snps,
                         "rs2417733",
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.locus,
                         "W3_WIN1",
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.snps,
                         "rs2795491|rs2795492|rs2250494",
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.ors,
                         0.9937,
                         "HapAssocReader doesn't work correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST or PLINK_TEST, "taking too long time to test")
    def test_merge_hap_assocs_2(self):
        """ test merging assoc.hap files with filtering criteria """

        self.init_test(self.current_func_name)
        filter_criteria = JOBS_SETUP_HAP_ASSOC_FILTER_PVALUE005
        filter_criteria += "," + JOBS_SETUP_HAP_ASSOC_FILTER_DISEASE_SNP
        jobs_setup_file = self.__create_jobs_setup_file(filter_criteria=filter_criteria)
        pl = PlinkPipeline(jobs_setup_file=jobs_setup_file)
        input_files = []
        input_files.append(join_path(self.data_dir,
                                     "input1.assoc.hap"))
        input_files.append(join_path(self.data_dir,
                                     "input2.assoc.hap"))
        locus_prefixs = ["W1", "W3"]
        out_file = join_path(self.working_dir,
                             "merged_hap_assoc.assoc.hap"
                             )
        pl.merge_hap_assocs(input_files, out_file, locus_prefixs=locus_prefixs)
        hap_assoc_reader = HapAssocReader(file_name=out_file)
        self.assertEqual(len(list(hap_assoc_reader)),
                         5,
                         "merge_hap_assoc doesn't work correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST or PLINK_TEST, "taking too long time to test")
    def test_extract_snps_stat_1(self):
        """ test extracting SNPs missing genotyping information """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(input_dna_regions="9:100911000-100946000",
                                                        )
        pl = PlinkPipeline(jobs_setup_file=jobs_setup_file)
        snps_stat_out = pl.extract_snps_stat(pl.plink_params.input_dna_regions[0])
        lmiss_reader = LMissReader(file_name=snps_stat_out)
        lmiss_rec = lmiss_reader.next()
        self.assertEqual(lmiss_rec.clst,
                         "CLST",
                         "function extract_snps_stat doesn't work correctly")
        lmiss_rec = lmiss_reader.next()
        self.assertEqual(lmiss_rec.snp,
                         "rs2795491",
                         "function extract_snps_stat doesn't work correctly")
        lmiss_rec = lmiss_reader.next()
        lmiss_rec = lmiss_reader.next()
        lmiss_rec = lmiss_reader.next()
        self.assertEqual(lmiss_rec.n_geno,
                         "200",
                         "function extract_snps_stat doesn't work correctly")
        lmiss_rec = lmiss_reader.next()
        self.assertEqual(lmiss_rec.clst,
                         "2",
                         "function extract_snps_stat doesn't work correctly")
        lmiss_rec = lmiss_reader.next()
        lmiss_rec = lmiss_reader.next()
        lmiss_rec = lmiss_reader.next()
        lmiss_rec = lmiss_reader.next()
        self.assertEqual(lmiss_rec.f_miss,
                         0.005,
                         "function extract_snps_stat doesn't work correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST or PLINK_TEST, "taking too long time to test")
    def test_extract_snps_pos_1(self):
        """ test extracting genotyping position """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(input_dna_regions="9:100911000-100946000",
                                                        )
        pl = PlinkPipeline(jobs_setup_file=jobs_setup_file)
        snps_pos_out = pl.extract_snps_pos(pl.plink_params.input_dna_regions[0])
        map_reader = MapReader(file_name=snps_pos_out)
        map_rec = map_reader.next()
        self.assertEqual(map_rec.snp,
                         "rs2795491",
                         "function extract_snps_pos doesn't work correctly")
        map_rec = map_reader.next()
        map_rec = map_reader.next()
        self.assertEqual(map_rec.chrom,
                         "9",
                         "function extract_snps_pos doesn't work correctly")
        map_rec = map_reader.next()
        self.assertEqual(map_rec.dist,
                         102.26,
                         "function extract_snps_pos doesn't work correctly")
        map_rec = map_reader.next()
        map_rec = map_reader.next()
        map_rec = map_reader.next()
        map_rec = map_reader.next()
        self.assertEqual(map_rec.pos,
                         "100930947",
                         "function extract_snps_pos doesn't work correctly")

    def tearDown(self):
        self.remove_working_dir()
