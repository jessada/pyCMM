import unittest
from os.path import join as join_path
from pycmm.template import SafeTester
from pycmm.cmmlib.plinklib import HapAssocReader
from pycmm.cmmlib.plinklib import FamReader
from pycmm.cmmlib.plinklib import TPedReader
from pycmm.cmmlib.plinklib import merge_hap_assocs
from pycmm.cmmlib.plinklib import merge_lmiss_map
from pycmm.cmmlib.plinklib import SnpInfoReader
from pycmm.utils import is_number


class TestHapAssocReader(SafeTester):

    def __init__(self, methodName):
        super(TestHapAssocReader, self).__init__(methodName=methodName,
                                                 test_module_name=__name__,
                                                 )

    def setUp(self):
        pass

    def test_reader_1(self):
        """ test basic reading without any other parameters """

        self.init_test(self.current_func_name)
        # test windowsize = 1
        input_file = join_path(self.data_dir,
                               "input1.assoc.hap")
        hap_assoc_reader = HapAssocReader(file_name=input_file)
        self.assertEqual(len(list(hap_assoc_reader)),
                         21,
                         "HapAssocReader doesn't work correctly")
        hap_assoc_reader = HapAssocReader(file_name=input_file)
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.chisq,
                         "CHISQ",
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.chisq,
                         3.613,
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.f_u,
                         0.6645,
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.snps,
                         "rs2250494",
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.locus,
                         "WIN4",
                         "HapAssocReader doesn't work correctly")
        # test windowsize = 3
        input_file = join_path(self.data_dir,
                               "input2.assoc.hap")
        hap_assoc_reader = HapAssocReader(file_name=input_file)
        self.assertEqual(len(list(hap_assoc_reader)),
                         38,
                         "HapAssocReader doesn't work correctly")
        hap_assoc_reader = HapAssocReader(file_name=input_file)
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.chisq,
                         0.5204,
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.snps,
                         "rs2795491|rs2795492|rs2250494",
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.pvalue,
                         0.001868,
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.haplotype,
                         "CAG",
                         "HapAssocReader doesn't work correctly")
        self.assertEqual(hap_assoc_rec.locus,
                         "WIN3",
                         "HapAssocReader doesn't work correctly")

    def test_reader_2(self):
        """ test basic reading with giving locus prefix """

        self.init_test(self.current_func_name)
        # test windowsize = 1
        input_file = join_path(self.data_dir,
                               "input1.assoc.hap")
        hap_assoc_reader = HapAssocReader(file_name=input_file)
        self.assertEqual(len(list(hap_assoc_reader)),
                         21,
                         "HapAssocReader doesn't work correctly")
        hap_assoc_reader = HapAssocReader(file_name=input_file,
                                          locus_prefix="S1")
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.chisq,
                         "CHISQ",
                         "HapAssocReader doesn't work correctly")
        self.assertEqual(hap_assoc_rec.locus,
                         "LOCUS",
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.locus,
                         "S1_WIN1",
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.chisq,
                         1.095,
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.f_u,
                         0.5622,
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.locus,
                         "S1_WIN3",
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.snps,
                         "rs3780457",
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.locus,
                         "S1_WIN5",
                         "HapAssocReader doesn't work correctly")
        # test windowsize = 3
        input_file = join_path(self.data_dir,
                               "input2.assoc.hap")
        hap_assoc_reader = HapAssocReader(file_name=input_file)
        self.assertEqual(len(list(hap_assoc_reader)),
                         38,
                         "HapAssocReader doesn't work correctly")
        hap_assoc_reader = HapAssocReader(file_name=input_file,
                                          locus_prefix="Win3")
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.locus,
                         "LOCUS",
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.chisq,
                         0.7346,
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.locus,
                         "Win3_WIN1",
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.snps,
                         "rs2795491|rs2795492|rs2250494",
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.locus,
                         "Win3_WIN2",
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.pvalue,
                         0.5963,
                         "HapAssocReader doesn't work correctly")
        hap_assoc_rec = hap_assoc_reader.next()
        self.assertEqual(hap_assoc_rec.haplotype,
                         "AAG",
                         "HapAssocReader doesn't work correctly")
        self.assertEqual(hap_assoc_rec.locus,
                         "Win3_WIN3",
                         "HapAssocReader doesn't work correctly")

    def tearDown(self):
        self.remove_working_dir()

class TestFamReader(SafeTester):

    def __init__(self, methodName):
        super(TestFamReader, self).__init__(methodName=methodName,
                                            test_module_name=__name__,
                                            )

    def setUp(self):
        pass

    def test_reader_1(self):
        """ test basic reading without any other parameters """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "input.tfam")
        fam_reader = FamReader(file_name=input_file)
        self.assertEqual(len(list(fam_reader)),
                         16,
                         "FamReader doesn't work correctly")
        fam_reader = FamReader(file_name=input_file)
        fam_rec = fam_reader.next()
        self.assertEqual(fam_rec.fam_id,
                         "1",
                         "FamReader doesn't work correctly")
        fam_rec = fam_reader.next()
        self.assertEqual(fam_rec.indv_id,
                         "old_fam1_shared_only",
                         "FamReader doesn't work correctly")

class TestTPedReader(SafeTester):

    def __init__(self, methodName):
        super(TestTPedReader, self).__init__(methodName=methodName,
                                             test_module_name=__name__,
                                             )

    def setUp(self):
        pass

    def test_reader_1(self):
        """ test basic reading without any other parameters """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "input.tped")
        tped_reader = TPedReader(file_name=input_file)
        self.assertEqual(len(list(tped_reader)),
                         21,
                         "TPedReader doesn't work correctly")
        tped_reader = TPedReader(file_name=input_file)
        tped_rec = tped_reader.next()
        self.assertEqual(tped_rec.chrom,
                         "9",
                         "TPedReader doesn't work correctly")
        tped_rec = tped_reader.next()
        self.assertEqual(tped_rec.snp,
                         "rs814027",
                         "TPedReader doesn't work correctly")
        tped_rec = tped_reader.next()
        self.assertEqual(tped_rec.dist,
                         105.04,
                         "TPedReader doesn't work correctly")
        tped_rec = tped_reader.next()
        self.assertEqual(tped_rec.pos,
                         "100919318",
                         "TPedReader doesn't work correctly")
        tped_rec = tped_reader.next()
        self.assertEqual(len(tped_rec.gts),
                         16,
                         "TPedReader doesn't work correctly")
        tped_rec = tped_reader.next()
        self.assertEqual(tped_rec.gts[2],
                         "G G",
                         "TPedReader doesn't work correctly")
        self.assertEqual(tped_rec.gts[3],
                         "A G",
                         "TPedReader doesn't work correctly")
        self.assertEqual(tped_rec.gts[6],
                         "G A",
                         "TPedReader doesn't work correctly")
        self.assertEqual(tped_rec.gts[7],
                         "0 0",
                         "TPedReader doesn't work correctly")
        self.assertEqual(tped_rec.gts[8],
                         "G G",
                         "TPedReader doesn't work correctly")

class TestPlinkLib(SafeTester):

    def __init__(self, methodName):
        super(TestPlinkLib, self).__init__(methodName=methodName,
                                           test_module_name=__name__,
                                           )

    def setUp(self):
        pass

    def test_merge_lmiss_map_1(self):
        """ test merging lmiss and map files """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        lmiss_file = join_path(self.data_dir,
                               "input.lmiss")
        map_file = join_path(self.data_dir,
                             "input.map")
        out_file = join_path(self.working_dir,
                             "merged_lmiss_map.snp.info"
                             )
        merge_lmiss_map(lmiss_file, map_file, out_file)
        snp_info_reader = SnpInfoReader(file_name=out_file)
        self.assertEqual(len(list(snp_info_reader)),
                         11,
                         "merge_lmiss_map doesn't work correctly")
        snp_info_reader = SnpInfoReader(file_name=out_file)
        snp_info_rec = snp_info_reader.next()
        self.assertEqual(snp_info_rec.f_miss_a,
                         "F_MISS_A",
                         "merge_lmiss_map doesn't work correctly")
        snp_info_rec = snp_info_reader.next()
        self.assertEqual(snp_info_rec.snp,
                         "rs2795491",
                         "merge_lmiss_map doesn't work correctly")
        snp_info_rec = snp_info_reader.next()
        self.assertEqual(snp_info_rec.pos,
                         "100913376",
                         "merge_lmiss_map doesn't work correctly")
        snp_info_rec = snp_info_reader.next()
        self.assertEqual(snp_info_rec.chrom,
                         "9",
                         "merge_lmiss_map doesn't work correctly")
        snp_info_rec = snp_info_reader.next()
        snp_info_rec = snp_info_reader.next()
        self.assertEqual(snp_info_rec.f_miss_a,
                         0.005,
                         "merge_lmiss_map doesn't work correctly")
        snp_info_rec = snp_info_reader.next()
        self.assertEqual(snp_info_rec.f_miss_u,
                         0,
                         "merge_lmiss_map doesn't work correctly")
        snp_info_rec = snp_info_reader.next()
        self.assertEqual(snp_info_rec.f_miss_a,
                         0.02,
                         "merge_lmiss_map doesn't work correctly")
        self.assertEqual(snp_info_rec.pos,
                         "100930253",
                         "merge_lmiss_map doesn't work correctly")
        snp_info_rec = snp_info_reader.next()
        snp_info_rec = snp_info_reader.next()
        self.assertEqual(snp_info_rec.f_miss_u,
                         0.005,
                         "merge_lmiss_map doesn't work correctly")

    def tearDown(self):
        self.remove_working_dir()
