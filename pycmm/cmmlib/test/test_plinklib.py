import unittest
from os.path import join as join_path
from pycmm.template import SafeTester
from pycmm.cmmlib.plinklib import HapAssocReader
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

class TestPlinkLib(SafeTester):

    def __init__(self, methodName):
        super(TestPlinkLib, self).__init__(methodName=methodName,
                                           test_module_name=__name__,
                                           )

    def setUp(self):
        pass

    def test_merge_hap_assocs_1(self):
        """ test merging assoc.hap files """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        input_files = []
        input_files.append(join_path(self.data_dir,
                                     "input1.assoc.hap"))
        input_files.append(join_path(self.data_dir,
                                     "input2.assoc.hap"))
        locus_prefixs = ["W1", "W3"]
        out_file = join_path(self.working_dir,
                             "merged_hap_assoc.assoc.hap"
                             )
        merge_hap_assocs(input_files, out_file, locus_prefixs=locus_prefixs)
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
