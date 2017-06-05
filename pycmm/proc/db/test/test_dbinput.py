import unittest
from os.path import join as join_path
from pycmm.template import SafeTester
from pycmm.proc.db.dbinput import AVDBReader
from pycmm.proc.db.dbinput import TAVcfReader
#from pycmm.utils import is_number


class TestAVDBReader(SafeTester):

    def __init__(self, methodName):
        super(TestAVDBReader, self).__init__(methodName=methodName,
                                             test_module_name=__name__,
                                             )

    def setUp(self):
        pass

    def test_reader_1(self):
        """ test reading avdb file without header """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "hg19_avsift.anno.txt")
        avdb_reader = AVDBReader(file_name=input_file)
        self.assertEqual(len(list(avdb_reader)),
                         10,
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_reader.header_cols[3],
                         "Ref",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_reader.header_cols[5],
                         "avsift.anno",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(len(avdb_reader.header_cols),
                         6,
                         "AVDBReader doesn't work correctly")
        avdb_reader = AVDBReader(file_name=input_file)
        avdb_rec = avdb_reader.next()
        self.assertEqual(avdb_rec[0],
                         "1",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_rec[1],
                         14696,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.next()
        avdb_rec = avdb_reader.next()
        self.assertEqual(avdb_rec[2],
                         14697,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.next()
        self.assertEqual(avdb_rec[5],
                         "0.1",
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.next()
        self.assertEqual(avdb_rec,
                         ("1", 14697, 14697, "T", "G", "0.86"),
                         "AVDBReader doesn't work correctly")

    def test_reader_2(self):
        """ test reading avdb file with header and one anno column """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "input.avdb")
        avdb_reader = AVDBReader(file_name=input_file, header_exist=True)
        self.assertEqual(len(list(avdb_reader)),
                         9,
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_reader.header_cols[4],
                         "Alt",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_reader.header_cols[5],
                         "OAF_FAMILIAL_CRCS_WT",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(len(avdb_reader.header_cols),
                         6,
                         "AVDBReader doesn't work correctly")
        avdb_reader = AVDBReader(file_name=input_file, header_exist=True)
        avdb_rec = avdb_reader.next()
        avdb_rec = avdb_reader.next()
        avdb_rec = avdb_reader.next()
        self.assertEqual(avdb_rec[0],
                         "1",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_rec[1],
                         10394,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.next()
        avdb_rec = avdb_reader.next()
        self.assertEqual(avdb_rec[2],
                         10440,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.next()
        self.assertEqual(avdb_rec[5],
                         "7",
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.next()
        self.assertEqual(avdb_rec,
                         ("1", 10440, 10440, "C", "-", "7"),
                         "AVDBReader doesn't work correctly")

    def test_reader_3(self):
        """ test reading avdb file with header and nine anno columns """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "input.avdb")
        avdb_reader = AVDBReader(file_name=input_file, header_exist=True)
        self.assertEqual(len(list(avdb_reader)),
                         11,
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_reader.header_cols[4],
                         "Alt",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_reader.header_cols[5],
                         "OAF_CHEK2_WT",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_reader.header_cols[9],
                         "OAF_CHEK2_NA",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(len(avdb_reader.header_cols),
                         14,
                         "AVDBReader doesn't work correctly")
        avdb_reader = AVDBReader(file_name=input_file, header_exist=True)
        avdb_rec = avdb_reader.next()
        avdb_rec = avdb_reader.next()
        avdb_rec = avdb_reader.next()
        self.assertEqual(avdb_rec[0],
                         "1",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_rec[1],
                         2050195,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.next()
        avdb_rec = avdb_reader.next()
        self.assertEqual(avdb_rec[4],
                         "-",
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.next()
        self.assertEqual(avdb_rec[5],
                         "14",
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.next()
        self.assertEqual(avdb_rec[3],
                         "-",
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.next()
        self.assertEqual(avdb_rec[13],
                         "0.8036",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(len(avdb_rec),
                         14,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.next()
        self.assertEqual(avdb_rec,
                         ("1", 2051493, 2051493, "T", "C", "8", "0", "1", "0", "19", "9", "0.3214", "0.1111", "0.0357"),
                         "AVDBReader doesn't work correctly")

    def tearDown(self):
        self.remove_working_dir()

class TestTAVcfReader(SafeTester):

    def __init__(self, methodName):
        super(TestTAVcfReader, self).__init__(methodName=methodName,
                                             test_module_name=__name__,
                                             )

    def setUp(self):
        pass

    def test_reader_1(self):
        """ test reading standard vcf file from annovar """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "input.vcf.gz")
        vcf_reader = TAVcfReader(file_name=input_file)
        self.assertEqual(len(list(vcf_reader)),
                         12,
                         "TAVcfReader doesn't work correctly")
        self.assertEqual(vcf_reader.header_cols[3],
                         "ALT",
                         "TAVcfReader doesn't work correctly")
        self.assertEqual(vcf_reader.header_cols[5],
                         "Gene.refGene",
                         "TAVcfReader doesn't work correctly")
        self.assertEqual(vcf_reader.header_cols[9],
                         "cytoBand",
                         "TAVcfReader doesn't work correctly")
        self.assertEqual(len(vcf_reader.header_cols),
                         11,
                         "TAVcfReader doesn't work correctly")
        vcf_reader = TAVcfReader(file_name=input_file)
        vcf_rec = vcf_reader.next()
        vcf_rec = vcf_reader.next()
        vcf_rec = vcf_reader.next()
        self.assertEqual(vcf_rec[0],
                         "10",
                         "TAVcfReader doesn't work correctly")
        self.assertEqual(vcf_rec[1],
                         89440372,
                         "TAVcfReader doesn't work correctly")
        vcf_rec = vcf_reader.next()
        vcf_rec = vcf_reader.next()
        self.assertEqual(vcf_rec[4],
                         "intronic",
                         "TAVcfReader doesn't work correctly")
        vcf_rec = vcf_reader.next()
        self.assertEqual(vcf_rec[5],
                         "PAPSS2",
                         "TAVcfReader doesn't work correctly")
        vcf_rec = vcf_reader.next()
        self.assertEqual(vcf_rec[3],
                         "A",
                         "TAVcfReader doesn't work correctly")
        vcf_rec = vcf_reader.next()
        self.assertEqual(vcf_rec[9],
                         "10q23.2",
                         "TAVcfReader doesn't work correctly")
        self.assertEqual(len(vcf_rec),
                         11,
                         "TAVcfReader doesn't work correctly")
        vcf_rec = vcf_reader.next()
        self.assertEqual(vcf_rec[:10],
                         ('10', 89442697, 'C', 'T', 'intronic', 'PAPSS2', '.', '.', '.', '10q23.2'),
                         "TAVcfReader doesn't work correctly")


    def tearDown(self):
        self.remove_working_dir()
