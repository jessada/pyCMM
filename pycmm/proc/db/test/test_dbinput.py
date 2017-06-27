import unittest
from os.path import join as join_path
from pycmm.template import SafeTester
from pycmm.proc.db.dbinput import AVDBReader
from pycmm.proc.db.dbinput import TAVcfInfoReader
from pycmm.proc.db.dbinput import TAVcfGTZReader
from pycmm.proc.db.dbinput import DATA_TYPE_AVDB_CMM_AF
from pycmm.proc.db.dbinput import DATA_TYPE_AVDB_PUBLIC_AF


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
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[0],
                         "1",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_rec[1],
                         14696,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[2],
                         14697,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[5],
                         "0.1",
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
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
                         "_OAF_FAMILIAL_CRCS_WT",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(len(avdb_reader.header_cols),
                         6,
                         "AVDBReader doesn't work correctly")
        avdb_reader = AVDBReader(file_name=input_file, header_exist=True)
        avdb_rec = avdb_reader.record_tuples.next()
        avdb_rec = avdb_reader.record_tuples.next()
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[0],
                         "1",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_rec[1],
                         10394,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[2],
                         10440,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[5],
                         "7",
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
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
                         "_OAF_CHEK2_WT",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_reader.header_cols[9],
                         "_OAF_CHEK2_NA",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(len(avdb_reader.header_cols),
                         14,
                         "AVDBReader doesn't work correctly")
        avdb_reader = AVDBReader(file_name=input_file,
                                 header_exist=True,
                                 data_type=DATA_TYPE_AVDB_CMM_AF)
        avdb_rec = avdb_reader.record_tuples.next()
        avdb_rec = avdb_reader.record_tuples.next()
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[0],
                         "1",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_rec[1],
                         2050195,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[4],
                         "-",
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[5],
                         14,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[3],
                         "-",
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[13],
                         0.8036,
                         "AVDBReader doesn't work correctly")
        self.assertEqual(len(avdb_rec),
                         14,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec,
                         ("1", 2051493, 2051493, "T", "C", 8, 0, 1, 0, 19, 9, 0.3214, 0.1111, 0.0357),
                         "AVDBReader doesn't work correctly")

    def test_reader_4(self):
        """ test reading avdb file with AF has both number an text """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "input.avdb")
        avdb_reader = AVDBReader(file_name=input_file,
                                 header_exist=True,
                                 data_type=DATA_TYPE_AVDB_CMM_AF)
        self.assertEqual(len(list(avdb_reader)),
                         23,
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_reader.header_cols[4],
                         "Alt",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_reader.header_cols[5],
                         "_OAF_EARLYONSET_WT",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(len(avdb_reader.header_cols),
                         14,
                         "AVDBReader doesn't work correctly")
        avdb_reader = AVDBReader(file_name=input_file,
                                 header_exist=True,
                                 data_type=DATA_TYPE_AVDB_CMM_AF)
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec,
                         ("10", 3181107, 3181107, "G", "A", 49, 1, 1, 0, 0, 51, 1.0, 0.0294, 0.0294),
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[0],
                         "10",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_rec[1],
                         8968113,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[4],
                         "T",
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[5],
                         0,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[3],
                         "C",
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[12],
                         "NA",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_rec[13],
                         0.0,
                         "AVDBReader doesn't work correctly")
        self.assertEqual(len(avdb_rec),
                         14,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec,
                         ("10", 8968857, 8968857, "G", "A", 0, 0, 0, 0, 51, 0, 0.0, "NA", 0.0),
                         "AVDBReader doesn't work correctly")

    def test_reader_5(self):
        """ test reading avdb file from gnomad """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "input.avdb")
        avdb_reader = AVDBReader(file_name=input_file,
                                 header_exist=True,
                                 data_type=DATA_TYPE_AVDB_PUBLIC_AF)
        self.assertEqual(len(list(avdb_reader)),
                         1010,
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_reader.header_cols[2],
                         "End",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_reader.header_cols[5],
                         "_gnomAD_genome_ALL",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(len(avdb_reader.header_cols),
                         13,
                         "AVDBReader doesn't work correctly")
        avdb_reader = AVDBReader(file_name=input_file,
                                 header_exist=True,
                                 data_type=DATA_TYPE_AVDB_PUBLIC_AF)
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec,
                         ('10', 896814, 896814, 'A', '-', 0.00003247, 0.0, 0.0, 0.0, 0.0, 0.0, 6.697e-05, 0.0),
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[0],
                         "10",
                         "AVDBReader doesn't work correctly")
        self.assertEqual(avdb_rec[1],
                         896816,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[5],
                         0.00003341,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec[12],
                         0,
                         "AVDBReader doesn't work correctly")
        self.assertEqual(len(avdb_rec),
                         13,
                         "AVDBReader doesn't work correctly")
        avdb_rec = avdb_reader.record_tuples.next()
        self.assertEqual(avdb_rec,
                         ('10', 896823, 896823, 'A', 'C', 3.232e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 6.661e-05, 0.0),
                         "AVDBReader doesn't work correctly")

    def tearDown(self):
        self.remove_working_dir()

class TestTAVcfInfoReader(SafeTester):

    def __init__(self, methodName):
        super(TestTAVcfInfoReader, self).__init__(methodName=methodName,
                                                  test_module_name=__name__,
                                                  )

    def setUp(self):
        pass

    def test_reader_1(self):
        """ test reading standard vcf file from annovar """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "input.vcf.gz")
        vcf_reader = TAVcfInfoReader(file_name=input_file)
        self.assertEqual(len(list(vcf_reader)),
                         12,
                         "TAVcfInfoReader doesn't work correctly")
        self.assertEqual(vcf_reader.header_cols[3],
                         "ALT",
                         "TAVcfInfoReader doesn't work correctly")
        self.assertEqual(vcf_reader.header_cols[5],
                         "_Gene_refGene",
                         "TAVcfInfoReader doesn't work correctly")
        self.assertEqual(vcf_reader.header_cols[9],
                         "_cytoBand",
                         "TAVcfInfoReader doesn't work correctly")
        self.assertEqual(len(vcf_reader.header_cols),
                         11,
                         "TAVcfInfoReader doesn't work correctly")
        vcf_reader = TAVcfInfoReader(file_name=input_file)
        vcf_rec = vcf_reader.record_tuples.next()
        vcf_rec = vcf_reader.record_tuples.next()
        vcf_rec = vcf_reader.record_tuples.next()
        self.assertEqual(vcf_rec[0],
                         "10",
                         "TAVcfInfoReader doesn't work correctly")
        self.assertEqual(vcf_rec[1],
                         89440372,
                         "TAVcfInfoReader doesn't work correctly")
        vcf_rec = vcf_reader.record_tuples.next()
        vcf_rec = vcf_reader.record_tuples.next()
        self.assertEqual(vcf_rec[4],
                         "intronic",
                         "TAVcfInfoReader doesn't work correctly")
        vcf_rec = vcf_reader.record_tuples.next()
        self.assertEqual(vcf_rec[5],
                         "PAPSS2",
                         "TAVcfInfoReader doesn't work correctly")
        vcf_rec = vcf_reader.record_tuples.next()
        self.assertEqual(vcf_rec[3],
                         "A",
                         "TAVcfInfoReader doesn't work correctly")
        vcf_rec = vcf_reader.record_tuples.next()
        self.assertEqual(vcf_rec[9],
                         "10q23.2",
                         "TAVcfInfoReader doesn't work correctly")
        self.assertEqual(len(vcf_rec),
                         11,
                         "TAVcfInfoReader doesn't work correctly")
        vcf_rec = vcf_reader.record_tuples.next()
        self.assertEqual(vcf_rec[:10],
                         ('10', 89442697, 'C', 'T', 'intronic', 'PAPSS2', '', '', '', '10q23.2'),
                         "TAVcfInfoReader doesn't work correctly")

    def tearDown(self):
        self.remove_working_dir()

class TestTAVcfGTZReader(SafeTester):

    def __init__(self, methodName):
        super(TestTAVcfGTZReader, self).__init__(methodName=methodName,
                                                 test_module_name=__name__,
                                                 )

    def setUp(self):
        pass

    def test_reader_1(self):
        """ test reading GTZ from vcf file without ANNOVAR """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "input.vcf.gz")
        vcf_reader = TAVcfGTZReader(file_name=input_file)
        self.assertEqual(len(list(vcf_reader)),
                         6,
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(vcf_reader.header_cols[3],
                         "ALT",
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(vcf_reader.header_cols[4],
                         "QUAL",
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(vcf_reader.header_cols[5],
                         "FILTER",
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(vcf_reader.header_cols[7],
                         "_1416_10D",
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(vcf_reader.header_cols[14],
                         "_2016_18116",
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(len(vcf_reader.header_cols),
                         15,
                         "TAVcfGTZReader doesn't work correctly")
        vcf_reader = TAVcfGTZReader(file_name=input_file)
        vcf_rec = vcf_reader.record_tuples.next()
        self.assertEqual(vcf_rec[0],
                         "10",
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(vcf_rec[1],
                         100000130,
                         "TAVcfGTZReader doesn't work correctly")
        vcf_rec = vcf_reader.record_tuples.next()
        self.assertEqual(vcf_rec[4],
                         308.9,
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(vcf_rec[5],
                         "VQSRTrancheINDEL99.00to99.90",
                         "TAVcfGTZReader doesn't work correctly")
        vcf_rec = vcf_reader.record_tuples.next()
        self.assertEqual(vcf_rec[5],
                         "PASS",
                         "TAVcfGTZReader doesn't work correctly")
        vcf_rec = vcf_reader.record_tuples.next()
        # Plese note that the GTZ test below is not accurate due to 
        # de-multiallelic upstream as some of 'oth' cannot be correctly
        # identified
        self.assertEqual(vcf_rec[6],
                         ".",
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(vcf_rec[7],
                         ".",
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(vcf_rec[8],
                         ".",
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(vcf_rec[9],
                         "het",
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(vcf_rec[10],
                         ".",
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(vcf_rec[11],
                         ".",
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(vcf_rec[12],
                         "oth",
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(vcf_rec[13],
                         ".",
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(vcf_rec[14],
                         ".",
                         "TAVcfGTZReader doesn't work correctly")
        self.assertEqual(len(vcf_rec),
                         15,
                         "TAVcfGTZReader doesn't work correctly")
        vcf_rec = vcf_reader.record_tuples.next()
        self.assertEqual(vcf_rec,
                         ('10', 100000554, 'A', 'ATT', 5045.91, 'PASS', 'het', 'het', '.', '.', 'hom', 'hom', 'het', 'het', 'het'),
                         "TAVcfGTZReader doesn't work correctly")

    def tearDown(self):
        self.remove_working_dir()
