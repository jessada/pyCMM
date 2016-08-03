import unittest
from os.path import join as join_path
from pycmm.template import SafeTester
from pycmm.cmmlib.fastqlib import FastqReader
from pycmm.cmmlib.fastqlib import Fastq
from pycmm.cmmlib.fastqlib import ENCODING_ILLUMINA_1_5
from pycmm.cmmlib.fastqlib import ENCODING_ILLUMINA_1_8

class TestFastqReader(SafeTester):

    def __init__(self, methodName):
        super(TestFastqReader, self).__init__(methodName=methodName,
                                              test_module_name=__name__,
                                              )

    def setUp(self):
        pass

    def test_reader_1(self):
        """ test basic reading without any other parameters """

        self.init_test(self.current_func_name)
        # Ill_Br concatted fastq (problem fixed in previous GATKBP run)
        input_file1 = join_path(self.data_dir,
                                "input1.fastq.gz")
        fastq_reader = FastqReader(file_name=input_file1)
        self.assertEqual(fastq_reader.next_phred_score().phred_score,
                         "@<@A=DDEAFFGFIJ@;FEHJJJIJJIJ<BFFG@FHH8DB=3BG@3*?B;B09*.8=@4C><@CGGGG;C>GGEEEE;=@3?BC;?A9?;@>C3@A@C::>",
                         "FastqReader doesn't work correctly")
        self.assertEqual(fastq_reader.next_phred_score().phred_score,
                         "@@?DFFADHFHH?GBA@F<EEE;ADCHHC>BAG?FDFGAFHGGGCFAB@<7@B7F3;<DH>)782A>B>?B>1(55:<@>>55@@44:(+2<CCA>:ACB?",
                         "FastqReader doesn't work correctly")
        # Ill_Br raw fastq (problematic in previous GATKBP run)
        input_file2 = join_path(self.data_dir,
                                "input2.fastq.gz")
        fastq_reader = FastqReader(file_name=input_file2)
        self.assertEqual(fastq_reader.next_phred_score().phred_score,
                         "_[_`\ccd`eefehi_Zedgiiihiihi[aeef_eggWca\Raf_RI^aZaOXIMW\_Sb][_bffffZb]ffddddZ\_R^abZ^`X^Z_]bR_`_bYY]",
                         "FastqReader doesn't work correctly")
        self.assertEqual(fastq_reader.next_phred_score().phred_score,
                         "__^cee`cgegg^fa`_e[dddZ`cbggb]a`f^ecef`egfffbe`a_[V_aVeRZ[cg]HVWQ`]a]^a]PGTTY[_]]TT__SSYGJQ[bb`]Y`ba^",
                         "FastqReader doesn't work correctly")
        # Ill_Br raw fastq (supposed to be good)
        input_file3 = join_path(self.data_dir,
                                "input3.fastq.gz")
        fastq_reader = FastqReader(file_name=input_file3)
        self.assertEqual(fastq_reader.next_phred_score().phred_score,
                         "#0;?#2@@?@@@@?@@@@???@?@????@@??@=@@@@????????????????????@@@@@@@>##--#######,,9==????>???##########",
                         "FastqReader doesn't work correctly")
        self.assertEqual(fastq_reader.next_phred_score().phred_score,
                         "BCBF#2ADHHFHHJJJJJIJJJJJJJJIJJJJJJJJJJIEHIJIIJIJJJJJJJJJJIJJJJIJJJJJJGHHHHHFFFFFDFEEEDEDDDDDDDDFDEDD",
                         "FastqReader doesn't work correctly")
        # Axeq chr3_6_14_18
        input_file4 = join_path(self.data_dir,
                                "input4.fastq.gz")
        fastq_reader = FastqReader(file_name=input_file4)
        self.assertEqual(fastq_reader.next_phred_score().phred_score,
                         "@@@DDA+A:ADA<ED@ABE>F43CFEE<<CAC@+AFCAEFCEFCE7??D*00?FF@:((88ACACG7.8=@@@D@)7@EA)7?7=?;?BDCCCAACBB;5;",
                         "FastqReader doesn't work correctly")
        self.assertEqual(fastq_reader.next_phred_score().phred_score,
                         "??@D=D:B:A<<CFEGEEBHGGCHIBHHHA>EG;CGGCHHFDCD@HHGI>CGGIIICGAHGCDGCECGH).@DC>>C;?;==@@CD<A=BC55(5@@8?##",
                         "FastqReader doesn't work correctly")

class TestFastq(SafeTester):

    def __init__(self, methodName):
        super(TestFastq, self).__init__(methodName=methodName,
                                        test_module_name=__name__,
                                        )

    def setUp(self):
        pass

    def test_encoding_1(self):
        """ to check if Phred score encoding version can be identified """

        self.init_test(self.current_func_name)
        # Ill_Br concatted fastq (problem fixed in previous GATKBP run)
        input_file1 = join_path(self.data_dir,
                                "input1.fastq.gz")
        fastq = Fastq(file_name=input_file1)
        self.assertEqual(fastq.encoding,
                         ENCODING_ILLUMINA_1_8,
                         "Encoding of fastq file cannot be guessed correctly")
        # Ill_Br raw fastq (problematic in previous GATKBP run)
        input_file2 = join_path(self.data_dir,
                                "input2.fastq.gz")
        fastq = Fastq(file_name=input_file2)
        self.assertEqual(fastq.encoding,
                         ENCODING_ILLUMINA_1_5,
                         "Encoding of fastq file cannot be guessed correctly")
        # Ill_Br raw fastq (supposed to be good)
        input_file3 = join_path(self.data_dir,
                                "input3.fastq.gz")
        fastq = Fastq(file_name=input_file3)
        self.assertEqual(fastq.encoding,
                         ENCODING_ILLUMINA_1_8,
                         "Encoding of fastq file cannot be guessed correctly")
        # Axeq chr3_6_14_18
        input_file4 = join_path(self.data_dir,
                                "input4.fastq.gz")
        fastq = Fastq(file_name=input_file4)
        self.assertEqual(fastq.encoding,
                         ENCODING_ILLUMINA_1_8,
                         "Encoding of fastq file cannot be guessed correctly")
