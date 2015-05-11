import filecmp
import unittest
from os.path import join as join_path
from os.path import dirname
from pycmm.template import SafeTester
from pycmm.db.stat import vcf_query


class TestStatFunctions(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            )

    def setUp(self):
        self.module_name = 'stat'

    def __create_db_instance(self):
        return None

    def test_vcf_query_1(self):
        """
        test if vcf_query work correctly if
          - filename
          - chr
          - start pos
          - end  pos
        are specified
        """

        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                            'input.vcf.gz')
        in_chr = "10"
        in_start_pos = "5931000"
        in_end_pos = "5933000"
        records = list(vcf_query(vcf_in_file=in_file,
                                 vcf_chr=in_chr,
                                 vcf_start_pos=in_start_pos,
                                 vcf_end_pos=in_end_pos,
                                 ))
        self.assertEqual(len(records),
                         6,
                         "Invalid number of queried records")
        record = records[0]
        self.assertEqual(record.CHROM,
                         '10',
                         "CHR of vcf query is not correct")
        self.assertEqual(record.POS,
                         5931230,
                         "POS of vcf query is not correct")
        self.assertEqual(record.REF,
                         'C',
                         "REF of vcf query is not correct")
        self.assertEqual(record.ALT,
                         ['T'],
                         "ALT of vcf query is not correct")
        self.assertEqual(record.genotype("Sample38").sample,
                         "Sample38",
                         "The name of samples[0] of vcf query is not correct")
        self.assertEqual(record.genotype("Sample38").gt_alleles,
                         ['0', '0'],
                         "genotype data of vcf query is not correct")
        self.assertEqual(record.genotype("sample23").gt_alleles,
                         [None, None],
                         "genotype data of vcf query is not correct")
        record = records[4]
        self.assertEqual(record.CHROM,
                         '10',
                         "CHR of vcf query is not correct")
        self.assertEqual(record.POS,
                         5932475,
                         "POS of vcf query is not correct")
        self.assertEqual(record.REF,
                         'T',
                         "REF of vcf query is not correct")
        self.assertEqual(record.ALT,
                         ['C'],
                         "ALT of vcf query is not correct")
        self.assertEqual(record.genotype("Sample38").gt_alleles,
                         ['0', '0'],
                         "genotype data of vcf query is not correct")
        self.assertEqual(record.genotype("S-43").gt_alleles,
                         ['1', '1'],
                         "genotype data of vcf query is not correct")
        self.assertEqual(record.genotype("6-s").gt_alleles,
                         [None, None],
                         "genotype data of vcf query is not correct")
        self.assertEqual(record.genotype("47s").gt_alleles,
                         ['0', '1'],
                         "genotype data of vcf query is not correct")

    def test_vcf_query_2(self):
        """
        test if vcf_query work correctly if
          - filename
          - chr
        are specified
        """

        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                            'input.vcf.gz')
        in_chr = "3"
        records = list(vcf_query(vcf_in_file=in_file,
                                 vcf_chr=in_chr,
                                 ))
        self.assertEqual(len(records),
                         13,
                         "Invalid number of queried records")
        record = records[6]
        self.assertEqual(record.CHROM,
                         '3',
                         "CHR of vcf query is not correct")
        self.assertEqual(record.POS,
                         142542683,
                         "POS of vcf query is not correct")
        self.assertEqual(record.REF,
                         'T',
                         "REF of vcf query is not correct")
        self.assertEqual(record.ALT,
                         ['C'],
                         "ALT of vcf query is not correct")

    def test_vcf_query_3(self):
        """
        test if vcf_query work correctly if filename is specified
        """

        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                            'input.vcf.gz')
        records = list(vcf_query(vcf_in_file=in_file,
                                 ))
        self.assertEqual(len(records),
                         65,
                         "Invalid number of queried records")
        record = records[46]
        self.assertEqual(record.CHROM,
                         '19',
                         "CHR of vcf query is not correct")
        self.assertEqual(record.POS,
                         57740608,
                         "POS of vcf query is not correct")
        self.assertEqual(record.REF,
                         'G',
                         "REF of vcf query is not correct")
        self.assertEqual(record.ALT,
                         ['A'],
                         "ALT of vcf query is not correct")
        self.assertEqual(record.genotype("Sample38").gt_alleles,
                         ['1', '1'],
                         "genotype data of vcf query is not correct")

    @unittest.skipUnless(settings.FULL_SYSTEM_TEST, "taking too long time to test")
    def test_vcf_query_4(self):
        """
        test if vcf_query work correctly with very big vcf file
        """

        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                            'input.vcf.gz')
        records = vcf_query(vcf_in_file=in_file)
        rec_count = 0
        for record in records:
            rec_count +=1
        print rec_count
        self.assertEqual(rec_count,
                         65,
                         "Invalid number of queried records")

