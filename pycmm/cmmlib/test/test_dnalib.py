import unittest
from pycmm.template import SafeTester
from pycmm.cmmlib.dnalib import DNARegion


class TestDNARegion(SafeTester):

    def __init__(self, methodName):
        super(TestDNARegion, self).__init__(methodName=methodName,
                                            test_module_name=__name__,
                                            )

    def setUp(self):
        pass

    def test_parse_region_1(self):
        """ test with only chromosome """

        self.init_test(self.current_func_name)
        dr = DNARegion("12")
        self.assertEqual(dr.chrom,
                         "12",
                         "DNARegion cannot correctly identify 'chromosome'")
        self.assertEqual(dr.start_pos,
                         None,
                         "DNARegion cannot correctly identify 'start position'")
        self.assertEqual(dr.end_pos,
                         None,
                         "DNARegion cannot correctly identify 'end position'")

    def test_parse_region_2(self):
        """ test with both chromosome and positions """

        self.init_test(self.current_func_name)
        dr = DNARegion("X:1345-56322")
        self.assertEqual(dr.chrom,
                         "X",
                         "DNARegion cannot correctly identify 'chromosome'")
        self.assertEqual(dr.start_pos,
                         "1345",
                         "DNARegion cannot correctly identify 'start position'")
        self.assertEqual(dr.end_pos,
                         "56322",
                         "DNARegion cannot correctly identify 'end position'")

    def tearDown(self):
        self.remove_working_dir()
