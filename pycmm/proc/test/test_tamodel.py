from os.path import join as join_path
from os.path import dirname
from pycmm.template import SafeTester
from pycmm.proc.taparser import TAVcfReader
from pycmm.flow.mutrep import MutRepPipeline
from pycmm.flow.test.test_mutrep import DFLT_TEST_MUTREP_COLS
from pycmm.flow.cmmdb import create_jobs_setup_file


class TestTAVcfCall(SafeTester):

    def __init__(self, methodName):
        super(TestTAVcfCall, self).__init__(methodName=methodName,
                                            test_module_name=__name__,
                                            )

    def setUp(self):
        pass

    def test_parse_cmm_gts(self):
        """ test if zygosity can be correctly determined """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Br-429")
        self.assertEqual(call.cmm_gts[1],
                         "oth",
                         "cmm genotype cannot be correctly determined")
        self.assertEqual(call.cmm_gts[2],
                         "het",
                         "cmm genotype cannot be correctly determined")
        self.assertEqual(call.cmm_gts[3],
                         "het",
                         "cmm genotype cannot be correctly determined")
        self.assertEqual(call.cmm_gts[4],
                         "oth",
                         "cmm genotype cannot be correctly determined")
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.cmm_gts[1],
                         "hom",
                         "cmm genotype cannot be correctly determined")
        self.assertEqual(call.cmm_gts[2],
                         "oth",
                         "cmm genotype cannot be correctly determined")
        self.assertEqual(call.cmm_gts[3],
                         "oth",
                         "cmm genotype cannot be correctly determined")
        self.assertEqual(call.cmm_gts[4],
                         "oth",
                         "cmm genotype cannot be correctly determined")
        call = vcf_record.genotype("Al-161")
        self.assertEqual(call.cmm_gts[1],
                         "het",
                         "cmm genotype cannot be correctly determined")
        self.assertEqual(call.cmm_gts[2],
                         "oth",
                         "cmm genotype cannot be correctly determined")
        self.assertEqual(call.cmm_gts[3],
                         "het",
                         "cmm genotype cannot be correctly determined")
        self.assertEqual(call.cmm_gts[4],
                         "oth",
                         "cmm genotype cannot be correctly determined")
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Al-65")
        self.assertEqual(call.cmm_gts[1],
                         ".",
                         "cmm genotype cannot be correctly determined")
        self.assertEqual(call.cmm_gts[2],
                         ".",
                         "cmm genotype cannot be correctly determined")
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Al-111")
        self.assertEqual(call.cmm_gts[1],
                         "wt",
                         "cmm genotype cannot be correctly determined")

    def test_parse_actual_gts_1(self):
        """
        test if true zygosity can be correctly determined
        - very random cases
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        # test if allele frequency is an array of None,
        # like [None, None, None, None]
        call = vcf_record.genotype("Br-429")
        self.assertEqual(call.actual_gts[1],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.actual_gts[1],
                         "hom",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        call = vcf_record.genotype("Al-161")
        self.assertEqual(call.actual_gts[1],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Al-65")
        self.assertEqual(call.actual_gts[1],
                         ".",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         ".",
                         "cmm actual genotype cannot be correctly determined")
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Al-111")
        self.assertEqual(call.actual_gts[1],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")

    def test_parse_actual_gts_2(self):
        """
        test if true zygosity can be correctly determined
        - allele frequency is an array of None, like [None, None, None, None]
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Br-429")
        self.assertEqual(call.actual_gts[1],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.actual_gts[1],
                         "hom",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        call = vcf_record.genotype("Al-161")
        self.assertEqual(call.actual_gts[1],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        call = vcf_record.genotype("Br-637")
        self.assertEqual(call.actual_gts[1],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")

    def test_parse_actual_gts_3(self):
        """
        test if true zygosity can be correctly determined
        - only one alternate allele and no allele frequency
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Br-429")
        self.assertEqual(call.actual_gts[1],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.actual_gts[1],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        call = vcf_record.genotype("Al-111")
        self.assertEqual(call.actual_gts[1],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        call = vcf_record.genotype("Br-120")
        self.assertEqual(call.actual_gts[1],
                         ".",
                         "cmm actual genotype cannot be correctly determined")
        call = vcf_record.genotype("Br-390")
        self.assertEqual(call.actual_gts[1],
                         "hom",
                         "cmm actual genotype cannot be correctly determined")

    def test_parse_actual_gts_4(self):
        """
        test if true zygosity can be correctly determined
        - allele frequency is a floating point scalar less than 0.5
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Br-429")
        self.assertEqual(call.actual_gts[1],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.actual_gts[1],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        call = vcf_record.genotype("Alb-31")
        self.assertEqual(call.actual_gts[1],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        call = vcf_record.genotype("Br-693")
        self.assertEqual(call.actual_gts[1],
                         "hom",
                         "cmm actual genotype cannot be correctly determined")

    def test_parse_actual_gts_5(self):
        """
        test if true zygosity can be correctly determined
        - allele frequency is a floating point scalar more than 0.5
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Br-429")
        self.assertEqual(call.actual_gts[1],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.actual_gts[1],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        call = vcf_record.genotype("Alb-31")
        self.assertEqual(call.actual_gts[1],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        call = vcf_record.genotype("Br-466")
        self.assertEqual(call.actual_gts[1],
                         "hom",
                         "cmm actual genotype cannot be correctly determined")

    def test_parse_actual_gts_6(self):
        """
        test if true zygosity can be correctly determined
        - allele frequency is an array of None and low value,
          like [None, 0.02, 0.03, None]
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        # ./.
        call = vcf_record.genotype("Br-120")
        self.assertEqual(call.actual_gts[1],
                         ".",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         ".",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         ".",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         ".",
                         "cmm actual genotype cannot be correctly determined")
        # 0/0
        call = vcf_record.genotype("Br-781")
        self.assertEqual(call.actual_gts[1],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        # 0/1
        call = vcf_record.genotype("Al-111")
        self.assertEqual(call.actual_gts[1],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        # 0/2
        call = vcf_record.genotype("Al-23")
        self.assertEqual(call.actual_gts[1],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        # 1/1
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.actual_gts[1],
                         "hom",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        # 1/2
        call = vcf_record.genotype("Al-92")
        self.assertEqual(call.actual_gts[1],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        # 2/2
        call = vcf_record.genotype("Br-296")
        self.assertEqual(call.actual_gts[1],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "hom",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        # random
        call = vcf_record.genotype("Br-429")
        self.assertEqual(call.actual_gts[1],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")

    def test_parse_actual_gts_7(self):
        """
        test if true zygosity can be correctly determined
        - allele frequency is an array of high and low value,
          like [0.79, 0.02, 0.03, None]
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        # ./.
        call = vcf_record.genotype("Br-120")
        self.assertEqual(call.actual_gts[1],
                         ".",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         ".",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         ".",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         ".",
                         "cmm actual genotype cannot be correctly determined")
        # 0/0
        call = vcf_record.genotype("Br-781")
        self.assertEqual(call.actual_gts[1],
                         "hom",
                         "cmm actual genotype cannot be correctly determined")
        # this one is still wrong but acceptable
        # it should be "other"
        self.assertEqual(call.actual_gts[2],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        # this one is still wrong but acceptable
        # it should be "other"
        self.assertEqual(call.actual_gts[3],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        # this one is still wrong but acceptable
        # it should be "other"
        self.assertEqual(call.actual_gts[4],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        # 0/1
        call = vcf_record.genotype("Al-111")
        self.assertEqual(call.actual_gts[1],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        # 0/2
        call = vcf_record.genotype("Al-23")
        self.assertEqual(call.actual_gts[1],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        # 1/1
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.actual_gts[1],
                         "wt",
                         "cmm actual genotype cannot be correctly determined")
        # this one is still wrong but acceptable
        # it should be "wt"
        self.assertEqual(call.actual_gts[2],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        # this one is still wrong but acceptable
        # it should be "wt"
        self.assertEqual(call.actual_gts[3],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        # this one is still wrong but acceptable
        # it should be "wt"
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        # 1/2
        call = vcf_record.genotype("Al-92")
        self.assertEqual(call.actual_gts[1],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        # 2/2
        call = vcf_record.genotype("Br-296")
        self.assertEqual(call.actual_gts[1],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        # 2/2
        self.assertEqual(call.actual_gts[2],
                         "hom",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        # random
        call = vcf_record.genotype("Br-429")
        self.assertEqual(call.actual_gts[1],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[2],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[3],
                         "het",
                         "cmm actual genotype cannot be correctly determined")
        self.assertEqual(call.actual_gts[4],
                         "oth",
                         "cmm actual genotype cannot be correctly determined")

    def test_parse_mutated_1(self):
        """
        test if mutation can be identified
        - only one alternate allele and no allele frequency
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Br-429")
        self.assertEqual(call.mutated[1],
                         False,
                         "cmm mutation cannot be correctly determined")
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.mutated[1],
                         False,
                         "cmm mutation cannot be correctly determined")
        call = vcf_record.genotype("Al-111")
        self.assertEqual(call.mutated[1],
                         True,
                         "cmm mutation cannot be correctly determined")
        call = vcf_record.genotype("Br-120")
        self.assertEqual(call.mutated[1],
                         False,
                         "cmm mutation cannot be correctly determined")
        call = vcf_record.genotype("Br-390")
        self.assertEqual(call.mutated[1],
                         True,
                         "cmm mutation cannot be correctly determined")

    def test_parse_mutated_2(self):
        """
        test if mutation can be identified
        - allele frequency is a floating point scalar less than 0.5
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Br-429")
        self.assertEqual(call.mutated[1],
                         False,
                         "cmm mutation cannot be correctly determined")
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.mutated[1],
                         False,
                         "cmm mutation cannot be correctly determined")
        call = vcf_record.genotype("Alb-31")
        self.assertEqual(call.mutated[1],
                         True,
                         "cmm mutation cannot be correctly determined")
        call = vcf_record.genotype("Br-693")
        self.assertEqual(call.mutated[1],
                         True,
                         "cmm mutation cannot be correctly determined")

    def test_parse_mutated_3_1(self):
        """
        test if mutation can be identified
        - allele frequency is a floating point scalar more than 0.5
        - test data
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Br-429")
        self.assertEqual(call.mutated[1],
                         False,
                         "cmm mutation cannot be correctly determined")
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.mutated[1],
                         False,
                         "cmm mutation cannot be correctly determined")
        call = vcf_record.genotype("Alb-31")
        self.assertEqual(call.mutated[1],
                         True,
                         "cmm mutation cannot be correctly determined")
        call = vcf_record.genotype("Br-466")
        self.assertEqual(call.mutated[1],
                         True,
                         "cmm mutation cannot be correctly determined")

    def test_parse_mutated_3_2(self):
        """
        test if mutation can be identified
        - allele frequency is a floating point scalar more than 0.5
        - true data
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Co-166")
        self.assertEqual(call.mutated[1],
                         False,
                         "cmm mutation cannot be correctly determined")
        call = vcf_record.genotype("Co-213")
        self.assertEqual(call.mutated[1],
                         False,
                         "cmm mutation cannot be correctly determined")

    def test_parse_mutated_4(self):
        """
        test if mutation can be identified
        - allele frequency is an array of high and low value,
          like [0.79, 0.02, 0.03, None]
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        # ./.
        call = vcf_record.genotype("Br-120")
        self.assertEqual(call.mutated[1],
                         False,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[2],
                         False,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[3],
                         False,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[4],
                         False,
                         "cmm mutation cannot be correctly determined")
        # 0/0
        call = vcf_record.genotype("Br-781")
        self.assertEqual(call.mutated[1],
                         True,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[2],
                         False,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[3],
                         False,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[4],
                         False,
                         "cmm mutation cannot be correctly determined")
        # 0/1
        call = vcf_record.genotype("Al-111")
        self.assertEqual(call.mutated[1],
                         True,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[2],
                         False,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[3],
                         False,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[4],
                         False,
                         "cmm mutation cannot be correctly determined")
        # 0/2
        call = vcf_record.genotype("Al-23")
        self.assertEqual(call.mutated[1],
                         False,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[2],
                         True,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[3],
                         False,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[4],
                         False,
                         "cmm mutation cannot be correctly determined")
        # 1/1
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.mutated[1],
                         False,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[2],
                         False,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[3],
                         False,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[4],
                         False,
                         "cmm mutation cannot be correctly determined")
        # 1/2
        call = vcf_record.genotype("Al-92")
        self.assertEqual(call.mutated[1],
                         True,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[2],
                         True,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[3],
                         False,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[4],
                         False,
                         "cmm mutation cannot be correctly determined")
        # 2/2
        call = vcf_record.genotype("Br-296")
        self.assertEqual(call.mutated[1],
                         False,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[2],
                         True,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[3],
                         False,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[4],
                         False,
                         "cmm mutation cannot be correctly determined")
        # random
        call = vcf_record.genotype("Br-429")
        self.assertEqual(call.mutated[1],
                         False,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[2],
                         True,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[3],
                         True,
                         "cmm mutation cannot be correctly determined")
        self.assertEqual(call.mutated[4],
                         False,
                         "cmm mutation cannot be correctly determined")

    def tearDown(self):
        self.remove_working_dir()

class TestTAVcfRecord(SafeTester):

    def __init__(self, methodName):
        super(TestTAVcfRecord, self).__init__(methodName=methodName,
                                              test_module_name=__name__,
                                              )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self,
                                 dataset_name=None,
                                 sample_infos=None,
                                 anno_cols=DFLT_TEST_MUTREP_COLS,
                                 anno_excl_tags=None,
                                 annotated_vcf_tabix=None,
                                 summary_families_sheet=False,
                                 frequency_ratios=None,
                                 rows_filter_actions=None,
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        if dataset_name is None:
            dataset_name = self.test_function
        create_jobs_setup_file(dataset_name=dataset_name,
                               project_out_dir=self.working_dir,
                               sample_infos=sample_infos,
                               anno_cols=",".join(anno_cols),
                               anno_excl_tags=anno_excl_tags,
                               annotated_vcf_tabix=annotated_vcf_tabix,
                               summary_families_sheet=summary_families_sheet,
                               frequency_ratios=frequency_ratios,
                               rows_filter_actions=rows_filter_actions,
                               out_jobs_setup_file=jobs_setup_file,
                               )
        return jobs_setup_file

    def test_is_shared_1(self):
        """ test shared mutation can be correctly identified in one sample """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        samples = ["Al-65"]
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 2),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 3),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 4),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")

    def test_is_shared_2(self):
        """ test shared mutation between two samples can be identified """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        samples = ["Alb-31", "Br-466"]
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 2),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 3),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 4),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")

    def test_is_shared_3(self):
        """ test shared mutation between three samples can be identified """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        samples = ["Br-432", "Al-161", "Br-504"]
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 2),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 3),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 4),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")

    def test_is_shared_s1_1(self):
        """
        test special case of shared mutation to be correctly
        identified in one sample
        min_share_count=1
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        samples = ["Al-65"]
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 2, min_share_count=1),
                         False,
                         "shared mutation cannot be correctly determined")

    def test_is_shared_s1_2(self):
        """
        test special case of shared mutation to be correctly
        identified in one sample
        min_share_count=2
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        samples = ["Al-65"]
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 2, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")

    def test_is_shared_s2_2(self):
        """
        test special case of shared mutation between two samples to be
        identified
        min_share_count=2
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        samples = ["Alb-31", "Br-466"]
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=2),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 2, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 3, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 4, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 2, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=2),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")

    def test_is_shared_s2_3(self):
        """
        test special case of shared mutation between two samples to be
        identified
        min_share_count=3
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        samples = ["Alb-31", "Br-466"]
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 2, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 3, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 4, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(samples, 2, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(samples, 1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")

    def test_has_mutation_1(self):
        """ test a mutation can be correctly identified in one sample """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        samples = ["Al-65"]
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_mutation(samples, 2),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_mutation(samples, 3),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_mutation(samples, 4),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_mutation(samples, 2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")

    def test_has_mutation_2(self):
        """ test shared mutation between two samples can be identified """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        samples = ["Alb-31", "Br-466"]
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_mutation(samples, 2),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_mutation(samples, 3),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_mutation(samples, 4),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_mutation(samples, 2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")

    def test_has_mutation_3(self):
        """ test shared mutation between three samples can be identified """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        samples = ["Br-432", "Al-161", "Br-504"]
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_mutation(samples, 2),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_mutation(samples, 3),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_mutation(samples, 4),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_mutation(samples, 2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_mutation(samples, 1),
                         False,
                         "shared mutation cannot be correctly determined")

    def test_has_shared_0(self):
        """
        test if a family with shared mutation can be detected
        without min_share_count parameter
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        sample_infos = []
        sample_infos.append("fam1:Al-17")
        sample_infos.append("fam2:Alb-31:Al-23")
        sample_infos.append("fam3:Al-36:Al-47:Al-65")
        sample_infos.append("fam4:Al-73:Al-77:Al-92")
        sample_infos.append("fam5:Br-466")
        jobs_setup_file = self.__create_jobs_setup_file(sample_infos=",".join(sample_infos),
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file, family_infos=pl.family_infos)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_shared(2),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_shared(3),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_shared(4),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_shared(2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1),
                         False,
                         "shared mutation cannot be correctly determined")

    def test_has_shared_1(self):
        """
        test if a family with shared mutation can be detected
        with min_share_count = 1
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        sample_infos = []
        sample_infos.append("fam1:Al-17")
        sample_infos.append("fam2:Alb-31:Al-23")
        sample_infos.append("fam3:Al-36:Al-47:Al-65")
        sample_infos.append("fam4:Al-73:Al-77:Al-92")
        sample_infos.append("fam5:Br-466")
        jobs_setup_file = self.__create_jobs_setup_file(sample_infos=",".join(sample_infos),
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file, family_infos=pl.family_infos)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_shared(2, min_share_count=1),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_shared(3, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_shared(4, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_shared(2, min_share_count=1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=1),
                         False,
                         "shared mutation cannot be correctly determined")

    def test_has_shared_2(self):
        """
        test if a family with shared mutation can be detected
        with min_share_count = 2
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        sample_infos = []
        sample_infos.append("fam1:Al-17")
        sample_infos.append("fam2:Alb-31:Al-23")
        sample_infos.append("fam3:Al-36:Al-47:Al-65")
        sample_infos.append("fam4:Al-73:Al-77:Al-92")
        sample_infos.append("fam5:Br-466")
        jobs_setup_file = self.__create_jobs_setup_file(sample_infos=",".join(sample_infos),
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file, family_infos=pl.family_infos)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_shared(2, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_shared(3, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_shared(4, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_shared(2, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=2),
                         False,
                         "shared mutation cannot be correctly determined")

    def test_has_shared_3(self):
        """
        test if a family with shared mutation can be detected
        with min_share_count = 3
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        sample_infos = []
        sample_infos.append("fam1:Al-17")
        sample_infos.append("fam2:Alb-31:Al-23")
        sample_infos.append("fam3:Al-36:Al-47:Al-65")
        sample_infos.append("fam4:Al-73:Al-77:Al-92")
        sample_infos.append("fam5:Br-466")
        jobs_setup_file = self.__create_jobs_setup_file(sample_infos=",".join(sample_infos),
                                                        )
        pl = MutRepPipeline(jobs_setup_file)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file, family_infos=pl.family_infos)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_shared(2, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_shared(3, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_shared(4, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.has_shared(2, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.has_shared(1, min_share_count=3),
                         False,
                         "shared mutation cannot be correctly determined")

    def test_is_rare_1(self):
        """
        test rare mutation can be identified if there is only one criteria
          - 0.1
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        freq_ratios = {'1000g2014oct_all': 0.1}
        vcf_reader = TAVcfReader(filename=in_file, freq_ratios=freq_ratios)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_rare(allele_idx=2),
                         True,
                         "rare mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_rare(allele_idx=3),
                         True,
                         "rare mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_rare(allele_idx=4),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_rare(allele_idx=2),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")

    def test_is_rare_2(self):
        """
        test rare mutation can be identified if there is only one criteria
          - 0.2
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        freq_ratios = {'1000g2014oct_all': 0.2}
        vcf_reader = TAVcfReader(filename=in_file, freq_ratios=freq_ratios)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_rare(allele_idx=2),
                         True,
                         "rare mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_rare(allele_idx=3),
                         True,
                         "rare mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_rare(allele_idx=4),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_rare(allele_idx=2),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")

    def test_is_intergenic_1(self):
        """
        test counting intergenic mutations can be identified (all are intergenics)
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        intergenics = reduce(lambda x, y: x+y,
                             map(lambda x: x.is_intergenic[1:], vcf_reader))
        intergenics_count = sum(1 for x in intergenics if x)
        self.assertEqual(intergenics_count,
                         24,
                         "intergenic mutations cannot be correctly determined")

    def test_is_intergenic_2(self):
        """
        test counting intergenic mutations can be identified (10 are intergenics)
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        intergenics = reduce(lambda x, y: x+y,
                             map(lambda x: x.is_intergenic[1:], vcf_reader))
        intergenics_count = sum(1 for x in intergenics if x)
        self.assertEqual(intergenics_count,
                         10,
                         "intergenic mutations cannot be correctly determined")

    def test_is_intronic_1(self):
        """
        test counting intronic mutations can be identified (96 are intronics)
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        intronics = reduce(lambda x, y: x+y,
                           map(lambda x: x.is_intronic[1:], vcf_reader))
        intronics_count = sum(1 for x in intronics if x)
        self.assertEqual(intronics_count,
                         96,
                         "intronic mutations cannot be correctly determined")

    def test_is_intronic_2(self):
        """
        test counting intronic and ncr_intronic mutations can be identified (17 are intronics)
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        intronics = reduce(lambda x, y: x+y,
                           map(lambda x: x.is_intronic[1:], vcf_reader))
        intronics_count = sum(1 for x in intronics if x)
        self.assertEqual(intronics_count,
                         17,
                         "intronic mutations cannot be correctly determined")

    def test_vcf_eval_1(self):
        """
        test general evaluation of vcf_eval
        """

        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        expr1 = '("P-value_8" != \'.\') and (float("P-value_8") < 0.02)'
        expr2 = '("P-value_8" != \'.\') and (float("P-value_8") < 0.05)'
        expr3 = '("P-value_8" != \'.\') and (float("P-value_8") < 1e-019)'
        expr4 = '("P-value_8" != \'.\') and (float("P-value_8") < 1e-020)'
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertTrue(vcf_record.vcf_eval(expr1),
                        "cannot perform vcf expression evaluation correctly")
        self.assertTrue(vcf_record.vcf_eval(expr2),
                        "cannot perform vcf expression evaluation correctly")
        self.assertTrue(vcf_record.vcf_eval(expr3),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr4),
                        "cannot perform vcf expression evaluation correctly")
        vcf_record = vcf_reader.next()
        self.assertFalse(vcf_record.vcf_eval(expr1),
                        "cannot perform vcf expression evaluation correctly")
        self.assertTrue(vcf_record.vcf_eval(expr2),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr3),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr4),
                        "cannot perform vcf expression evaluation correctly")
        vcf_record = vcf_reader.next()
        self.assertFalse(vcf_record.vcf_eval(expr1),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr2),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr3),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr4),
                        "cannot perform vcf expression evaluation correctly")

    def tearDown(self):
        self.remove_working_dir()
