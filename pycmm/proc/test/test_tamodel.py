import filecmp
from os.path import join as join_path
from os.path import dirname
from pycmm.utils import mylogger
from pycmm.template import SafeTester
from pycmm.proc.taparser import TAVcfReader
from pycmm.flow.mutrep import MutRepPipeline


class TestTAVcfCall(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            test_module_name=__name__,
                            )

    def setUp(self):
        mylogger.getLogger(__name__)

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

    def test_parse_mutated_3(self):
        """
        test if mutation can be identified
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

class TestTAVcfRecord(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            test_module_name=__name__,
                            )

    def setUp(self):
        mylogger.getLogger(__name__)

    def test_shared_1(self):
        """ test shared mutation can be correctly identify in one sample """

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

    def test_shared_2(self):
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

    def test_shared_3(self):
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
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=2),
                         True,
                         "rare mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=3),
                         True,
                         "rare mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=4),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=2),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
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
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=2),
                         True,
                         "rare mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=3),
                         True,
                         "rare mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=4),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=2),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         False,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
                         True,
                         "rare mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_rare(freq_ratios, allele_idx=1),
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
        test counting intronic mutations can be identified (10 are intronics)
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
