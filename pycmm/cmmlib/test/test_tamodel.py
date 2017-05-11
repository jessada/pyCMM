from os.path import join as join_path
from os.path import dirname
from pycmm.template import SafeTester
from pycmm.settings import KG2014OCT_ALL_COL_NAME
from pycmm.settings import EST_KVOT_EARLYONSET_VS_BRC_COL_NAME
from pycmm.settings import EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME
from pycmm.settings import EST_KVOT_EARLYONSET_VS_KG_EUR_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_EXP_SYN_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_N_SYN_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_SYN_Z_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_EXP_MIS_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_N_MIS_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_MIS_Z_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_EXP_LOF_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_N_LOF_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_PLI_COL_NAME
from pycmm.settings import MAX_REF_MAF_COL_NAME
from pycmm.settings import WES294_OAF_EARLYONSET_AF_COL_NAME
from pycmm.settings import WES294_OAF_BRCS_AF_COL_NAME
from pycmm.cmmlib.taparser import TAVcfReader
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

        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file,
                                 freq_ratios={KG2014OCT_ALL_COL_NAME:0.2},
                                )
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

        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file,
                                 freq_ratios={KG2014OCT_ALL_COL_NAME:0.2},
                                )
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
                                 project_name=None,
                                 sample_info=None,
                                 anno_cols=DFLT_TEST_MUTREP_COLS,
                                 anno_excl_tags=None,
                                 annotated_vcf_tabix=None,
                                 summary_families_sheet=False,
                                 frequency_ratios=None,
                                 rows_filter_actions=None,
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        if project_name is None:
            project_name = self.test_function
        create_jobs_setup_file(project_name=project_name,
                               project_out_dir=self.working_dir,
                               sample_info=sample_info,
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

        self.init_test(self.current_func_name)
        sample_infos = []
        sample_infos.append("fam1:Al-17")
        sample_infos.append("fam2:Alb-31:Al-23")
        sample_infos.append("fam3:Al-36:Al-47:Al-65")
        sample_infos.append("fam4:Al-73:Al-77:Al-92")
        sample_infos.append("fam5:Br-466")
        jobs_setup_file = self.__create_jobs_setup_file(sample_info=",".join(sample_infos),
                                                        )
        pl = MutRepPipeline(jobs_setup_file=jobs_setup_file)
        in_file = join_path(self.data_dir,
                            'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file, family_infos=pl.families_info)
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

        self.init_test(self.current_func_name)
        sample_infos = []
        sample_infos.append("fam1:Al-17")
        sample_infos.append("fam2:Alb-31:Al-23")
        sample_infos.append("fam3:Al-36:Al-47:Al-65")
        sample_infos.append("fam4:Al-73:Al-77:Al-92")
        sample_infos.append("fam5:Br-466")
        jobs_setup_file = self.__create_jobs_setup_file(sample_info=",".join(sample_infos),
                                                        )
        pl = MutRepPipeline(jobs_setup_file=jobs_setup_file)
        in_file = join_path(self.data_dir,
                            'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file, family_infos=pl.families_info)
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

        self.init_test(self.current_func_name)
        sample_infos = []
        sample_infos.append("fam1:Al-17")
        sample_infos.append("fam2:Alb-31:Al-23")
        sample_infos.append("fam3:Al-36:Al-47:Al-65")
        sample_infos.append("fam4:Al-73:Al-77:Al-92")
        sample_infos.append("fam5:Br-466")
        jobs_setup_file = self.__create_jobs_setup_file(sample_info=",".join(sample_infos),
                                                        )
        pl = MutRepPipeline(jobs_setup_file=jobs_setup_file)
        in_file = join_path(self.data_dir,
                            'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file, family_infos=pl.families_info)
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

        self.init_test(self.current_func_name)
        sample_infos = []
        sample_infos.append("fam1:Al-17")
        sample_infos.append("fam2:Alb-31:Al-23")
        sample_infos.append("fam3:Al-36:Al-47:Al-65")
        sample_infos.append("fam4:Al-73:Al-77:Al-92")
        sample_infos.append("fam5:Br-466")
        jobs_setup_file = self.__create_jobs_setup_file(sample_info=",".join(sample_infos),
                                                        )
        pl = MutRepPipeline(jobs_setup_file=jobs_setup_file)
        in_file = join_path(self.data_dir,
                            'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file, family_infos=pl.families_info)
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

    def test_is_pass_vqsr_1(self):
        """
        test if is_pass_vqsr can differentiate between "PASS", "VQSRTranch",
        and "."
        """

        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                            'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertTrue(vcf_record.is_pass_vqsr(allele_idx=1),
                        "VQSR cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertFalse(vcf_record.is_pass_vqsr(allele_idx=1),
                         "VQSR cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertFalse(vcf_record.is_pass_vqsr(allele_idx=1),
                         "VQSR cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertTrue(vcf_record.is_pass_vqsr(allele_idx=1),
                        "VQSR cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertTrue(vcf_record.is_pass_vqsr(allele_idx=1),
                        "VQSR cannot be correctly determined")

    def test_is_intergenic_1(self):
        """
        test counting intergenic mutations can be identified (all are intergenics)
        """

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
        expr1 = '("P-value_8" != \'\') and (float("P-value_8") < 0.02)'
        expr2 = '("P-value_8" != \'\') and (float("P-value_8") < 0.05)'
        expr3 = '("P-value_8" != \'\') and (float("P-value_8") < 1e-019)'
        expr4 = '("P-value_8" != \'\') and (float("P-value_8") < 1e-020)'
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        # allele idx equal to 1 mean the first alternate allele
        # allele idx equal to 2 mean the second alternate allele
        # so on
        self.assertTrue(vcf_record.vcf_eval(expr1, 1),
                        "cannot perform vcf expression evaluation correctly")
        self.assertTrue(vcf_record.vcf_eval(expr2, 1),
                        "cannot perform vcf expression evaluation correctly")
        self.assertTrue(vcf_record.vcf_eval(expr3, 1),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr4, 1),
                        "cannot perform vcf expression evaluation correctly")
        vcf_record = vcf_reader.next()
        self.assertFalse(vcf_record.vcf_eval(expr1, 1),
                        "cannot perform vcf expression evaluation correctly")
        self.assertTrue(vcf_record.vcf_eval(expr2, 1),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr3, 1),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr4, 1),
                        "cannot perform vcf expression evaluation correctly")
        vcf_record = vcf_reader.next()
        self.assertFalse(vcf_record.vcf_eval(expr1, 1),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr2, 1),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr3, 1),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr4, 1),
                        "cannot perform vcf expression evaluation correctly")

    def test_vcf_eval_2(self):
        """
        test vcf_eval with multi alleleic INFO
        """

        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                            'input.vcf.gz')
        expr1 = '("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" == \'NA\')'
        expr2 = '("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" == \'INF\')'
        expr3 = '("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" == \'\')'
        expr4 = '("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'NA\')'
        expr4 += ' and ("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'INF\')'
        expr4 += ' and ("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'\')'
        expr4 += ' and (float("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '") < 0.5)'
        expr5 = '("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'NA\')'
        expr5 += ' and ("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'INF\')'
        expr5 += ' and ("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'\')'
        expr5 += ' and (float("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '") < 0.8)'
        expr6 = '("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'NA\')'
        expr6 += ' and ("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'INF\')'
        expr6 += ' and ("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'\')'
        expr6 += ' and (float("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '") < 2)'
        expr7 = '("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" == \'NA\')'
        expr8 = '("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" == \'INF\')'
        expr9 = '("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" == \'\')'
        expr10 = '("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'NA\')'
        expr10 += ' and ("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'INF\')'
        expr10 += ' and ("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'\')'
        expr10 += ' and (float("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '") < 0.5)'
        expr11 = '("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'NA\')'
        expr11 += ' and ("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'INF\')'
        expr11 += ' and ("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'\')'
        expr11 += ' and (float("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '") < 0.8)'
        expr12 = '("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'NA\')'
        expr12 += ' and ("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'INF\')'
        expr12 += ' and ("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'\')'
        expr12 += ' and (float("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '") < 2)'
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        vcf_record = vcf_reader.next()
        # allele idx equal to 1 mean the first alternate allele
        # allele idx equal to 2 mean the second alternate allele
        # so on
        self.assertTrue(vcf_record.vcf_eval(expr1, 1),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr2, 1),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr3, 1),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr4, 1),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr5, 1),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr6, 1),
                         "cannot perform vcf expression evaluation correctly")
        vcf_record = vcf_reader.next()
        vcf_record = vcf_reader.next()
        vcf_record = vcf_reader.next()
        self.assertFalse(vcf_record.vcf_eval(expr7, 1),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr8, 1),
                         "cannot perform vcf expression evaluation correctly")
        self.assertTrue(vcf_record.vcf_eval(expr9, 1),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr10, 1),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr11, 1),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr12, 1),
                         "cannot perform vcf expression evaluation correctly")
        vcf_record = vcf_reader.next()
        self.assertFalse(vcf_record.vcf_eval(expr1, 2),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr2, 2),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr3, 2),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr4, 2),
                         "cannot perform vcf expression evaluation correctly")
        self.assertTrue(vcf_record.vcf_eval(expr5, 2),
                        "cannot perform vcf expression evaluation correctly")
        self.assertTrue(vcf_record.vcf_eval(expr6, 2),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr1, 3),
                         "cannot perform vcf expression evaluation correctly")
        self.assertTrue(vcf_record.vcf_eval(expr2, 3),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr3, 3),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr4, 3),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr5, 3),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr6, 3),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr7, 3),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr8, 3),
                         "cannot perform vcf expression evaluation correctly")
        self.assertTrue(vcf_record.vcf_eval(expr9, 3),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr10, 3),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr11, 3),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr12, 3),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr1, 4),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr2, 4),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr3, 4),
                         "cannot perform vcf expression evaluation correctly")
        self.assertTrue(vcf_record.vcf_eval(expr4, 4),
                        "cannot perform vcf expression evaluation correctly")
        self.assertTrue(vcf_record.vcf_eval(expr5, 4),
                        "cannot perform vcf expression evaluation correctly")
        self.assertTrue(vcf_record.vcf_eval(expr6, 4),
                        "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr1, 5),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr2, 5),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr3, 5),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr4, 5),
                         "cannot perform vcf expression evaluation correctly")
        self.assertFalse(vcf_record.vcf_eval(expr5, 5),
                         "cannot perform vcf expression evaluation correctly")
        self.assertTrue(vcf_record.vcf_eval(expr6, 5),
                        "cannot perform vcf expression evaluation correctly")

    def test_pathogenic_count_1(self):
        """
        test if harmful pathogenic variants can be correctly count
        """

        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                            'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.pathogenic_count(allele_idx=1),
                         1,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.pathogenic_count(allele_idx=1),
                         9,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.pathogenic_count(allele_idx=2),
                         0,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        self.assertEqual(vcf_record.pathogenic_count(allele_idx=3),
                         1,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        vcf_record = vcf_reader.next()
        vcf_record = vcf_reader.next()
        vcf_record = vcf_reader.next()
        vcf_record = vcf_reader.next()
        vcf_record = vcf_reader.next()
        vcf_record = vcf_reader.next()
        vcf_record = vcf_reader.next()
        vcf_record = vcf_reader.next()
        vcf_record = vcf_reader.next()
        vcf_record = vcf_reader.next()
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.pathogenic_count(allele_idx=1),
                         1,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        self.assertEqual(vcf_record.pathogenic_count(allele_idx=2),
                         9,
                         "number of harmful pathogenic prediction cannot be correctly determined")

    def test_parse_exac_constraint_1(self):
        """
        test parsing exac constraint for the annotated vcf
        """

        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                            'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        vcf_record = vcf_reader.next()
        self.assertTrue(vcf_record.get_info(EXAC03_CONSTRAINT_EXP_SYN_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        self.assertTrue(vcf_record.get_info(EXAC03_CONSTRAINT_N_SYN_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        self.assertTrue(vcf_record.get_info(EXAC03_CONSTRAINT_SYN_Z_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        self.assertTrue(vcf_record.get_info(EXAC03_CONSTRAINT_EXP_MIS_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        self.assertTrue(vcf_record.get_info(EXAC03_CONSTRAINT_N_MIS_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        self.assertTrue(vcf_record.get_info(EXAC03_CONSTRAINT_MIS_Z_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        self.assertTrue(vcf_record.get_info(EXAC03_CONSTRAINT_EXP_LOF_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        self.assertTrue(vcf_record.get_info(EXAC03_CONSTRAINT_N_LOF_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        self.assertTrue(vcf_record.get_info(EXAC03_CONSTRAINT_PLI_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.get_info(EXAC03_CONSTRAINT_EXP_SYN_COL_NAME),
                         '45.7717977506',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(vcf_record.get_info(EXAC03_CONSTRAINT_N_SYN_COL_NAME),
                         '40',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(vcf_record.get_info(EXAC03_CONSTRAINT_SYN_Z_COL_NAME),
                         '0.528885612026385',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(vcf_record.get_info(EXAC03_CONSTRAINT_EXP_MIS_COL_NAME),
                         '115.696060891',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(vcf_record.get_info(EXAC03_CONSTRAINT_N_MIS_COL_NAME),
                         '34',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(vcf_record.get_info(EXAC03_CONSTRAINT_MIS_Z_COL_NAME),
                         '3.71499650338098',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(vcf_record.get_info(EXAC03_CONSTRAINT_EXP_LOF_COL_NAME),
                         '15.533928734',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(vcf_record.get_info(EXAC03_CONSTRAINT_N_LOF_COL_NAME),
                         '1',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(vcf_record.get_info(EXAC03_CONSTRAINT_PLI_COL_NAME),
                         '0.975506865848027',
                         "values of ExAC constraint cannot be correctly determined")

    def test_parse_exac_constraint_2(self):
        """
        test parsing exac constraint for vcf without the annotation
        """

        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                            'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertTrue(vcf_record.get_info(EXAC03_CONSTRAINT_EXP_SYN_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")

    def test_parse_exac_constraint_4(self):
        """
        test parsing exac constraint for vcf without the annotation
        """

        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                            'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.get_info(EXAC03_CONSTRAINT_N_LOF_COL_NAME, 1),
                         '',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(vcf_record.get_info(EXAC03_CONSTRAINT_PLI_COL_NAME, 2),
                         '',
                         "values of ExAC constraint cannot be correctly determined")

    def test_max_ref_maf_1(self):
        """
        test finding maximum reference minor allele frequency
        """

        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                            'input.vcf.gz')
        vcf_reader = TAVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.get_info(MAX_REF_MAF_COL_NAME, 1),
                         0.0008,
                         "values of MAX_REF_MAF constraint cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.get_info(MAX_REF_MAF_COL_NAME, 1),
                         0.095,
                         "values of MAX_REF_MAF constraint cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.get_info(MAX_REF_MAF_COL_NAME, 1),
                         0.00119808,
                         "values of MAX_REF_MAF constraint cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.get_info(MAX_REF_MAF_COL_NAME, 1),
                         0.0035,
                         "values of MAX_REF_MAF constraint cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.get_info(MAX_REF_MAF_COL_NAME, 1),
                         0.0125,
                         "values of MAX_REF_MAF constraint cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.get_info(MAX_REF_MAF_COL_NAME, 1),
                         0.0377,
                         "values of MAX_REF_MAF constraint cannot be correctly determined")

    def tearDown(self):
        self.remove_working_dir()
