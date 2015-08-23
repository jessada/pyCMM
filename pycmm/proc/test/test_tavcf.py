import filecmp
from os.path import join as join_path
from os.path import dirname
from pycmm.utils import mylogger
from pycmm.template import SafeTester
from pycmm.proc.tavcf import TableAnnovarVcfReader
from pycmm.flow.mutrep import MutRepPipeline


class TestTableAnnovarVcfReader(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            test_module_name=__name__,
                            )

    def setUp(self):
        mylogger.getLogger(__name__)

    def test_annovar_info_1(self):
        """ to check if all the variables annotated by ANNOVAR are correctly listed"""

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
        self.assertEqual(len(vcf_reader.annovar_infos),
                         16,
                         "TableAnnovarVcfReader cannot locate ANNOVAR info")

    def test_annovar_info_2(self):
        """ to check if the property work correctly if no info annotated by ANNOVAR """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
        self.assertEqual(len(vcf_reader.annovar_infos),
                         0,
                         "TableAnnovarVcfReader cannot locate ANNOVAR info")

    def test_parse_info_1(self):
        """
        to ensure that comma-separated string annotated by ANNOVAR
        shall be parsed as string
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.INFO['Gene.refGene'],
                         "ZNF264,AURKC",
                         "TableAnnovarVcfReader cannot correctly parse comma-separated info annotated by ANNOVAR")

    def test_parse_info_2(self):
        """ all the hex character shall be converted back to normal character """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
        vcf_reader.next()
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.INFO['GeneDetail.refGene'],
                         "dist=6394;dist=1769",
                         "TableAnnovarVcfReader cannot correctly parse hex-encoded info annotated by ANNOVAR")

    def test_parse_info_3(self):
        """
        if there are more than one alternate alleles,
        shall be parsed as string
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.INFO['AXEQ_BR_AF'],
                         "0.1170",
                         "TableAnnovarVcfReader cannot identify value for info entries annotated by ANNOVAR")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.INFO['AXEQ_BR_AF'],
                         ['0.3068', '0.1818', '0.2614', '0.0341'],
                         "TableAnnovarVcfReader cannot identify value for info entries annotated by ANNOVAR")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.INFO['AXEQ_BR_AF'],
                         ['0.2692', '0.2308'],
                         "TableAnnovarVcfReader cannot identify value for info entries annotated by ANNOVAR")

class TestCmmVcfCall(SafeTester):

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
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
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
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
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
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
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
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
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
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
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
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
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
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
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
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
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
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
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
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
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
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
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
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
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

class TestCmmVcfRecord(SafeTester):

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
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(2, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(3, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(4, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(2, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")

    def test_shared_2(self):
        """ test shared mutation between two samples can be identified """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        samples = ["Alb-31", "Br-466"]
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(2, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(3, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(4, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(2, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")

    def test_shared_3(self):
        """ test shared mutation between three samples can be identified """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        samples = ["Br-432", "Al-161", "Br-504"]
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         True,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(2, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(3, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(4, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(vcf_record.is_shared(2, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
                         False,
                         "shared mutation cannot be correctly determined")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.is_shared(1, samples),
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
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
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
    def test_is_rare_1(self):
        """
        test rare mutation can be identified if there is only one criteria
          - 0.2
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        freq_ratios = {'1000g2014oct_all': 0.2}
        vcf_reader = TableAnnovarVcfReader(filename=in_file)
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
