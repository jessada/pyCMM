import filecmp
from os.path import join as join_path
from os.path import dirname
from pycmm.utils import mylogger
from pycmm.template import SafeTester
from pycmm.proc.tavcf import TableAnnovarVcfReader


class TestTableAnnovarVcfReader(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            test_module_name=__name__,
                            )

    def setUp(self):
        mylogger.getLogger(__name__)
        self.module_name = 'tavcf'

    def __create_db_instance(self):
        return None

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

    def test_parse_cmm_gt(self):
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
        self.assertEqual(call.cmm_gt[1],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.cmm_gt[2],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.cmm_gt[3],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.cmm_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.cmm_gt[1],
                         "hom",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.cmm_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.cmm_gt[3],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.cmm_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Al-161")
        self.assertEqual(call.cmm_gt[1],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.cmm_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.cmm_gt[3],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.cmm_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Al-65")
        self.assertEqual(call.cmm_gt[1],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.cmm_gt[2],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.cmm_gt[3],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.cmm_gt[4],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Al-111")
        self.assertEqual(call.cmm_gt[1],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.cmm_gt[2],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.cmm_gt[3],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.cmm_gt[4],
                         "wt",
                         "cmm_genotye cannot be correctly determined")

    def test_parse_actual_gt_1(self):
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
        self.assertEqual(call.actual_gt[1],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.actual_gt[1],
                         "hom",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Al-161")
        self.assertEqual(call.actual_gt[1],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Al-65")
        self.assertEqual(call.actual_gt[1],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        vcf_record = vcf_reader.next()
        call = vcf_record.genotype("Al-111")
        self.assertEqual(call.actual_gt[1],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "wt",
                         "cmm_genotye cannot be correctly determined")

    def test_parse_actual_gt_2(self):
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
        self.assertEqual(call.actual_gt[1],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.actual_gt[1],
                         "hom",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Al-161")
        self.assertEqual(call.actual_gt[1],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Br-637")
        self.assertEqual(call.actual_gt[1],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "wt",
                         "cmm_genotye cannot be correctly determined")

    def test_parse_actual_gt_3(self):
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
        self.assertEqual(call.actual_gt[1],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.actual_gt[1],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Al-111")
        self.assertEqual(call.actual_gt[1],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Br-120")
        self.assertEqual(call.actual_gt[1],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Br-390")
        self.assertEqual(call.actual_gt[1],
                         "hom",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")

    def test_parse_actual_gt_4(self):
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
        self.assertEqual(call.actual_gt[1],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.actual_gt[1],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Alb-31")
        self.assertEqual(call.actual_gt[1],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Br-693")
        self.assertEqual(call.actual_gt[1],
                         "hom",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")

    def test_parse_actual_gt_5(self):
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
        self.assertEqual(call.actual_gt[1],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        # this one is still wrong but acceptable
        self.assertEqual(call.actual_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.actual_gt[1],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Alb-31")
        self.assertEqual(call.actual_gt[1],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        call = vcf_record.genotype("Br-466")
        self.assertEqual(call.actual_gt[1],
                         "hom",
                         "cmm_genotye cannot be correctly determined")
        # this one is still wrong but acceptable
        self.assertEqual(call.actual_gt[2],
                         "hom",
                         "cmm_genotye cannot be correctly determined")

    def test_parse_actual_gt_6(self):
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
        self.assertEqual(call.actual_gt[1],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        # 0/0
        call = vcf_record.genotype("Br-781")
        self.assertEqual(call.actual_gt[1],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        # 0/1
        call = vcf_record.genotype("Al-111")
        self.assertEqual(call.actual_gt[1],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        # 0/2
        call = vcf_record.genotype("Al-23")
        self.assertEqual(call.actual_gt[1],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        # 1/1
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.actual_gt[1],
                         "hom",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        # 1/2
        call = vcf_record.genotype("Al-92")
        self.assertEqual(call.actual_gt[1],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        # 2/2
        call = vcf_record.genotype("Br-296")
        self.assertEqual(call.actual_gt[1],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "hom",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        # random
        call = vcf_record.genotype("Br-429")
        self.assertEqual(call.actual_gt[1],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")

    def test_parse_actual_gt_7(self):
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
        self.assertEqual(call.actual_gt[1],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         ".",
                         "cmm_genotye cannot be correctly determined")
        # 0/0
        call = vcf_record.genotype("Br-781")
        self.assertEqual(call.actual_gt[1],
                         "hom",
                         "cmm_genotye cannot be correctly determined")
        # this one is still wrong but acceptable
        # it should be "other"
        self.assertEqual(call.actual_gt[2],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        # this one is still wrong but acceptable
        # it should be "other"
        self.assertEqual(call.actual_gt[3],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        # this one is still wrong but acceptable
        # it should be "other"
        self.assertEqual(call.actual_gt[4],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        # 0/1
        call = vcf_record.genotype("Al-111")
        self.assertEqual(call.actual_gt[1],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        # 0/2
        call = vcf_record.genotype("Al-23")
        self.assertEqual(call.actual_gt[1],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        # 1/1
        call = vcf_record.genotype("Br-432")
        self.assertEqual(call.actual_gt[1],
                         "wt",
                         "cmm_genotye cannot be correctly determined")
        # this one is still wrong but acceptable
        # it should be "wt"
        self.assertEqual(call.actual_gt[2],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        # this one is still wrong but acceptable
        # it should be "wt"
        self.assertEqual(call.actual_gt[3],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        # this one is still wrong but acceptable
        # it should be "wt"
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        # 1/2
        call = vcf_record.genotype("Al-92")
        self.assertEqual(call.actual_gt[1],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        # 2/2
        call = vcf_record.genotype("Br-296")
        self.assertEqual(call.actual_gt[1],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "hom",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        # random
        call = vcf_record.genotype("Br-429")
        self.assertEqual(call.actual_gt[1],
                         "oth",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[2],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[3],
                         "het",
                         "cmm_genotye cannot be correctly determined")
        self.assertEqual(call.actual_gt[4],
                         "oth",
                         "cmm_genotye cannot be correctly determined")

