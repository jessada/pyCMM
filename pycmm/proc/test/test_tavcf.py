import filecmp
from os.path import join as join_path
from os.path import dirname
from pycmm.utils import mylogger
from pycmm.template import SafeTester
from pycmm.proc.tavcf import TableAnnovarVCFReader


class TestTableAnnovarVCFReader(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            test_module_name=__name__,
                            )

    def setUp(self):
        self.module_name = 'tavcf'

    def __create_db_instance(self):
        return None

    def test_annovar_info_1(self):
        """ to check if all the variables annotated by ANNOVAR are correctly listed"""

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TableAnnovarVCFReader(filename=in_file)
        self.assertEqual(len(vcf_reader.annovar_infos),
                         16,
                         "TableAnnovarVCFReader cannot locate ANNOVAR info")

    def test_annovar_info_2(self):
        """ to check if the property work correctly if no info annotated by ANNOVAR """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TableAnnovarVCFReader(filename=in_file)
        self.assertEqual(len(vcf_reader.annovar_infos),
                         0,
                         "TableAnnovarVCFReader cannot locate ANNOVAR info")

    def test_parse_info_1(self):
        """
        to ensure that comma-separated string annotated by ANNOVAR
        shall be parsed as string
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TableAnnovarVCFReader(filename=in_file)
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.INFO['Gene.refGene'],
                         "ZNF264,AURKC",
                         "TableAnnovarVCFReader cannot correctly parse comma-separated info annotated by ANNOVAR")

    def test_parse_info_2(self):
        """ all the hex character shall be converted back to normal character """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TableAnnovarVCFReader(filename=in_file)
        vcf_reader.next()
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.INFO['GeneDetail.refGene'],
                         "dist=6394;dist=1769",
                         "TableAnnovarVCFReader cannot correctly parse hex-encoded info annotated by ANNOVAR")

    def test_parse_info_3(self):
        """
        if there are more than one alternate alleles,
        shall be parsed as string
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        in_file = join_path(self.data_dir,
                               'input.vcf.gz')
        vcf_reader = TableAnnovarVCFReader(filename=in_file)
        vcf_reader.next()
        vcf_reader.next()
        vcf_reader.next()
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.INFO['AXEQ_BR_AF'],
                         "0.1170",
                         "TableAnnovarVCFReader cannot identify value for info entries annotated by ANNOVAR")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.INFO['AXEQ_BR_AF'],
                         ['0.3068', '0.1818', '0.2614', '0.0341'],
                         "TableAnnovarVCFReader cannot identify value for info entries annotated by ANNOVAR")
        vcf_record = vcf_reader.next()
        self.assertEqual(vcf_record.INFO['AXEQ_BR_AF'],
                         ['0.2692', '0.2308'],
                         "TableAnnovarVCFReader cannot identify value for info entries annotated by ANNOVAR")
