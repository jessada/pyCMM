import filecmp
from os.path import join as join_path
from os.path import dirname
from pycmm.utils import mylogger
from pycmm.template import SafeTester
#from pycmm.proc.tavcf import TableAnnovarVcfReader
from pycmm.flow.mutrep import MutRepPipeline
from pycmm.proc.annovarlib import PredictionTranslator
from pycmm.proc.annovarlib import SIFT_PRED_COL
from pycmm.proc.annovarlib import POLYPHEN2_HDIV_PRED_COL
from pycmm.proc.annovarlib import POLYPHEN2_HVAR_PRED_COL
from pycmm.proc.annovarlib import LRT_PRED_COL
from pycmm.proc.annovarlib import MUTATIONTASTER_PRED_COL
from pycmm.proc.annovarlib import MUTATIONASSESSOR_PRED_COL
from pycmm.proc.annovarlib import FATHMM_PRED_COL
from pycmm.proc.annovarlib import RADIALSVM_PRED_COL
from pycmm.proc.annovarlib import LR_PRED_COL


class TestPredictionTranslator(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            test_module_name=__name__,
                            )

    def setUp(self):
        mylogger.getLogger(__name__)

    def test_description_1(self):
        """ to test if normal description of effect predictor can be shown correctly """

        self.init_test(self.current_func_name)
        pred = PredictionTranslator()
        self.assertEqual(pred.get_prediction_info(SIFT_PRED_COL, "D").description,
                         "Deleterious",
                         "PredictionTranslator cannot describe "+SIFT_PRED_COL)
        self.assertEqual(pred.get_prediction_info(SIFT_PRED_COL, "T").description,
                         "Tolerated",
                         "PredictionTranslator cannot describe "+SIFT_PRED_COL)
        self.assertEqual(pred.get_prediction_info(POLYPHEN2_HDIV_PRED_COL, "D").description,
                         "Probably damaging",
                         "PredictionTranslator cannot describe "+POLYPHEN2_HDIV_PRED_COL)
        self.assertEqual(pred.get_prediction_info(POLYPHEN2_HDIV_PRED_COL, "P").description,
                         "Possibly damaging",
                         "PredictionTranslator cannot describe "+POLYPHEN2_HDIV_PRED_COL)
        self.assertEqual(pred.get_prediction_info(POLYPHEN2_HDIV_PRED_COL, "B").description,
                         "Benign",
                         "PredictionTranslator cannot describe "+POLYPHEN2_HDIV_PRED_COL)
        self.assertEqual(pred.get_prediction_info(POLYPHEN2_HVAR_PRED_COL, "D").description,
                         "Probably damaging",
                         "PredictionTranslator cannot describe "+POLYPHEN2_HVAR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(POLYPHEN2_HVAR_PRED_COL, "P").description,
                         "Possibly damaging",
                         "PredictionTranslator cannot describe "+POLYPHEN2_HVAR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(POLYPHEN2_HVAR_PRED_COL, "B").description,
                         "Benign",
                         "PredictionTranslator cannot describe "+POLYPHEN2_HVAR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(LRT_PRED_COL, "D").description,
                         "Deleterious",
                         "PredictionTranslator cannot describe "+LRT_PRED_COL)
        self.assertEqual(pred.get_prediction_info(LRT_PRED_COL, "N").description,
                         "Neutral",
                         "PredictionTranslator cannot describe "+LRT_PRED_COL)
        self.assertEqual(pred.get_prediction_info(LRT_PRED_COL, "U").description,
                         "Unknown",
                         "PredictionTranslator cannot describe "+LRT_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONTASTER_PRED_COL, "A").description,
                         "Disease causing automatic",
                         "PredictionTranslator cannot describe "+MUTATIONTASTER_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONTASTER_PRED_COL, "D").description,
                         "Disease causing",
                         "PredictionTranslator cannot describe "+MUTATIONTASTER_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONTASTER_PRED_COL, "N").description,
                         "Polymorphism",
                         "PredictionTranslator cannot describe "+MUTATIONTASTER_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONTASTER_PRED_COL, "P").description,
                         "Polymorphism automatic",
                         "PredictionTranslator cannot describe "+MUTATIONTASTER_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONASSESSOR_PRED_COL, "H").description,
                         "High",
                         "PredictionTranslator cannot describe "+MUTATIONASSESSOR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONASSESSOR_PRED_COL, "M").description,
                         "Medium",
                         "PredictionTranslator cannot describe "+MUTATIONASSESSOR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONASSESSOR_PRED_COL, "L").description,
                         "Low",
                         "PredictionTranslator cannot describe "+MUTATIONASSESSOR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONASSESSOR_PRED_COL, "N").description,
                         "Neutral",
                         "PredictionTranslator cannot describe "+MUTATIONASSESSOR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONASSESSOR_PRED_COL, "H/M").description,
                         "Functional",
                         "PredictionTranslator cannot describe "+MUTATIONASSESSOR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONASSESSOR_PRED_COL, "L/N").description,
                         "Non-functional",
                         "PredictionTranslator cannot describe "+MUTATIONASSESSOR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(FATHMM_PRED_COL, "D").description,
                         "Deleterious",
                         "PredictionTranslator cannot describe "+FATHMM_PRED_COL)
        self.assertEqual(pred.get_prediction_info(FATHMM_PRED_COL, "T").description,
                         "Tolerated",
                         "PredictionTranslator cannot describe "+FATHMM_PRED_COL)
        self.assertEqual(pred.get_prediction_info(RADIALSVM_PRED_COL, "D").description,
                         "Deleterious",
                         "PredictionTranslator cannot describe "+RADIALSVM_PRED_COL)
        self.assertEqual(pred.get_prediction_info(RADIALSVM_PRED_COL, "T").description,
                         "Tolerated",
                         "PredictionTranslator cannot describe "+RADIALSVM_PRED_COL)
        self.assertEqual(pred.get_prediction_info(LR_PRED_COL, "D").description,
                         "Deleterious",
                         "PredictionTranslator cannot describe "+LR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(LR_PRED_COL, "T").description,
                         "Tolerated",
                         "PredictionTranslator cannot describe "+LR_PRED_COL)

    def test_description_2(self):
        """ to test specials case of the description of effect predictor code """
        self.init_test(self.current_func_name)
        pred = PredictionTranslator()
        self.assertEqual(pred.get_prediction_info('ABC predictor', "K").description,
                         "NA(K)",
                         "PredictionTranslator cannot work correctly if the predictor is unknown")
        self.assertEqual(pred.get_prediction_info(LR_PRED_COL, "J").description,
                         "NA(J)",
                         "PredictionTranslator cannot work correctly if the prediction code is unknown")
        self.assertEqual(pred.get_prediction_info(LR_PRED_COL, ".").description,
                         ".",
                         "PredictionTranslator cannot work correctly if there is no prediction code")

    def test_is_harmful_1(self):
        """ to if the predictor can normally tell the harmfulness """

        self.init_test(self.current_func_name)
        pred = PredictionTranslator()
        self.assertEqual(pred.get_prediction_info(SIFT_PRED_COL, "D").harmful,
                         True,
                         "PredictionTranslator cannot tell harmfulness of "+SIFT_PRED_COL)
        self.assertEqual(pred.get_prediction_info(SIFT_PRED_COL, "T").harmful,
                         False,
                         "PredictionTranslator cannot tell harmfulness of "+SIFT_PRED_COL)
        self.assertEqual(pred.get_prediction_info(POLYPHEN2_HDIV_PRED_COL, "D").harmful,
                         True,
                         "PredictionTranslator cannot tell harmfulness of "+POLYPHEN2_HDIV_PRED_COL)
        self.assertEqual(pred.get_prediction_info(POLYPHEN2_HDIV_PRED_COL, "P").harmful,
                         True,
                         "PredictionTranslator cannot tell harmfulness of "+POLYPHEN2_HDIV_PRED_COL)
        self.assertEqual(pred.get_prediction_info(POLYPHEN2_HDIV_PRED_COL, "B").harmful,
                         False,
                         "PredictionTranslator cannot tell harmfulness of "+POLYPHEN2_HDIV_PRED_COL)
        self.assertEqual(pred.get_prediction_info(POLYPHEN2_HVAR_PRED_COL, "D").harmful,
                         True,
                         "PredictionTranslator cannot tell harmfulness of "+POLYPHEN2_HVAR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(POLYPHEN2_HVAR_PRED_COL, "P").harmful,
                         True,
                         "PredictionTranslator cannot tell harmfulness of "+POLYPHEN2_HVAR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(POLYPHEN2_HVAR_PRED_COL, "B").harmful,
                         False,
                         "PredictionTranslator cannot tell harmfulness of "+POLYPHEN2_HVAR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(LRT_PRED_COL, "D").harmful,
                         True,
                         "PredictionTranslator cannot tell harmfulness of "+LRT_PRED_COL)
        self.assertEqual(pred.get_prediction_info(LRT_PRED_COL, "N").harmful,
                         False,
                         "PredictionTranslator cannot tell harmfulness of "+LRT_PRED_COL)
        self.assertEqual(pred.get_prediction_info(LRT_PRED_COL, "U").harmful,
                         False,
                         "PredictionTranslator cannot tell harmfulness of "+LRT_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONTASTER_PRED_COL, "A").harmful,
                         True,
                         "PredictionTranslator cannot tell harmfulness of "+MUTATIONTASTER_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONTASTER_PRED_COL, "D").harmful,
                         True,
                         "PredictionTranslator cannot tell harmfulness of "+MUTATIONTASTER_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONTASTER_PRED_COL, "N").harmful,
                         False,
                         "PredictionTranslator cannot tell harmfulness of "+MUTATIONTASTER_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONTASTER_PRED_COL, "P").harmful,
                         False,
                         "PredictionTranslator cannot tell harmfulness of "+MUTATIONTASTER_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONASSESSOR_PRED_COL, "H").harmful,
                         True,
                         "PredictionTranslator cannot tell harmfulness of "+MUTATIONASSESSOR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONASSESSOR_PRED_COL, "M").harmful,
                         True,
                         "PredictionTranslator cannot tell harmfulness of "+MUTATIONASSESSOR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONASSESSOR_PRED_COL, "L").harmful,
                         False,
                         "PredictionTranslator cannot tell harmfulness of "+MUTATIONASSESSOR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONASSESSOR_PRED_COL, "N").harmful,
                         False,
                         "PredictionTranslator cannot tell harmfulness of "+MUTATIONASSESSOR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONASSESSOR_PRED_COL, "H/M").harmful,
                         True,
                         "PredictionTranslator cannot tell harmfulness of "+MUTATIONASSESSOR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(MUTATIONASSESSOR_PRED_COL, "L/N").harmful,
                         False,
                         "PredictionTranslator cannot tell harmfulness of "+MUTATIONASSESSOR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(FATHMM_PRED_COL, "D").harmful,
                         True,
                         "PredictionTranslator cannot tell harmfulness of "+FATHMM_PRED_COL)
        self.assertEqual(pred.get_prediction_info(FATHMM_PRED_COL, "T").harmful,
                         False,
                         "PredictionTranslator cannot tell harmfulness of "+FATHMM_PRED_COL)
        self.assertEqual(pred.get_prediction_info(RADIALSVM_PRED_COL, "D").harmful,
                         True,
                         "PredictionTranslator cannot tell harmfulness of "+RADIALSVM_PRED_COL)
        self.assertEqual(pred.get_prediction_info(RADIALSVM_PRED_COL, "T").harmful,
                         False,
                         "PredictionTranslator cannot tell harmfulness of "+RADIALSVM_PRED_COL)
        self.assertEqual(pred.get_prediction_info(LR_PRED_COL, "D").harmful,
                         True,
                         "PredictionTranslator cannot tell harmfulness of "+LR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(LR_PRED_COL, "T").harmful,
                         False,
                         "PredictionTranslator cannot tell harmfulness of "+LR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(LR_PRED_COL, "K").harmful,
                         False,
                         "PredictionTranslator cannot tell harmfulness of "+LR_PRED_COL)

    def test_is_harmful_1(self):
        """ to if the predictor can tell the harmfulness in special scenarios"""
        self.init_test(self.current_func_name)
        pred = PredictionTranslator()
        self.assertEqual(pred.get_prediction_info('ABC predictor', "K").harmful,
                         False,
                         "PredictionTranslator cannot work correctly if the predictor is unknown")
        self.assertEqual(pred.get_prediction_info(LR_PRED_COL, "J").harmful,
                         False,
                         "PredictionTranslator cannot work correctly if the prediction code is unknown")
        self.assertEqual(pred.get_prediction_info(LR_PRED_COL, ".").harmful,
                         False,
                         "PredictionTranslator cannot work correctly if there is no prediction code")
#    def test_annovar_info_2(self):
#        """ to check if the property work correctly if no info annotated by ANNOVAR """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        in_file = join_path(self.data_dir,
#                               'input.vcf.gz')
#        vcf_reader = TableAnnovarVcfReader(filename=in_file)
#        self.assertEqual(len(vcf_reader.annovar_infos),
#                         0,
#                         "TableAnnovarVcfReader cannot locate ANNOVAR info")
#
#    def test_parse_info_1(self):
#        """
#        to ensure that comma-separated string annotated by ANNOVAR
#        shall be parsed as string
#        """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        in_file = join_path(self.data_dir,
#                               'input.vcf.gz')
#        vcf_reader = TableAnnovarVcfReader(filename=in_file)
#        vcf_reader.next()
#        vcf_reader.next()
#        vcf_reader.next()
#        vcf_record = vcf_reader.next()
#        self.assertEqual(vcf_record.INFO['Gene.refGene'],
#                         "ZNF264,AURKC",
#                         "TableAnnovarVcfReader cannot correctly parse comma-separated info annotated by ANNOVAR")
#
#    def test_parse_info_2(self):
#        """ all the hex character shall be converted back to normal character """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        in_file = join_path(self.data_dir,
#                               'input.vcf.gz')
#        vcf_reader = TableAnnovarVcfReader(filename=in_file)
#        vcf_reader.next()
#        vcf_reader.next()
#        vcf_record = vcf_reader.next()
#        self.assertEqual(vcf_record.INFO['GeneDetail.refGene'],
#                         "dist=6394;dist=1769",
#                         "TableAnnovarVcfReader cannot correctly parse hex-encoded info annotated by ANNOVAR")
#
#    def test_parse_info_3(self):
#        """
#        if there are more than one alternate alleles,
#        shall be parsed as string
#        """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        in_file = join_path(self.data_dir,
#                               'input.vcf.gz')
#        vcf_reader = TableAnnovarVcfReader(filename=in_file)
#        vcf_reader.next()
#        vcf_reader.next()
#        vcf_reader.next()
#        vcf_record = vcf_reader.next()
#        self.assertEqual(vcf_record.INFO['AXEQ_BR_AF'],
#                         "0.1170",
#                         "TableAnnovarVcfReader cannot identify value for info entries annotated by ANNOVAR")
#        vcf_record = vcf_reader.next()
#        self.assertEqual(vcf_record.INFO['AXEQ_BR_AF'],
#                         ['0.3068', '0.1818', '0.2614', '0.0341'],
#                         "TableAnnovarVcfReader cannot identify value for info entries annotated by ANNOVAR")
#        vcf_record = vcf_reader.next()
#        self.assertEqual(vcf_record.INFO['AXEQ_BR_AF'],
#                         ['0.2692', '0.2308'],
#                         "TableAnnovarVcfReader cannot identify value for info entries annotated by ANNOVAR")
