import filecmp
import unittest
from os.path import join as join_path
from os.path import dirname
from os.path import basename
from pycmm.settings import FULL_SYSTEM_TEST
from pycmm.settings import DFLT_ANV_DB_DIR
from pycmm.settings import DFLT_ANV_DB_NAMES
from pycmm.settings import DFLT_ANV_DB_OPS
from pycmm.template import SafeTester
from pycmm.utils import count_lines
from pycmm.flow.mutrep import MutRepPipeline
from pycmm.cmmlib.annovarlib import Annovar
from pycmm.cmmlib.annovarlib import PredictionTranslator
from pycmm.cmmlib.annovarlib import SIFT_PRED_COL
from pycmm.cmmlib.annovarlib import POLYPHEN2_HDIV_PRED_COL
from pycmm.cmmlib.annovarlib import POLYPHEN2_HVAR_PRED_COL
from pycmm.cmmlib.annovarlib import LRT_PRED_COL
from pycmm.cmmlib.annovarlib import MUTATIONTASTER_PRED_COL
from pycmm.cmmlib.annovarlib import MUTATIONASSESSOR_PRED_COL
from pycmm.cmmlib.annovarlib import FATHMM_PRED_COL
from pycmm.cmmlib.annovarlib import RADIALSVM_PRED_COL
from pycmm.cmmlib.annovarlib import LR_PRED_COL

DFLT_ANV_TEST_DB_DIR = DFLT_ANV_DB_DIR
DFLT_ANV_TEST_DB_NAMES = "refGene"
DFLT_ANV_TEST_DB_OPS = "g"
DFLT_ANV_TEST_DB_NAMES += ",cytoBand"
DFLT_ANV_TEST_DB_OPS += ",r"
DFLT_ANV_TEST_DB_NAMES += ",genomicSuperDups"
DFLT_ANV_TEST_DB_OPS += ",r"


class TestAnnovar(SafeTester):

    def __init__(self, methodName):
        super(TestAnnovar, self).__init__(methodName=methodName,
                                          test_module_name=__name__,
                                          )

    def setUp(self):
        pass

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_table_annovar_offline_1(self):
        """ test offline version (w/o slurm) of table_annovar """

        self.init_test(self.current_func_name)
        dataset_name = self.current_func_name
        input_file = join_path(self.data_dir,
                               "input.vcf.gz")
        data_out_folder = self.working_dir
        av = Annovar(dataset_name=dataset_name,
                        input_file=input_file,
                        db_dir=DFLT_ANV_TEST_DB_DIR,
                        buildver="hg19",
                        protocols=DFLT_ANV_TEST_DB_NAMES,
                        operations=DFLT_ANV_TEST_DB_OPS,
                        nastring=".",
                        data_out_folder=data_out_folder,
                        )
        av.run_table_annovar()
        self.assertEqual(count_lines(av.out_annotated_vcf),
                         151,
                         "table annovar result is incorrect ")
        exp_annotated_vcf = dataset_name + "_annotated.vcf"
        self.assertEqual(basename(av.out_annotated_vcf),
                         exp_annotated_vcf,
                         "invalid expected output annatated vcf file name")

    def tearDown(self):
        self.remove_working_dir()

    @unittest.skipUnless(FULL_SYSTEM_TEST, "taking too long time to test")
    def test_table_annovar_offline_2(self):
        """ test offline version (w/o slurm) of table_annovar together with custom cal_stat """

        self.init_test(self.current_func_name)
        dataset_name = self.current_func_name
        input_file = join_path(self.data_dir,
                               "input.vcf.gz")
        annovar_db_names = DFLT_ANV_TEST_DB_NAMES
        annovar_db_names += ",test_pyCMM"
        annovar_db_ops = DFLT_ANV_TEST_DB_OPS
        annovar_db_ops += ",f"
        data_out_folder = self.working_dir
        av = Annovar(dataset_name=dataset_name,
                        input_file=input_file,
                        db_dir=DFLT_ANV_TEST_DB_DIR,
                        buildver="hg19",
                        protocols=annovar_db_names,
                        operations=annovar_db_ops,
                        nastring=".",
                        data_out_folder=data_out_folder,
                        )
        av.run_table_annovar()
        self.assertEqual(count_lines(av.out_annotated_vcf),
                         158,
                         "table annoar result is incorrect ")

class TestPredictionTranslator(SafeTester):

    def __init__(self, methodName):
        super(TestPredictionTranslator, self).__init__(methodName=methodName,
                                                       test_module_name=__name__,
                                                       )

    def setUp(self):
        pass

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
                         False,
                         "PredictionTranslator cannot tell harmfulness of "+POLYPHEN2_HDIV_PRED_COL)
        self.assertEqual(pred.get_prediction_info(POLYPHEN2_HDIV_PRED_COL, "B").harmful,
                         False,
                         "PredictionTranslator cannot tell harmfulness of "+POLYPHEN2_HDIV_PRED_COL)
        self.assertEqual(pred.get_prediction_info(POLYPHEN2_HVAR_PRED_COL, "D").harmful,
                         True,
                         "PredictionTranslator cannot tell harmfulness of "+POLYPHEN2_HVAR_PRED_COL)
        self.assertEqual(pred.get_prediction_info(POLYPHEN2_HVAR_PRED_COL, "P").harmful,
                         False,
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
                         True,
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

    def test_is_harmful_2(self):
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
                         True,
                         "PredictionTranslator cannot work correctly if there is no prediction code")

    def tearDown(self):
        self.remove_working_dir()
