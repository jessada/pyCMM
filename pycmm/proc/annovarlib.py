import sys
from pycmm.template import pyCMMBase
from pycmm.settings import LJB_SIFT_PREDICTION_COL_NAME as SIFT_PRED_COL
from pycmm.settings import LJB_POLYPHEN2_HDIV_PREDICTION_COL_NAME as POLYPHEN2_HDIV_PRED_COL
from pycmm.settings import LJB_POLYPHEN2_HVAR_PREDICTION_COL_NAME as POLYPHEN2_HVAR_PRED_COL
from pycmm.settings import LJB_LRT_PREDICTION_COL_NAME as LRT_PRED_COL
from pycmm.settings import LJB_MUTATIONTASTER_PREDICTION_COL_NAME as MUTATIONTASTER_PRED_COL
from pycmm.settings import LJB_MUTATIONASSESSOR_PREDICTION_COL_NAME as MUTATIONASSESSOR_PRED_COL
from pycmm.settings import LJB_FATHMM_PREDICTION_COL_NAME as FATHMM_PRED_COL
from pycmm.settings import LJB_RADIALSVM_PREDICTION_COL_NAME as RADIALSVM_PRED_COL
from pycmm.settings import LJB_LR_PREDICTION_COL_NAME as LR_PRED_COL

ANNOVAR_PARAMS_INPUT_FILE_KEY = "input_file"
ANNOVAR_PARAMS_DB_FOLDER_KEY = "db_folder"
ANNOVAR_PARAMS_BUILDVER_KEY = "buildver"
ANNOVAR_PARAMS_OUT_PREFIX_KEY = "out_prefix"
ANNOVAR_PARAMS_DB_NAMES_KEY = "db_names"
ANNOVAR_PARAMS_DB_OPS_KEY = "db_ops"
ANNOVAR_PARAMS_NASTRING_KEY = "nastring"

class Annovar(pyCMMBase):
    """ A structure to parse and keep sample information """

    def __init__(self,
                 annovar_params,
                 ):
        self.__annovar_params = annovar_params

    @property
    def input_file(self):
        return self.__annovar_params[ANNOVAR_PARAMS_INPUT_FILE_KEY]

    @property
    def db_folder(self):
        return self.__annovar_params[ANNOVAR_PARAMS_DB_FOLDER_KEY]

    @property
    def buildver(self):
        return self.__annovar_params[ANNOVAR_PARAMS_BUILDVER_KEY]

    @property
    def out_prefix(self):
        return self.__annovar_params[ANNOVAR_PARAMS_OUT_PREFIX_KEY]

    @property
    def protocols(self):
        return self.__annovar_params[ANNOVAR_PARAMS_DB_NAMES_KEY]

    @property
    def operations(self):
        return self.__annovar_params[ANNOVAR_PARAMS_DB_OPS_KEY]

    @property
    def nastring(self):
        return self.__annovar_params[ANNOVAR_PARAMS_NASTRING_KEY]

    @property
    def annotated_vcf(self):
        return self.out_prefix + "." + self.buildver + "_multianno.vcf"

    @property
    def table_annovar_cmd(self):
        cmd = "table_annovar.pl"
        cmd += " " + self.input_file
        cmd += " " + self.db_folder
        cmd += " -buildver " + self.buildver
        cmd += " -out " + self.out_prefix
        cmd += " -remove"
        cmd += " -protocol " + self.protocols
        cmd += " -operation " + self.operations
        cmd += " -nastring " + self.nastring
        cmd += " -vcfinput"
        return cmd

_DESCRIPTION = "Description"
_HARMFUL = "Harmful"


class PredictionInfo(pyCMMBase):
    """
    To encapsulate the effect prediction information
    """

    def __init__(self,
                 code,
                 description,
                 harmful,
                 ):
        self.__code = code
        self.__description = description
        self.__harmful = harmful

    @property
    def code(self):
        return self.__code

    @property
    def description(self):
        if self.__code == ".":
            return self.__code
        if self.__description == ".":
            return "NA(" + self.__code + ")"
        return self.__description

    @property
    def harmful(self):
        return self.__harmful


class PredictionTranslator(pyCMMBase):
    """
    A class to translate codes from effect predictors by using informaiton
    from http://annovar.openbioinformatics.org/en/latest/user-guide/filter/
    As 20 August 2015 the translation is as below

    SIFT (sift)                        D: Deleterious (sift<=0.05); T: tolerated (sift>0.tolerated05)
    PolyPhen 2 HDIV (pp2_hdiv)         D: Probably damaging (>=0.957),damaging P: possibly damaging (0.453<=pp2_hdiv<=0.956); B: benign (pp2_hdiv<=0.452)
    PolyPhen 2 HVar (pp2_hvar)         D: Probably damaging (>=0.909), P: possibly damaging (0.447<=pp2_hdiv<=0.909); B: benign (pp2_hdiv<=0.446)
    LRT (lrt)                          D: Deleterious; N: Neutral; U: Unknown
    MutationTaster (mt)                A" ("disease_causing_automatic"); "D" ("disease_causing"); "N" ("polymorphism"); "P" ("polymorphism_automatic")
    MutationAssessor (ma)              H: high; M: medium; L: low; N: neutral. H/M means functional and L/N   means non-functional
    FATHMM (fathmm)                    D: Deleterious; T: Tolerated
    MetaSVM (metasvm)                  D: Deleterious; T: Tolerated
    MetaLR (metasvmalr)                D: Deleterious; T: Tolerated
    GERP++ (gerp++)                    higher scores are more deleterious
    PhyloP (phylop)                    higher scores are more deleterious
    SiPhy (siphy)                      higher scores are more deleterious
    """

    def __init__(self):
        self.__set_prediction_info()
        self.__null_prediction = PredictionInfo(code='.',
                                                description='.',
                                                harmful=False,
                                                )

    def __set_prediction_info(self):
        self.__pred_info = {}

        self.__pred_info[SIFT_PRED_COL] = {}
        pred_info = PredictionInfo(code='D',
                                   description='Deleterious',
                                   harmful=True,
                                   )
        self.__pred_info[SIFT_PRED_COL]['D'] = pred_info
        pred_info = PredictionInfo(code='T',
                                   description='Tolerated',
                                   harmful=False,
                                   )
        self.__pred_info[SIFT_PRED_COL]['T'] = pred_info

        self.__pred_info[POLYPHEN2_HDIV_PRED_COL] = {}
        pred_info = PredictionInfo(code='D',
                                   description='Probably damaging',
                                   harmful=True,
                                   )
        self.__pred_info[POLYPHEN2_HDIV_PRED_COL]['D'] = pred_info
        pred_info = PredictionInfo(code='P',
                                   description='Possibly damaging',
                                   harmful=True,
                                   )
        self.__pred_info[POLYPHEN2_HDIV_PRED_COL]['P'] = pred_info
        pred_info = PredictionInfo(code='B',
                                   description='Benign',
                                   harmful=False,
                                   )
        self.__pred_info[POLYPHEN2_HDIV_PRED_COL]['B'] = pred_info

        self.__pred_info[POLYPHEN2_HVAR_PRED_COL] = {}
        pred_info = PredictionInfo(code='D',
                                   description='Probably damaging',
                                   harmful=True,
                                   )
        self.__pred_info[POLYPHEN2_HVAR_PRED_COL]['D'] = pred_info
        pred_info = PredictionInfo(code='P',
                                   description='Possibly damaging',
                                   harmful=True,
                                   )
        self.__pred_info[POLYPHEN2_HVAR_PRED_COL]['P'] = pred_info
        pred_info = PredictionInfo(code='B',
                                   description='Benign',
                                   harmful=False,
                                   )
        self.__pred_info[POLYPHEN2_HVAR_PRED_COL]['B'] = pred_info

        self.__pred_info[LRT_PRED_COL] = {}
        pred_info = PredictionInfo(code='D',
                                   description='Deleterious',
                                   harmful=True,
                                   )
        self.__pred_info[LRT_PRED_COL]['D'] = pred_info
        pred_info = PredictionInfo(code='N',
                                   description='Neutral',
                                   harmful=False,
                                   )
        self.__pred_info[LRT_PRED_COL]['N'] = pred_info
        pred_info = PredictionInfo(code='U',
                                   description='Unknown',
                                   harmful=False,
                                   )
        self.__pred_info[LRT_PRED_COL]['U'] = pred_info

        self.__pred_info[MUTATIONTASTER_PRED_COL] = {}
        pred_info = PredictionInfo(code='A',
                                   description='Disease causing automatic',
                                   harmful=True,
                                   )
        self.__pred_info[MUTATIONTASTER_PRED_COL]['A'] = pred_info
        pred_info = PredictionInfo(code='D',
                                   description='Disease causing',
                                   harmful=True,
                                   )
        self.__pred_info[MUTATIONTASTER_PRED_COL]['D'] = pred_info
        pred_info = PredictionInfo(code='N',
                                   description='Polymorphism',
                                   harmful=False,
                                   )
        self.__pred_info[MUTATIONTASTER_PRED_COL]['N'] = pred_info
        pred_info = PredictionInfo(code='P',
                                   description='Polymorphism automatic',
                                   harmful=False,
                                   )
        self.__pred_info[MUTATIONTASTER_PRED_COL]['P'] = pred_info

        self.__pred_info[MUTATIONASSESSOR_PRED_COL] = {}
        pred_info = PredictionInfo(code='H',
                                   description='High',
                                   harmful=True,
                                   )
        self.__pred_info[MUTATIONASSESSOR_PRED_COL]['H'] = pred_info
        pred_info = PredictionInfo(code='M',
                                   description='Medium',
                                   harmful=True,
                                   )
        self.__pred_info[MUTATIONASSESSOR_PRED_COL]['M'] = pred_info
        pred_info = PredictionInfo(code='L',
                                   description='Low',
                                   harmful=False,
                                   )
        self.__pred_info[MUTATIONASSESSOR_PRED_COL]['L'] = pred_info
        pred_info = PredictionInfo(code='N',
                                   description='Neutral',
                                   harmful=False,
                                   )
        self.__pred_info[MUTATIONASSESSOR_PRED_COL]['N'] = pred_info
        pred_info = PredictionInfo(code='H/M',
                                   description='Functional',
                                   harmful=True,
                                   )
        self.__pred_info[MUTATIONASSESSOR_PRED_COL]['H/M'] = pred_info
        pred_info = PredictionInfo(code='L/N',
                                   description='Non-functional',
                                   harmful=False,
                                   )
        self.__pred_info[MUTATIONASSESSOR_PRED_COL]['L/N'] = pred_info

        self.__pred_info[FATHMM_PRED_COL] = {}
        pred_info = PredictionInfo(code='D',
                                   description='Deleterious',
                                   harmful=True,
                                   )
        self.__pred_info[FATHMM_PRED_COL]['D'] = pred_info
        pred_info = PredictionInfo(code='T',
                                   description='Tolerated',
                                   harmful=False,
                                   )
        self.__pred_info[FATHMM_PRED_COL]['T'] = pred_info

        self.__pred_info[RADIALSVM_PRED_COL] = {}
        pred_info = PredictionInfo(code='D',
                                   description='Deleterious',
                                   harmful=True,
                                   )
        self.__pred_info[RADIALSVM_PRED_COL]['D'] = pred_info
        pred_info = PredictionInfo(code='T',
                                   description='Tolerated',
                                   harmful=False,
                                   )
        self.__pred_info[RADIALSVM_PRED_COL]['T'] = pred_info

        self.__pred_info[LR_PRED_COL] = {}
        pred_info = PredictionInfo(code='D',
                                   description='Deleterious',
                                   harmful=True,
                                   )
        self.__pred_info[LR_PRED_COL]['D'] = pred_info
        pred_info = PredictionInfo(code='T',
                                   description='Tolerated',
                                   harmful=False,
                                   )
        self.__pred_info[LR_PRED_COL]['T'] = pred_info

    def get_prediction_info(self,
                            predictor_name,
                            code,
                            ):
        # if no prediction
        if code == ".":
            return self.__null_prediction
        # if predictor is not yet in the system
        if predictor_name not in self.__pred_info.keys():
            warning_msg = "!! Unknown predictor '" + predictor_name + "'"
            self.warning(warning_msg)
            # create a dummy for the new predictor entity in the system
            self.__pred_info[predictor_name] = {}
            pred_info = PredictionInfo(code=code,
                                       description='.',
                                       harmful=False,
                                       )
            self.__pred_info[predictor_name][code] = pred_info
        # if prediction code is not yet in the system
        elif code not in self.__pred_info[predictor_name].keys():
            warning_msg = "!! Unknown prediction code '" + code + "'"
            warning_msg += " for predictor '" + predictor_name + "'"
            self.warning(warning_msg)
            # create a dummy for the new prediction code entity in the system
            pred_info = PredictionInfo(code=code,
                                       description='.',
                                       harmful=False,
                                       )
            self.__pred_info[predictor_name][code] = pred_info
        return self.__pred_info[predictor_name][code]

    @property
    def predictor_list(self):
        return self.__pred_info.keys()
